// src/main.rs
use clap::Parser;
use csv::ReaderBuilder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use anyhow::{Result, Context};
use chrono::Local;

#[derive(Parser)]
#[command(name = "convert-to-vcf")]
#[command(about = "Convert CSV genotype to VCF format")]
struct Args {
    #[arg(short = 'i', long)]
    input: PathBuf,
    
    #[arg(short = 'm', long)]
    map: PathBuf,
    
    #[arg(short = 'o', long, default_value = "output.vcf")]
    output: PathBuf,
    
    #[arg(long, default_value = "T")]
    ref_allele: String,
    
    #[arg(long, default_value = "C")]
    alt: String,
    
    #[arg(long, default_value = ".")]
    qual: String,
    
    #[arg(long, default_value = "PASS")]
    filter: String,
    
    #[arg(long, default_value = ".")]
    info: String,
    
    #[arg(long, default_value = "GT")]
    format_field: String,
    
    #[arg(long)]
    split_by_chr: bool,
    
    #[arg(long)]
    parallel: bool,
    
    #[arg(long, default_value = "1")]
    ncores: usize,
    
    #[arg(long)]
    merge_after_split: bool,
    
    #[arg(long)]
    keep_chr_files: bool,
}

#[derive(Debug)]
struct SnpInfo {
    chr: String,
    pos: String,
    id: String,
}

fn read_map(path: &PathBuf) -> Result<HashMap<String, SnpInfo>> {
    let mut map = HashMap::new();
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;
    
    for result in rdr.records() {
        let record = result?;
        if record.len() >= 3 {
            let id = record[0].to_string();
            map.insert(id.clone(), SnpInfo {
                chr: record[1].to_string(),
                pos: record[2].to_string(),
                id,
            });
        }
    }
    Ok(map)
}

fn encode_gt(val: &str) -> &'static str {
    match val {
        "0" => "0/0",
        "1" => "0/1",
        "2" => "1/1",
        _ => "./.",
    }
}

fn write_vcf_header<W: Write>(writer: &mut W, samples: &[String]) -> Result<()> {
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##fileDate={}", Local::now().format("%Y%m%d"))?;
    writeln!(writer, "##source=AlphaSimR-to-VCF-converter-rust")?;
    writeln!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    write!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample in samples {
        write!(writer, "\t{}", sample)?;
    }
    writeln!(writer)?;
    Ok(())
}

fn process_chromosome(
    chr: &str,
    input: &PathBuf,
    map: &HashMap<String, SnpInfo>,
    args: &Args,
    output: &PathBuf,
) -> Result<usize> {
    let mut rdr = ReaderBuilder::new().from_path(input)?;
    let headers = rdr.headers()?.clone();
    
    let meta_cols = ["ID", "individual", "generation", "population_type", "individual_id"];
    let snp_indices: Vec<_> = headers.iter().enumerate()
        .filter(|(_, h)| !meta_cols.contains(&h))
        .filter(|(_, h)| map.get(*h).map_or(false, |s| s.chr == chr))
        .collect();

    if snp_indices.is_empty() {
        return Ok(0);
    }
    
    let snp_cols: Vec<_> = snp_indices.iter().map(|(_, h)| h.to_string()).collect();
    
    let mut samples = Vec::new();
    let mut genotypes = Vec::new();
    
    for result in rdr.records() {
        let record = result?;
        samples.push(record[0].to_string());
        
        let mut gts = Vec::new();
        for (idx, _) in &snp_indices {
            gts.push(record[*idx].to_string());
        }
        genotypes.push(gts);
    }
    
    let mut file = BufWriter::new(File::create(output)?);
    write_vcf_header(&mut file, &samples)?;
    
    let mut count = 0;
    for (snp_idx, snp_id) in snp_cols.iter().enumerate() {
        if let Some(info) = map.get(snp_id) {
            // VCF spec requires position >= 1; auto-correct position=0 to position=1
            let adjusted_pos = if info.pos == "0" { "1" } else { &info.pos };
            
            write!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                info.chr, adjusted_pos, info.id,
                args.ref_allele, args.alt, args.qual,
                args.filter, args.info, args.format_field)?;
            
            for ind_gts in &genotypes {
                write!(file, "\t{}", encode_gt(&ind_gts[snp_idx]))?;
            }
            writeln!(file)?;
            count += 1;
        }
    }
    
    Ok(count)
}

fn merge_vcfs(chr_files: &[(String, PathBuf)], output: &PathBuf, keep: bool) -> Result<()> {
    let mut out = BufWriter::new(File::create(output)?);
    let mut header_written = false;
    
    for (chr, path) in chr_files {
        let content = std::fs::read_to_string(path)?;
        let mut var_count = 0;
        
        for line in content.lines() {
            if line.starts_with("##") {
                if !header_written {
                    writeln!(out, "{}", line)?;
                }
            } else if line.starts_with("#CHROM") {
                if !header_written {
                    writeln!(out, "{}", line)?;
                    header_written = true;
                }
            } else {
                writeln!(out, "{}", line)?;
                var_count += 1;
            }
        }
        
        eprintln!("Merged chr{}: {} variants", chr, var_count);
        
        if !keep {
            std::fs::remove_file(path)?;
        }
    }
    
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    
    let map = read_map(&args.map).context("Failed to read map file")?;
    
    if args.split_by_chr {
        let mut chromosomes: Vec<String> = map.values()
            .map(|s| s.chr.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        
        chromosomes.sort_by_key(|c| c.parse::<u32>().unwrap_or(u32::MAX));
        
        if args.parallel && chromosomes.len() > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.ncores)
                .build_global()?;
            
            let chr_files: Vec<_> = chromosomes.par_iter()
                .map(|chr| {
                    let out_path = args.output.with_file_name(
                        format!("{}_chr{}.vcf", 
                                args.output.file_stem().unwrap().to_str().unwrap(),
                                chr)
                    );
                    let count = process_chromosome(chr, &args.input, &map, &args, &out_path)
                        .unwrap_or(0);
                    eprintln!("Chr{}: {} variants", chr, count);
                    (chr.clone(), out_path)
                })
                .collect();
            
            if args.merge_after_split {
                merge_vcfs(&chr_files, &args.output, args.keep_chr_files)?;
                eprintln!("Merged {} chromosomes -> {}", chr_files.len(), args.output.display());
            }
        } else {
            for chr in chromosomes {
                let out_path = args.output.with_file_name(
                    format!("{}_chr{}.vcf",
                            args.output.file_stem().unwrap().to_str().unwrap(),
                            chr)
                );
                let count = process_chromosome(&chr, &args.input, &map, &args, &out_path)?;
                eprintln!("Chr{}: {} variants", chr, count);
            }
        }
    } else {
        let mut rdr = ReaderBuilder::new().from_path(&args.input)?;
        let headers = rdr.headers()?.clone();
        
        let meta_cols = ["ID", "individual", "generation", "population_type", "individual_id"];
        let snp_cols: Vec<_> = headers.iter()
            .filter(|h| !meta_cols.contains(h))
            .map(|h| h.to_string())
            .collect();
        
        let mut samples = Vec::new();
        let mut genotypes = Vec::new();
        
        for result in rdr.records() {
            let record = result?;
            samples.push(record[0].to_string());
            
            let gts: Vec<_> = snp_cols.iter()
                .enumerate()
                .map(|(i, _)| record.get(i + 1).unwrap_or(".").to_string())
                .collect();
            genotypes.push(gts);
        }
        
        let mut file = BufWriter::new(File::create(&args.output)?);
        write_vcf_header(&mut file, &samples)?;
        
        for (snp_idx, snp_id) in snp_cols.iter().enumerate() {
            if let Some(info) = map.get(snp_id) {
                // VCF spec requires position >= 1; auto-correct position=0 to position=1
                let adjusted_pos = if info.pos == "0" { "1" } else { &info.pos };
                
                write!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    info.chr, adjusted_pos, info.id,
                    args.ref_allele, args.alt, args.qual,
                    args.filter, args.info, args.format_field)?;
                
                for ind_gts in &genotypes {
                    write!(file, "\t{}", encode_gt(&ind_gts[snp_idx]))?;
                }
                writeln!(file)?;
            }
        }
    }
    
    Ok(())
}