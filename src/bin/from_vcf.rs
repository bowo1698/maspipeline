// src/from_vcf.rs
use clap::Parser;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::{File, create_dir_all};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use anyhow::{Result, Context};

#[derive(Parser)]
#[command(name = "convert-from-vcf")]
#[command(about = "Convert VCF to genotype/haplotype format")]
struct Args {
    #[arg(short = 'i', long, required = true, num_args = 1..)]
    input: Vec<PathBuf>,
    
    #[arg(short = 'm', long, default_value = "-9999")]
    missing: String,
    
    #[arg(long, default_value = "map_new.txt")]
    map: PathBuf,
    
    #[arg(long, default_value = "geno")]
    genofolder: PathBuf,
    
    #[arg(long, default_value = "hap")]
    hapfolder: PathBuf,
    
    #[arg(long, default_value = "chr")]
    genoprefix: String,
    
    #[arg(long, default_value = "chr")]
    happrefix: String,
    
    #[arg(long)]
    nosort: bool,
    
    #[arg(long, default_value = "1")]
    interval: usize,
    
    #[arg(short = 'V', long)]
    verbose: bool,
}

fn sort_files(files: &mut Vec<PathBuf>, nosort: bool) {
    if !nosort {
        let re = Regex::new(r"\d+").unwrap();
        files.sort_by_key(|f| {
            f.to_str()
                .and_then(|s| re.find(s))
                .and_then(|m| m.as_str().parse::<u32>().ok())
                .unwrap_or(u32::MAX)
        });
    }
}

fn decode_snp(s: &str, missing: &str) -> String {
    let alleles: Vec<&str> = s.split(|c| c == '|' || c == '/').collect();
    
    if alleles.len() == 2 {
        if alleles.contains(&".") {
            missing.to_string()
        } else if alleles[0] == alleles[1] {
            match alleles[0] {
                "1" => "2".to_string(),
                "0" => "0".to_string(),
                _ => {
                    eprintln!("ERROR: {}, setting as missing", s);
                    missing.to_string()
                }
            }
        } else {
            "1".to_string()
        }
    } else {
        eprintln!("ERROR: cannot split {}, setting as missing", s);
        missing.to_string()
    }
}

fn decode_hap(s: &str, missing: &str) -> String {
    let alleles: Vec<&str> = s.split('|').collect();
    
    if alleles.len() == 2 {
        let mut out = Vec::new();
        for a in alleles {
            match a {
                "0" => out.push("1"),
                "1" => out.push("2"),
                _ => out.push(missing),
            }
        }
        out.join("\t")
    } else {
        if s.contains('/') {
            eprintln!("ERROR: {}, found slash separator - is VCF phased?", s);
        }
        format!("{}\t{}", missing, missing)
    }
}

fn main() -> Result<()> {
    let mut args = Args::parse();
    
    eprintln!("VCF to Genotype/Haplotype Converter\n");
    
    create_dir_all(&args.genofolder)?;
    create_dir_all(&args.hapfolder)?;
    
    sort_files(&mut args.input, args.nosort);
    
    let mut map_data: HashMap<String, HashMap<u32, String>> = HashMap::new();
    
    for vcf_file in &args.input {
        if args.verbose {
            eprintln!("Processing: {}", vcf_file.display());
        }
        
        let file = File::open(vcf_file)
            .with_context(|| format!("Cannot open {}", vcf_file.display()))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        
        // Skip header lines
        let mut individuals = Vec::new();

        for line in lines.by_ref() {
            let line = line?;
            if line.starts_with("##") {
                continue;
            } else if line.starts_with("#CHROM") {
                let parts: Vec<&str> = line.split('\t').collect();
                individuals = parts[9..].iter().map(|s| s.to_string()).collect();
                break;
            }
        }
        
        let num_ind = individuals.len();
        
        // Read map data
        let file = File::open(vcf_file)?;
        let reader = BufReader::new(file);
        let mut map_data_local = Vec::new();
        
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 3 {
                let snp_id = parts[2].to_string();
                let chr = parts[0].to_string();
                let pos: u32 = parts[1].parse().unwrap_or(0);
                
                map_data_local.push((snp_id.clone(), chr.clone(), pos));
                map_data.entry(chr.clone())
                    .or_insert_with(HashMap::new)
                    .insert(pos, snp_id);
            }
        }
        
        let chromosomes: HashSet<String> = map_data_local.iter()
            .map(|(_, chr, _)| chr.clone())
            .collect();
        
        let interval = (num_ind / args.interval).max(1);
        
        if args.verbose {
            eprintln!("  Individuals: {}", num_ind);
            eprintln!("  Chromosomes: {}", chromosomes.len());
            eprintln!("  Interval: {} individuals\n", interval);
        }
        
        // Open output files
        let mut hap_files: HashMap<String, BufWriter<File>> = HashMap::new();
        let mut geno_files: HashMap<String, BufWriter<File>> = HashMap::new();
        
        for chr in &chromosomes {
            let hap_path = args.hapfolder.join(format!("{}{}", args.happrefix, chr));
            let geno_path = args.genofolder.join(format!("{}{}", args.genoprefix, chr));
            
            if args.verbose {
                eprintln!("  Opening {} for writing", hap_path.display());
                eprintln!("  Opening {} for writing", geno_path.display());
            }
            
            let hap_file = BufWriter::new(File::create(&hap_path)?);
            hap_files.insert(chr.clone(), hap_file);
            
            let mut geno_file = BufWriter::new(File::create(&geno_path)?);
            
            // Write geno header
            let snp_header: Vec<String> = map_data_local.iter()
                .filter(|(_, c, _)| c == chr)
                .map(|(snp, _, _)| snp.clone())
                .collect();
            writeln!(geno_file, "ID\t{}", snp_header.join("\t"))?;
            
            geno_files.insert(chr.clone(), geno_file);
        }
        
        if args.verbose {
            eprintln!("\nConverting...");
        }
        
        // Process in chunks
        let mut pass_num = 0;
        for chunk_start in (9..num_ind + 9).step_by(interval) {
            let chunk_end = (chunk_start + interval).min(num_ind + 9);
            if chunk_start >= chunk_end {
                break;
            }
            
            pass_num += 1;
            if args.verbose {
                eprintln!("  Pass {}: individuals {}-{}", 
                         pass_num, 
                         chunk_start - 9, 
                         chunk_end - 9);
            }
            
            let file = File::open(vcf_file)?;
            let reader = BufReader::new(file);
            
            let ind_start = chunk_start - 9;
            let ind_end = chunk_end - 9;
            let chunk_inds = &individuals[ind_start..ind_end];
            
            let mut chr_data: HashMap<String, Vec<Vec<String>>> = HashMap::new();
            for chr in &chromosomes {
                chr_data.insert(chr.clone(), vec![chunk_inds.to_vec()]);
            }
            
            for line in reader.lines() {
                let line = line?;
                if line.starts_with('#') {
                    continue;
                }
                
                let parts: Vec<&str> = line.split('\t').collect();
                let chr = parts[0];
                let genotypes: Vec<String> = parts[chunk_start..chunk_end]
                    .iter()
                    .map(|s| s.to_string())
                    .collect();
                
                if let Some(data) = chr_data.get_mut(chr) {
                    data.push(genotypes);
                }
            }
            
            // Write output
            for chr in &chromosomes {
                if let Some(data) = chr_data.get(chr) {
                    if data.len() <= 1 {
                        continue;
                    }
                    
                    let transposed = transpose(data);
                    
                    let hap_file = hap_files.get_mut(chr).unwrap();
                    let geno_file = geno_files.get_mut(chr).unwrap();
                    
                    for row in transposed {
                        let ind_id = &row[0];
                        
                        let hap_vals: Vec<String> = row[1..].iter()
                            .map(|g| decode_hap(g, &args.missing))
                            .collect();
                        writeln!(hap_file, "{}\t{}", ind_id, hap_vals.join("\t"))?;
                        
                        let geno_vals: Vec<String> = row[1..].iter()
                            .map(|g| decode_snp(g, &args.missing))
                            .collect();
                        writeln!(geno_file, "{}\t{}", ind_id, geno_vals.join("\t"))?;
                    }
                }
            }
        }
        
        if args.verbose {
            eprintln!("  Done processing {}\n", vcf_file.display());
        }
    }
    
    // Write map file
    if args.verbose {
        eprintln!("Writing map file: {}", args.map.display());
    }
    
    let mut map_file = BufWriter::new(File::create(&args.map)?);
    writeln!(map_file, "SNPID\tChr\tPosition")?;
    
    let mut sorted_chrs: Vec<_> = map_data.keys().collect();
    sorted_chrs.sort_by_key(|c| c.parse::<u32>().unwrap_or(u32::MAX));
    
    for chr in sorted_chrs {
        let mut sorted_pos: Vec<_> = map_data[chr].keys().collect();
        sorted_pos.sort();
        
        for pos in sorted_pos {
            writeln!(map_file, "{}\t{}\t{}", map_data[chr][pos], chr, pos)?;
        }
    }
    
    eprintln!("\n✓ Conversion completed successfully!");
    
    Ok(())
}

fn transpose(data: &[Vec<String>]) -> Vec<Vec<String>> {
    if data.is_empty() {
        return Vec::new();
    }
    
    let rows = data.len();
    let cols = data[0].len();
    
    (0..cols)
        .map(|col| (0..rows).map(|row| data[row][col].clone()).collect())
        .collect()
}