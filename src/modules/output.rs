// src/modules/output.rs
use std::collections::{HashMap, BTreeMap, HashSet};
use std::fs::{self, File};
use std::io::{Write, BufRead};
use std::path::{Path, PathBuf};
use anyhow::Result;
use rayon::prelude::*;

use super::types::*;
use super::io::*;
use super::scoring::calculate_allele_metrics;

pub fn write_output_files(
    all_blocks: &BTreeMap<i32, Vec<Block>>,
    output_dir: &str,
    verbose: bool,
) -> Result<usize> {
    fs::create_dir_all(output_dir)?;

    let mut block_counter = 1;
    let mut total_blocks = 0;

    for (chrom, blocks) in all_blocks {
        if blocks.is_empty() {
            continue;
        }

        let hap_block_file = format!("{}/hap_block_{}", output_dir, chrom);
        let hap_info_file = format!("{}/hap_block_info_{}", output_dir, chrom);

        let mut hb = File::create(&hap_block_file)?;
        let mut hi = File::create(&hap_info_file)?;

        for block in blocks {
            let block_id = format!("blk{}", block_counter);
            writeln!(hb, "{}\t {} {}", block_id, block.start_idx, block.end_idx)?;

            let indices: Vec<String> = (block.start_idx..=block.end_idx)
                .map(|i| i.to_string())
                .collect();
            writeln!(hi, "{}\t {}", block_id, indices.join("\t"))?;

            block_counter += 1;
        }

        total_blocks += blocks.len();

        if verbose {
            println!("  Chr {}: {} blocks → {}", chrom, blocks.len(), 
                     Path::new(&hap_block_file).file_name().unwrap().to_str().unwrap());
        }
    }

    Ok(total_blocks)
}

pub fn generate_haplotype_genotypes(
    all_blocks: &BTreeMap<i32, Vec<Block>>,
    hap_files: &[PathBuf],
    output_dir: &Path,
    missing_code: Option<&str>,
    missing_value: &str,
    noheader: bool,
    verbose: bool,
    allele_map_dir: Option<&Path>,
) -> Result<()> {
    fs::create_dir_all(output_dir)?;

    if verbose {
        println!("\nGenerating haplotype genotypes...");
    }

    let missing_code_owned = missing_code.map(|s| s.to_string());
    let missing_value_owned = missing_value.to_string();
    let output_dir_owned = output_dir.to_path_buf();
    let allele_map_dir_owned = allele_map_dir.map(|p| p.to_path_buf());

    let results: Vec<Result<()>> = hap_files
        .par_iter()
        .map(|hap_file| {
            let chr_num = get_chromosome_number(hap_file);

            if !all_blocks.contains_key(&chr_num) || all_blocks[&chr_num].is_empty() {
                if verbose {
                    println!("  Chr {}: No blocks defined, skipping", chr_num);
                }
                return Ok(());
            }

            let blocks = &all_blocks[&chr_num];

            if verbose {
                println!("  Chr {}: Processing {} blocks", chr_num, blocks.len());
            }

            let file = File::open(hap_file)?;
            let reader = std::io::BufReader::new(file);
            let lines: Vec<String> = reader
                .lines()
                .collect::<std::result::Result<Vec<_>, std::io::Error>>()?;

            let start_idx = if noheader { 0 } else { 1 };
            let haps: Vec<Vec<String>> = lines
                .iter()
                .skip(start_idx)
                .filter_map(|line: &String| {
                    let parts: Vec<String> = line
                        .split_whitespace()
                        .map(|s: &str| s.to_string())
                        .collect();
                    if parts.len() >= 2 { Some(parts) } else { None }
                })
                .collect();

            let mut codes = Vec::new();
            for block in blocks {
                let mut label = 1;
                let mut hap_dict = HashMap::new();

                for hh in &haps {
                    let h: String = hh[1..].join("");
                    let start_pos = block.start_idx * 2;
                    let end_pos = block.end_idx * 2 + 2;
                    let h_part = &h[start_pos..end_pos];

                    let parent1: String = h_part.chars().step_by(2).collect();
                    let parent2: String = h_part.chars().skip(1).step_by(2).collect();

                    for p in [parent1, parent2] {
                        if let Some(ref mc) = missing_code_owned {
                            if p.contains(mc.as_str()) {
                                hap_dict.insert(p, missing_value_owned.clone());
                                continue;
                            }
                        }
                        if !hap_dict.contains_key(&p) {
                            hap_dict.insert(p, label.to_string());
                            label += 1;
                        }
                    }
                }
                codes.push(hap_dict);
            }

            // Save allele map
            if let Some(ref map_dir) = allele_map_dir_owned {
                fs::create_dir_all(map_dir)?;
                let map_file = map_dir.join(format!("allele_map_chr{}.csv", chr_num));
                let mut mf = File::create(&map_file)?;
                writeln!(mf, "block_id,haplotype_string,allele_label")?;
                let mut global_block_counter = 1usize;
                for (c, _block) in blocks.iter().enumerate() {
                    let block_id = format!("blk{}", global_block_counter);
                    for (hap_str, label) in &codes[c] {
                        writeln!(mf, "{},{},{}", block_id, hap_str, label)?;
                    }
                    global_block_counter += 1;
                }
                if verbose {
                    println!("    Allele map saved: {}", map_file.display());
                }
            }

            let outfile = output_dir_owned.join(format!("hap_geno_{}", chr_num));
            let mut o = File::create(&outfile)?;

            let mut header = vec!["ID".to_string()];
            for b in 1..=blocks.len() {
                header.push(format!("hap_{}_{}\thap_{}_{}", chr_num, b, chr_num, b));
            }
            writeln!(o, "{}", header.join("\t"))?;

            for hh in &haps {
                let ind = &hh[0];
                let h: String = hh[1..].join("");
                let mut output = vec![ind.clone()];

                for (c, block) in blocks.iter().enumerate() {
                    let start_pos = block.start_idx * 2;
                    let end_pos = block.end_idx * 2 + 2;
                    let h_part = &h[start_pos..end_pos];

                    let parent1: String = h_part.chars().step_by(2).collect();
                    let parent2: String = h_part.chars().skip(1).step_by(2).collect();

                    for p in [parent1, parent2] {
                        output.push(codes[c][&p].clone());
                    }
                }
                writeln!(o, "{}", output.join("\t"))?;
            }

            if verbose {
                println!("    Written: {}", outfile.display());
            }

            Ok(())
        })
        .collect();

    // propagate error jika ada
    for result in results {
        result?;
    }

    Ok(())
}

pub fn write_csv_outputs(
    all_blocks: &BTreeMap<i32, Vec<Block>>,
    hap_files: &[PathBuf],
    snp_map: &[SnpInfo],
    output_dir: &str,
    method_name: &str,
    noheader: bool,
    verbose: bool,
) -> Result<()> {
    let stats_dir = format!("{}/stats", output_dir);
    fs::create_dir_all(&stats_dir)?;

    if verbose {
        println!("\nGenerating analysis CSV files...");
        println!("  Output directory: {}/", stats_dir);
    }

    // === 1. MICROHAPLOTYPE COORDINATES CSV ===
    let mut coords_data = Vec::new();
    let mut block_counter = 1;

    for (chr_num, blocks) in all_blocks {
        for block in blocks {
            let mut row = vec![
                format!("blk{}", block_counter),
                chr_num.to_string(),
                block.start_pos.to_string(),
                block.end_pos.to_string(),
                block.n_snps.to_string(),
                (block.end_pos - block.start_pos).to_string(),
            ];

            if let Some(ld) = block.mean_ld_r2 {
                row.push(ld.to_string());
            }
            if let Some(score) = block.criterion_b_score {
                row.push(score.to_string());
            }

            coords_data.push(row);
            block_counter += 1;
        }
    }

    let coords_file = format!("{}/microhaplotype_coordinates.csv", stats_dir);
    let mut f = File::create(&coords_file)?;
    writeln!(f, "block_id,chr,start_pos,end_pos,n_snps,physical_span_bp,pic")?;
    for (row, block) in coords_data.iter().zip(all_blocks.values().flatten()) {
        let pic_str = block.pic.map_or("NA".to_string(), |p| format!("{:.6}", p));
        writeln!(f, "{},{}", row.join(","), pic_str)?;
    }

    if verbose {
        println!("  ✓ Coordinates: {}", coords_file);
    }

    // === 2. ALLELE FREQUENCY & DIVERSITY METRICS ===
    let mut allele_data: Vec<AlleleMetricsRow> = Vec::new();
    let mut all_freqs = Vec::new(); // For frequency spectrum
    
    let mut global_block_id = 0;
    for hap_file in hap_files {
        let chr_num = get_chromosome_number(hap_file);
        
        if !all_blocks.contains_key(&chr_num) || all_blocks[&chr_num].is_empty() {
            continue;
        }
        
        let hap_matrix = match read_haplotype_file(hap_file, noheader, false) {
            Ok(mat) => mat,
            Err(_) => continue,
        };
        
        let blocks = &all_blocks[&chr_num];
        
        for block in blocks {
            if let Some(metrics) = calculate_allele_metrics(&hap_matrix, block.start_idx, block.end_idx) {
                global_block_id += 1;
                
                allele_data.push(AlleleMetricsRow {
                    block_id: format!("blk{}", global_block_id),
                    chr: chr_num,
                    metrics,
                });
                
                // Collect frequencies for spectrum
                let snp_indices: Vec<usize> = (block.start_idx..=block.end_idx).collect();
                let block_haps = hap_matrix.select(ndarray::Axis(1), &snp_indices);
                
                let mut unique_haps = Vec::new();
                for row in block_haps.rows() {
                    let has_nan = row.iter().any(|v| v.is_nan());
                    if !has_nan {
                        let hap_tuple: Vec<i32> = row.iter().map(|v| *v as i32).collect();
                        unique_haps.push(hap_tuple);
                    }
                }
                
                if !unique_haps.is_empty() {
                    let mut counts = HashMap::new();
                    for hap in &unique_haps {
                        *counts.entry(hap.clone()).or_insert(0) += 1;
                    }
                    
                    let total = unique_haps.len() as f64;
                    let freqs: Vec<f64> = counts.values().map(|&c| c as f64 / total).collect();
                    all_freqs.extend(freqs);
                }
            }
        }
    }
    
    // Write allele metrics CSV
    if !allele_data.is_empty() {
        let allele_file = format!("{}/allele_frequency_metrics.csv", stats_dir);
        let mut f = File::create(&allele_file)?;
        writeln!(f, "block_id,chr,n_alleles,n_predictable,min_freq,max_freq,mean_freq,freq_entropy,expected_heterozygosity,observed_heterozygosity,pic,effective_alleles,rare_alleles_prop")?;
        
        for row in &allele_data {
            writeln!(
                f,
                "{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                row.block_id,
                row.chr,
                row.metrics.n_alleles,
                row.metrics.n_predictable,
                row.metrics.min_freq,
                row.metrics.max_freq,
                row.metrics.mean_freq,
                row.metrics.freq_entropy,
                row.metrics.expected_heterozygosity,
                row.metrics.observed_heterozygosity,
                row.metrics.pic,
                row.metrics.effective_alleles,
                row.metrics.rare_alleles_prop
            )?;
        }
        
        if verbose {
            println!("  ✓ Allele metrics: {}", allele_file);
        }
    }

    // === 3. CHROMOSOME SUMMARY STATS ===
    let mut chr_summary = Vec::new();
    for (chr_num, blocks) in all_blocks {
        if blocks.is_empty() {
            continue;
        }

        let chr_snps: Vec<&SnpInfo> = snp_map.iter().filter(|s| s.chr == *chr_num).collect();
        let total_snps = chr_snps.len();
        let chr_length = chr_snps.iter().map(|s| s.position).max().unwrap() 
                       - chr_snps.iter().map(|s| s.position).min().unwrap();

        let covered_bp: i64 = blocks.iter().map(|b| b.end_pos - b.start_pos).sum();
        let coverage_pct = if chr_length > 0 {
            (covered_bp as f64 / chr_length as f64) * 100.0
        } else {
            0.0
        };

        chr_summary.push(vec![
            chr_num.to_string(),
            blocks.len().to_string(),
            total_snps.to_string(),
            chr_length.to_string(),
            covered_bp.to_string(),
            format!("{:.2}", coverage_pct),
        ]);
    }

    let chr_file = format!("{}/chromosome_summary.csv", stats_dir);
    let mut f = File::create(&chr_file)?;
    writeln!(f, "chr,n_blocks,total_snps,chr_length_bp,coverage_bp,coverage_pct")?;
    for row in chr_summary {
        writeln!(f, "{}", row.join(","))?;
    }

    if verbose {
        println!("  ✓ Chromosome summary: {}", chr_file);
    }

    // === 4. ALLELE FREQUENCY SPECTRUM ===
    if !all_freqs.is_empty() {
        let bins = vec![0.0, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 1.0];
        let bin_labels = vec![
            "0.00-0.01", "0.01-0.05", "0.05-0.10", "0.10-0.15",
            "0.15-0.20", "0.20-0.25", "0.25-0.30", "0.30-0.40",
            "0.40-0.50", "0.50-1.00"
        ];
        
        let mut hist = vec![0; bins.len() - 1];
        for &freq in &all_freqs {
            for i in 0..bins.len() - 1 {
                if freq >= bins[i] && freq < bins[i + 1] {
                    hist[i] += 1;
                    break;
                }
                // Handle edge case for exactly 1.0
                if i == bins.len() - 2 && freq == bins[i + 1] {
                    hist[i] += 1;
                    break;
                }
            }
        }
        
        let spectrum_file = format!("{}/allele_frequency_spectrum.csv", stats_dir);
        let mut f = File::create(&spectrum_file)?;
        writeln!(f, "freq_bin,method,count")?;
        
        for (i, label) in bin_labels.iter().enumerate() {
            writeln!(f, "{},{},{}", label, method_name, hist[i])?;
        }
        
        if verbose {
            println!("  ✓ Frequency spectrum: {}", spectrum_file);
        }
    }

    // === 5. OVERALL SUMMARY ===
    let total_blocks: usize = all_blocks.values().map(|v| v.len()).sum();
    let total_snps: usize = all_blocks.values()
        .flat_map(|v| v.iter())
        .map(|b| b.n_snps)
        .sum();
    let mean_snps: f64 = total_snps as f64 / total_blocks as f64;
    
    let physical_spans: Vec<i64> = all_blocks.values()
        .flat_map(|v| v.iter())
        .map(|b| b.end_pos - b.start_pos)
        .collect();
    let mean_physical_span = if !physical_spans.is_empty() {
        physical_spans.iter().sum::<i64>() as f64 / physical_spans.len() as f64
    } else {
        0.0
    };
    
    let mut summary_data = vec![
        ("method".to_string(), method_name.to_string()),
        ("total_chromosomes".to_string(), all_blocks.len().to_string()),
        ("total_blocks".to_string(), total_blocks.to_string()),
        ("total_snps_covered".to_string(), total_snps.to_string()),
        ("mean_snps_per_block".to_string(), format!("{:.2}", mean_snps)),
        ("mean_physical_span".to_string(), format!("{:.2}", mean_physical_span)),
    ];
    
    if !allele_data.is_empty() {
        let mean_alleles = allele_data.iter().map(|r| r.metrics.n_alleles as f64).sum::<f64>() / allele_data.len() as f64;
        let mean_expected_het = allele_data.iter().map(|r| r.metrics.expected_heterozygosity).sum::<f64>() / allele_data.len() as f64;
        let mean_pic = allele_data.iter().map(|r| r.metrics.pic).sum::<f64>() / allele_data.len() as f64;
        
        summary_data.push(("mean_alleles_per_block".to_string(), format!("{:.2}", mean_alleles)));
        summary_data.push(("mean_expected_het".to_string(), format!("{:.6}", mean_expected_het)));
        summary_data.push(("mean_pic".to_string(), format!("{:.6}", mean_pic)));
    }

    let summary_file = format!("{}/method_summary.csv", stats_dir);
    let mut f = File::create(&summary_file)?;
    writeln!(f, "metric,value")?;
    for (metric, value) in summary_data {
        writeln!(f, "{},{}", metric, value)?;
    }

    if verbose {
        println!("  ✓ Method summary: {}", summary_file);
    }

    if verbose {
        let n_files = if allele_data.is_empty() { 3 } else { 5 };
        println!("\n  Generated {} CSV files in {}/", n_files, stats_dir);
    }

    Ok(())
}

// Helper struct for allele data
struct AlleleMetricsRow {
    block_id: String,
    chr: i32,
    metrics: AlleleMetrics,
}

pub fn write_snp_selection_file(
    all_blocks: &BTreeMap<i32, Vec<Block>>,
    snp_map: &[SnpInfo],
    output_dir: &str,
    verbose: bool,
) -> Result<()> {
    let stats_dir = format!("{}/stats", output_dir);
    fs::create_dir_all(&stats_dir)?;

    let mut snp_data = Vec::new();
    let mut block_counter = 1;

    for (chr_num, blocks) in all_blocks {
        let chr_snps: Vec<&SnpInfo> = snp_map.iter().filter(|s| s.chr == *chr_num).collect();

        for block in blocks {
            let block_id = format!("blk{}", block_counter);

            for idx in block.start_idx..=block.end_idx {
                if idx < chr_snps.len() {
                    let snp_info = chr_snps[idx];
                    snp_data.push(vec![
                        block_id.clone(),
                        chr_num.to_string(),
                        idx.to_string(),
                        snp_info.snpid.clone(),
                        snp_info.position.to_string(),
                    ]);
                }
            }

            block_counter += 1;
        }
    }

    let snp_file = format!("{}/snp_selection_detailed.csv", stats_dir);
    let mut f = File::create(&snp_file)?;
    writeln!(f, "block_id,chr,snp_index,snp_id,position")?;
    for row in snp_data.iter() {
        writeln!(f, "{}", row.join(","))?;
    }

    if verbose {
        println!("  ✓ SNP selection: {}", snp_file);
        println!("    Total SNPs selected: {}", snp_data.len());
        let unique_snps: HashSet<_> = snp_data.iter().map(|r| &r[3]).collect();
        println!("    Unique SNPs: {}", unique_snps.len());
    }

    Ok(())
}

pub fn calculate_genome_coverage(
    all_blocks: &BTreeMap<i32, Vec<Block>>,
    snp_map: &[SnpInfo],
    verbose: bool,
) -> Result<CoverageStats> {
    let mut total_genome_length = 0i64;
    let mut covered_length = 0i64;
    let mut chromosomes = BTreeMap::new();

    for (chr_num, blocks) in all_blocks {
        let chr_snps: Vec<&SnpInfo> = snp_map.iter().filter(|s| s.chr == *chr_num).collect();
        
        if chr_snps.is_empty() {
            continue;
        }

        let chr_length = chr_snps.iter().map(|s| s.position).max().unwrap()
                       - chr_snps.iter().map(|s| s.position).min().unwrap();
        let covered: i64 = blocks.iter().map(|b| b.end_pos - b.start_pos).sum();

        let coverage_pct = if chr_length > 0 {
            (covered as f64 / chr_length as f64) * 100.0
        } else {
            0.0
        };

        chromosomes.insert(*chr_num, ChromosomeCoverage {
            chr: *chr_num,
            n_blocks: blocks.len(),
            length_bp: chr_length,
            covered_bp: covered,
            coverage_pct,
        });

        total_genome_length += chr_length;
        covered_length += covered;
    }

    let coverage_pct = if total_genome_length > 0 {
        (covered_length as f64 / total_genome_length as f64) * 100.0
    } else {
        0.0
    };

    let total_blocks: usize = all_blocks.values().map(|v| v.len()).sum();

    if verbose {
        println!("\n{}", "=".repeat(70));
        println!("GENOME COVERAGE ANALYSIS");
        println!("{}", "=".repeat(70));
        println!("Total genome length:     {} bp", total_genome_length);
        println!("Covered by blocks:       {} bp", covered_length);
        println!("Coverage:                {:.2}%", coverage_pct);
        println!("{}", "=".repeat(70));
    }

    Ok(CoverageStats {
        total_genome_length_bp: total_genome_length,
        covered_length_bp: covered_length,
        coverage_pct,
        total_blocks,
        chromosomes,
    })
}

pub fn save_coverage_stats(
    stats: &CoverageStats,
    output_dir: &str,
    method_label: &str,
    verbose: bool,
) -> Result<()> {
    let stats_dir = format!("{}/stats", output_dir);
    fs::create_dir_all(&stats_dir)?;

    let coverage_file = format!("{}/genome_coverage.csv", stats_dir);
    let mut f = File::create(&coverage_file)?;

    writeln!(f, "method,genome_length_mb,covered_length_kb,coverage_pct")?;
    writeln!(
        f,
        "{},{:.2},{:.2},{:.2}",
        method_label,
        stats.total_genome_length_bp as f64 / 1e6,
        stats.covered_length_bp as f64 / 1e3,
        stats.coverage_pct
    )?;

    if verbose {
        println!("  ✓ Coverage stats saved: {}", coverage_file);
    }

    Ok(())
}

pub fn write_criterion_b_summary(
    all_blocks: &BTreeMap<i32, Vec<Block>>,
    output_dir: &str,
    verbose: bool,
) -> Result<()> {
    let stats_dir = format!("{}/stats", output_dir);
    fs::create_dir_all(&stats_dir)?;

    let mut scores: Vec<f64> = all_blocks.values()
        .flat_map(|blocks| blocks.iter())
        .filter_map(|b| b.criterion_b_score)
        .filter(|s| s.is_finite())
        .collect();

    let n_infinite = all_blocks.values()
        .flat_map(|b| b.iter())
        .filter_map(|b| b.criterion_b_score)
        .filter(|s| s.is_infinite())
        .count();

    if scores.is_empty() {
        return Ok(());
    }

    scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = scores.len();
    let mean = scores.iter().sum::<f64>() / n as f64;
    let median = scores[n / 2];
    let p25  = scores[n / 4];
    let p75  = scores[3 * n / 4];
    let p90  = scores[(n as f64 * 0.90) as usize];

    let bins   = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, f64::INFINITY];
    let labels = ["0.00-0.05","0.05-0.10","0.10-0.15","0.15-0.20",
                  "0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40",
                  "0.40-0.50",">0.50"];
    let mut hist = vec![0usize; labels.len()];
    for &s in &scores {
        for i in 0..bins.len()-1 {
            if s >= bins[i] && s < bins[i+1] { hist[i] += 1; break; }
        }
    }

    let summary_file = format!("{}/criterion_b_summary.txt", stats_dir);
    let mut f = File::create(&summary_file)?;
    writeln!(f, "CRITERION-B SCORE SUMMARY")?;
    writeln!(f, "=========================")?;
    writeln!(f, "Total blocks (finite):  {}", n)?;
    writeln!(f, "Blocks with inf score:  {}", n_infinite)?;
    writeln!(f, "Min:    {:.6}", scores[0])?;
    writeln!(f, "P25:    {:.6}", p25)?;
    writeln!(f, "Median: {:.6}", median)?;
    writeln!(f, "Mean:   {:.6}", mean)?;
    writeln!(f, "P75:    {:.6}", p75)?;
    writeln!(f, "P90:    {:.6}", p90)?;
    writeln!(f, "Max:    {:.6}", scores[n-1])?;
    writeln!(f, "\nDISTRIBUTION")?;
    writeln!(f, "------------")?;
    for (label, count) in labels.iter().zip(hist.iter()) {
        let pct = *count as f64 / n as f64 * 100.0;
        let bar: String = "#".repeat((*count as f64 / n as f64 * 40.0) as usize);
        writeln!(f, "{:12}  {:5}  ({:5.1}%)  {}", label, count, pct, bar)?;
    }
    writeln!(f, "\nPER-CHROMOSOME MEAN")?;
    writeln!(f, "-------------------")?;
    for (chr, blocks) in all_blocks {
        let chr_scores: Vec<f64> = blocks.iter()
            .filter_map(|b| b.criterion_b_score)
            .filter(|s| s.is_finite())
            .collect();
        if !chr_scores.is_empty() {
            let chr_mean = chr_scores.iter().sum::<f64>() / chr_scores.len() as f64;
            writeln!(f, "Chr {:2}:  {:.6}  (n={})", chr, chr_mean, chr_scores.len())?;
        }
    }

    if verbose {
        println!("  ✓ Criterion-B summary: {}", summary_file);
    }
    Ok(())
}