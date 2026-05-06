// src/modules/io.rs
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::collections::HashSet;
use anyhow::{Context, Result};
use ndarray::Array2;
use regex::Regex;

use super::types::SnpInfo;

pub fn read_map_file(path: &str) -> Result<Vec<SnpInfo>> {
    let file = File::open(path).context("Failed to open map file")?;
    let reader = BufReader::new(file);
    let mut snps = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i == 0 && (line.starts_with("SNPID") || line.starts_with("SNP")) {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            snps.push(SnpInfo {
                snpid: parts[0].to_string(),
                chr: parts[1].parse().context("Invalid chromosome number")?,
                position: parts[2].parse().context("Invalid position")?,
            });
        }
    }

    Ok(snps)
}

pub fn expand_and_sort_files(patterns: &[String]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    
    for pattern in patterns {
        if pattern.contains('{') {
            let re = Regex::new(r"(.*)\{(\d+)\.\.(\d+)\}(.*)").unwrap();
            if let Some(caps) = re.captures(pattern) {
                let prefix = &caps[1];
                let start: i32 = caps[2].parse()?;
                let end: i32 = caps[3].parse()?;
                let suffix = &caps[4];
                
                for i in start..=end {
                    let path = format!("{}{}{}", prefix, i, suffix);
                    if Path::new(&path).exists() {
                        files.push(PathBuf::from(path));
                    }
                }
            }
        } else {
            if Path::new(pattern).exists() {
                files.push(PathBuf::from(pattern));
            }
        }
    }

    files.sort_by_key(|f| get_chromosome_number(f));
    Ok(files)
}

pub fn get_chromosome_number(path: &Path) -> i32 {
    let filename = path.file_name().unwrap().to_str().unwrap();
    let re = Regex::new(r"chr(\d+)").unwrap();
    
    if let Some(caps) = re.captures(filename) {
        return caps[1].parse().unwrap_or(0);
    }
    
    let nums: String = filename.chars().filter(|c| c.is_numeric()).collect();
    nums.parse().unwrap_or(0)
}

pub fn read_haplotype_file(path: &Path, noheader: bool, verbose: bool) -> Result<Array2<f64>> {
    if verbose {
        println!("    Reading: {}", path.display());
    }

    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<Vec<_>, _>>()?;

    if lines.is_empty() {
        anyhow::bail!("File is empty");
    }

    let start_idx = if noheader { 0 } else { 1 };
    
    if start_idx >= lines.len() {
        anyhow::bail!("No data after header");
    }

    let first_parts: Vec<&str> = lines[start_idx].split_whitespace().collect();
    let n_hap_cols = first_parts.len() - 1;
    let n_snps = n_hap_cols / 2;

    if verbose {
        println!("    Total haplotype columns: {}", n_hap_cols);
        println!("    Inferred n_snps: {}", n_snps);
        println!("    First 10 values: {:?}", &first_parts[1..11.min(first_parts.len())]);
    }

    let mut all_hap1 = Vec::new();
    let mut all_hap2 = Vec::new();

    for (line_num, line) in lines.iter().skip(start_idx).enumerate() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        
        if parts.len() < 2 {
            continue;
        }

        let alleles = &parts[1..];
        
        if alleles.len() != n_hap_cols {
            if verbose {
                println!("    Warning: Line {} has {} cols, expected {}", 
                         line_num + start_idx + 1, alleles.len(), n_hap_cols);
            }
            continue;
        }

        let mut hap1 = Vec::new();
        let mut hap2 = Vec::new();

        for i in (0..alleles.len()).step_by(2) {
            if i + 1 < alleles.len() {
                let a1 = alleles[i].parse::<i32>().ok();
                let a2 = alleles[i + 1].parse::<i32>().ok();

                match (a1, a2) {
                    (Some(v1), Some(v2)) => {
                        hap1.push(if v1 == 1 { 0.0 } else { 1.0 });
                        hap2.push(if v2 == 1 { 0.0 } else { 1.0 });
                    }
                    _ => {
                        hap1.push(f64::NAN);
                        hap2.push(f64::NAN);
                    }
                }
            }
        }

        all_hap1.push(hap1);
        all_hap2.push(hap2);
    }

    if all_hap1.is_empty() {
        anyhow::bail!("No valid haplotype data found");
    }

    let mut all_haplotypes = Vec::new();
    for (h1, h2) in all_hap1.iter().zip(all_hap2.iter()) {
        all_haplotypes.push(h1.clone());
        all_haplotypes.push(h2.clone());
    }

    let n_haplotypes = all_haplotypes.len();
    let hap_array = Array2::from_shape_vec(
        (n_haplotypes, n_snps),
        all_haplotypes.into_iter().flatten().collect(),
    )?;

    if verbose {
        println!("    Final shape: {} × {} (n_haplotypes × n_snps)", 
                 hap_array.nrows(), hap_array.ncols());
        
        let unique_vals: HashSet<i32> = hap_array
            .iter()
            .filter(|v| !v.is_nan())
            .map(|v| *v as i32)
            .collect();
        println!("    Unique values: {:?}", unique_vals);
    }

    Ok(hap_array)
}