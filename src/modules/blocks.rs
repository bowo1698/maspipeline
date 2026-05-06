// src/modules/blocks.rs
use std::collections::{HashSet, BTreeMap};
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::Result;
use rayon::prelude::*;

use super::types::*;
use super::io::*;
use super::ld::*;
use super::scoring::*;

pub fn define_microhaplotype_blocks(
    hap_files: &[PathBuf],
    snp_map: &[SnpInfo],
    config: &BlockDefinitionConfig,
) -> Result<BTreeMap<i32, Vec<Block>>> {

    let results: Vec<Result<(i32, Vec<Block>)>> = hap_files
        .par_iter()
        .map(|hap_file| {
            let chr_num = get_chromosome_number(hap_file);

            let chr_snps: Vec<SnpInfo> = snp_map
                .iter()
                .filter(|s| s.chr == chr_num)
                .cloned()
                .collect();

            if chr_snps.is_empty() {
                if config.verbose {
                    println!("  Warning: No SNPs found for chromosome {}", chr_num);
                }
                return Ok((chr_num, vec![]));
            }

            let mut hap_matrix = None;
            if config.min_ld.is_some() || matches!(config.method, Method::LdHaploblock) {
                match read_haplotype_file(hap_file, config.noheader, config.verbose) {
                    Ok(mat) => {
                        if config.verbose {
                            println!("  Loaded: {} haplotypes × {} SNPs", mat.nrows(), mat.ncols());
                        }
                        if mat.ncols() != chr_snps.len() {
                            if config.verbose {
                                println!("  WARNING: Dimension mismatch!");
                                println!("    Haplotype file has {} SNPs", mat.ncols());
                                println!("    Map file has {} SNPs for chr {}", chr_snps.len(), chr_num);
                            }
                        }
                        hap_matrix = Some(mat);
                    }
                    Err(e) => {
                        if config.verbose {
                            println!("  Error loading haplotypes: {}", e);
                        }
                    }
                }
            }

            if config.verbose {
                println!("  Map SNPs: {}", chr_snps.len());
            }

            let mut blocks = Vec::new();

            match config.method {
                Method::SnpCountSimple => {
                    let n_snps = chr_snps.len();
                    let sd = config.window_snps;
                    let rem = n_snps % sd;
                    let blks = (n_snps - rem) / sd;

                    for i in (0..blks * sd).step_by(sd) {
                        blocks.push(Block {
                            chr: chr_num,
                            start_pos: chr_snps[i].position,
                            end_pos: chr_snps[i + sd - 1].position,
                            start_idx: i,
                            end_idx: i + sd - 1,
                            n_snps: sd,
                            mean_ld_r2: None,
                            criterion_b_score: None,
                            selection_type: None,
                            original_block_size: None,
                            physical_span: None,
                            split_type: None,
                            rare_allele_count: None,
                            pic: None,
                        });
                    }

                    if rem != 0 {
                        let start_idx = n_snps - rem;
                        blocks.push(Block {
                            chr: chr_num,
                            start_pos: chr_snps[start_idx].position,
                            end_pos: chr_snps[n_snps - 1].position,
                            start_idx,
                            end_idx: n_snps - 1,
                            n_snps: rem,
                            mean_ld_r2: None,
                            criterion_b_score: None,
                            selection_type: None,
                            original_block_size: None,
                            physical_span: None,
                            split_type: None,
                            rare_allele_count: None,
                            pic: None,
                        });
                    }
                }

                Method::LdHaploblock => {
                    let hap_mat = if let Some(ref mat) = hap_matrix {
                        mat
                    } else {
                        hap_matrix = Some(read_haplotype_file(hap_file, config.noheader, config.verbose)?);
                        hap_matrix.as_ref().unwrap()
                    };

                    let snp_positions: Vec<i64> = chr_snps.iter().map(|s| s.position).collect();

                    if config.verbose {
                        println!("  STAGE 1: Building LD haploblocks...");
                    }

                    let haploblocks = build_ld_haploblocks(
                        hap_mat,
                        Some(&snp_positions),
                        config.d_prime_threshold,
                        config.min_snps,
                        config.verbose,
                    );

                    // determine effective max snps: adaptive or from config
                    let effective_max_snps = if config.adaptive_max_snps {
                        let (r2_thresh, adaptive) = estimate_adaptive_max_snps(
                            hap_mat,
                            config.min_snps,
                            config.verbose,
                        );
                        if config.verbose {
                            println!("  Adaptive max_snps: {} (r²_threshold: {:.3})", adaptive, r2_thresh);
                        }
                        adaptive
                    } else {
                        config.max_snps
                    };

                    if config.verbose {
                        println!("  Built {} LD haploblocks", haploblocks.len());
                        match config.haplotype_type {
                            HaplotypeType::Pure => {
                                println!("  STAGE 2: Selecting best {}-SNP haplotypes from each block...", config.window_snps);
                            }
                            HaplotypeType::Micro => {
                                println!("  STAGE 2: Splitting blocks into {}bp windows (microhaplotypes)...", config.window_bp);
                            }
                        }
                    }

                    let mut n_output = 0;
                    let mut n_filtered_b = 0;
                    for haploblock in &haploblocks {
                        match config.haplotype_type {
                            HaplotypeType::Pure => {
                                let selected = select_best_haplotype_from_haploblock(
                                    &haploblock,
                                    hap_mat,
                                    config.window_snps,
                                    config.aft,
                                    config.md,
                                );
                                let sel_indices = &selected.snp_indices;
                                blocks.push(Block {
                                    chr: chr_num,
                                    start_pos: chr_snps[sel_indices[0]].position,
                                    end_pos: chr_snps[*sel_indices.last().unwrap()].position,
                                    start_idx: sel_indices[0],
                                    end_idx: *sel_indices.last().unwrap(),
                                    n_snps: selected.n_snps,
                                    mean_ld_r2: None,
                                    criterion_b_score: Some(selected.criterion_b_score),
                                    selection_type: Some(selected.selection_type),
                                    original_block_size: Some(selected.original_block_size),
                                    physical_span: None,
                                    split_type: None,
                                    rare_allele_count: None,
                                    pic: None,
                                });
                                n_output += 1;
                            }
                            HaplotypeType::Micro => {
                                let microhaplotypes = split_ld_block_by_physical_window(
                                    &haploblock,
                                    hap_mat,
                                    &chr_snps,
                                    config.window_bp,
                                    config.min_snps,
                                    effective_max_snps,
                                    config.aft,
                                    config.md,
                                );
                                for micro in microhaplotypes {
                                    if let Some(max_b) = config.max_criterion_b {
                                        if micro.criterion_b_score > max_b {
                                            n_filtered_b += 1;
                                            continue;
                                        }
                                    }
                                    let sel_indices = &micro.snp_indices;
                                    blocks.push(Block {
                                        chr: chr_num,
                                        start_pos: chr_snps[sel_indices[0]].position,
                                        end_pos: chr_snps[*sel_indices.last().unwrap()].position,
                                        start_idx: sel_indices[0],
                                        end_idx: *sel_indices.last().unwrap(),
                                        n_snps: micro.n_snps,
                                        mean_ld_r2: None,
                                        criterion_b_score: Some(micro.criterion_b_score),
                                        selection_type: None,
                                        original_block_size: None,
                                        physical_span: Some(micro.physical_span),
                                        split_type: Some(micro.split_type),
                                        rare_allele_count: None,
                                        pic: None,
                                    });
                                    n_output += 1;
                                }
                            }
                        }
                    }

                    if config.verbose {
                        let output_label = match config.haplotype_type {
                            HaplotypeType::Pure => "haplotypes",
                            HaplotypeType::Micro => "microhaplotypes",
                        };
                        println!("  Generated {} {} from {} LD blocks", n_output, output_label, haploblocks.len());
                        if config.max_criterion_b.is_some() {
                            println!("  Criterion-B filter: {} removed, {} retained",
                                n_filtered_b, n_output);
                        }
                        let scores: Vec<f64> = blocks
                            .iter()
                            .filter_map(|b| b.criterion_b_score)
                            .filter(|s| !s.is_infinite())
                            .collect();
                        if !scores.is_empty() {
                            let mean_score = scores.iter().sum::<f64>() / scores.len() as f64;
                            println!("  Mean Criterion-B score: {:.6}", mean_score);
                        }
                        match config.haplotype_type {
                            HaplotypeType::Pure => {
                                let mut type_counts = std::collections::HashMap::new();
                                for block in &blocks {
                                    if let Some(ref sel_type) = block.selection_type {
                                        *type_counts.entry(sel_type.clone()).or_insert(0) += 1;
                                    }
                                }
                                println!("  Selection types: {:?}", type_counts);
                            }
                            HaplotypeType::Micro => {
                                let spans: Vec<i64> = blocks.iter().filter_map(|b| b.physical_span).collect();
                                if !spans.is_empty() {
                                    let mean_span = spans.iter().sum::<i64>() as f64 / spans.len() as f64;
                                    let mut sorted_spans = spans.clone();
                                    sorted_spans.sort();
                                    let median_span = sorted_spans[sorted_spans.len() / 2];
                                    println!("  Physical span - mean: {:.1} bp, median: {:.0} bp", mean_span, median_span);
                                }
                                let mut type_counts = std::collections::HashMap::new();
                                for block in &blocks {
                                    if let Some(ref split_type) = block.split_type {
                                        *type_counts.entry(split_type.clone()).or_insert(0) += 1;
                                    }
                                }
                                println!("  Split types: {:?}", type_counts);
                            }
                        }
                    }
                }

                Method::FixedKb => {
                    let window_bp = config.window_bp as i64;
                    let mut i = 0;
                    while i < chr_snps.len() {
                        let start_pos = chr_snps[i].position;
                        let end_pos = start_pos + window_bp;

                        let snps_in_window: Vec<usize> = (i..chr_snps.len())
                            .take_while(|&j| chr_snps[j].position < end_pos)
                            .collect();

                        let n = snps_in_window.len();
                        if n >= config.min_snps {
                            blocks.push(Block {
                                chr: chr_num,
                                start_pos,
                                end_pos: chr_snps[i + n - 1].position,
                                start_idx: i,
                                end_idx: i + n - 1,
                                n_snps: n,
                                mean_ld_r2: None,
                                criterion_b_score: None,
                                selection_type: None,
                                original_block_size: None,
                                physical_span: Some(chr_snps[i + n - 1].position - start_pos),
                                split_type: Some("fixed_kb".to_string()),
                                rare_allele_count: None,
                                pic: None,
                            });
                        }
                        i += n.max(1);
                    }
                }
            }

            if config.verbose {
                println!("  Defined blocks: {}", blocks.len());
            }

            Ok((chr_num, blocks))
        })
        .collect();

    // collect results to BTreeMap
    let mut all_blocks = BTreeMap::new();
    for result in results {
        let (chr_num, blocks) = result?;
        all_blocks.insert(chr_num, blocks);
    }

    Ok(all_blocks)
}

pub fn deduplicate_blocks(
    mut all_blocks: BTreeMap<i32, Vec<Block>>,
    verbose: bool,
) -> BTreeMap<i32, Vec<Block>> {
    let mut total_removed = 0;

    for (chr_num, blocks) in all_blocks.iter_mut() {
        if blocks.is_empty() {
            continue;
        }

        let mut seen_ranges = HashSet::new();
        let mut unique_blocks = Vec::new();
        let mut duplicates = 0;

        for block in blocks.iter() {
            let snp_range: Vec<usize> = (block.start_idx..=block.end_idx).collect();
            let range_key = format!("{:?}", snp_range);

            if !seen_ranges.contains(&range_key) {
                seen_ranges.insert(range_key);
                unique_blocks.push(block.clone());
            } else {
                duplicates += 1;
            }
        }

        let original_len = blocks.len();
        *blocks = unique_blocks;
        total_removed += duplicates;

        if verbose && duplicates > 0 {
            println!("  Chr {}: removed {} duplicate block(s) ({} → {})",
                     chr_num, duplicates, original_len, blocks.len());
        }
    }

    if verbose && total_removed > 0 {
        println!("\nTotal duplicates removed: {}", total_removed);
    }

    all_blocks
}

pub fn ld_prune_blocks(
    all_blocks: BTreeMap<i32, Vec<Block>>,
    hap_files: &[std::path::PathBuf],
    config: &LdPruneConfig,
    noheader: bool,
    verbose: bool,
) -> BTreeMap<i32, Vec<Block>> {
    let mut pruned_all = BTreeMap::new();

    for (chr_num, blocks) in &all_blocks {
        if blocks.is_empty() {
            pruned_all.insert(*chr_num, vec![]);
            continue;
        }

        // Load haplotype matrix for this chromosome
        let hap_file = hap_files.iter().find(|f| get_chromosome_number(f) == *chr_num);
        let hap_matrix = match hap_file {
            Some(f) => match read_haplotype_file(f, noheader, false) {
                Ok(m) => m,
                Err(_) => {
                    pruned_all.insert(*chr_num, blocks.clone());
                    continue;
                }
            },
            None => {
                pruned_all.insert(*chr_num, blocks.clone());
                continue;
            }
        };

        // Compute dosage vectors and rare allele counts for all blocks
        let dosages: Vec<Vec<f64>> = blocks.iter()
            .map(|b| {
                let indices: Vec<usize> = (b.start_idx..=b.end_idx).collect();
                mh_to_dosage(&hap_matrix, &indices)
            })
            .collect();

        let rare_counts: Vec<usize> = blocks.iter()
            .map(|b| {
                let indices: Vec<usize> = (b.start_idx..=b.end_idx).collect();
                count_rare_alleles(&hap_matrix, &indices, config.rare_freq_threshold)
            })
            .collect();

        // Greedy pruning: sort by rare_count desc, then PIC/n_snps desc
        let mut priority: Vec<usize> = (0..blocks.len()).collect();
        priority.sort_by(|&a, &b| {
            rare_counts[b].cmp(&rare_counts[a])
                .then(blocks[b].n_snps.cmp(&blocks[a].n_snps))
        });

        let mut kept = vec![false; blocks.len()];
        let mut pruned_count = 0;

        for &i in &priority {
            if kept[i] { continue; }

            // Check r² with already-kept blocks within window
            let mut redundant = false;
            for j in 0..blocks.len() {
                if !kept[j] { continue; }

                // Only compare within genomic window
                let dist = (blocks[i].start_pos - blocks[j].start_pos).abs();
                if dist > config.window_bp { continue; }

                let r2 = compute_r2(&dosages[i], &dosages[j]);
                if r2 > config.r2_threshold {
                    redundant = true;
                    break;
                }
            }

            if !redundant {
                kept[i] = true;
            } else {
                pruned_count += 1;
            }
        }

        let pruned_blocks: Vec<Block> = blocks.iter().enumerate()
            .filter(|(i, _)| kept[*i])
            .map(|(_, b)| b.clone())
            .collect();

        if verbose {
            println!("  Chr {}: {} → {} blocks ({} pruned, r²>{:.2})",
                chr_num, blocks.len(), pruned_blocks.len(), pruned_count, config.r2_threshold);
        }

        pruned_all.insert(*chr_num, pruned_blocks);
    }

    pruned_all
}

pub fn load_block_definitions(
    block_dir: &Path,
    snp_map: &[SnpInfo],
    verbose: bool,
) -> Result<BTreeMap<i32, Vec<Block>>> {
    let mut all_blocks: BTreeMap<i32, Vec<Block>> = BTreeMap::new();

    let mut entries: Vec<_> = std::fs::read_dir(block_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.file_name().to_str()
                .map(|s| s.starts_with("hap_block_") && !s.contains("info"))
                .unwrap_or(false)
        })
        .collect();
    entries.sort_by_key(|e| get_chromosome_number(&e.path()));

    for entry in entries {
        let path = entry.path();
        let chr_num = get_chromosome_number(&path);
        let chr_snps: Vec<&SnpInfo> = snp_map.iter().filter(|s| s.chr == chr_num).collect();

        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        let mut blocks = Vec::new();

        for line in reader.lines() {
            let line = line?;
            let parts: Vec<&str> = line.splitn(2, '\t').collect();
            if parts.len() < 2 { continue; }

            let indices: Vec<usize> = parts[1].split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if indices.is_empty() { continue; }

            let start_idx = *indices.first().unwrap();
            let end_idx   = *indices.last().unwrap();

            blocks.push(Block {
                chr: chr_num,
                start_pos: chr_snps.get(start_idx).map_or(0, |s| s.position),
                end_pos:   chr_snps.get(end_idx).map_or(0, |s| s.position),
                start_idx,
                end_idx,
                n_snps: indices.len(),
                mean_ld_r2: None,
                criterion_b_score: None,
                selection_type: None,
                original_block_size: None,
                physical_span: None,
                split_type: Some("reused".to_string()),
                rare_allele_count: None,
                pic: None,
            });
        }

        if verbose {
            println!("  Chr {}: loaded {} blocks from {}", chr_num, blocks.len(), path.display());
        }
        all_blocks.insert(chr_num, blocks);
    }

    Ok(all_blocks)
}