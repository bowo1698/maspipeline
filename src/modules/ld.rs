// src/modules/ld.rs
use ndarray::Array2;
use super::types::Haploblock;

pub fn calculate_d_prime(snp1: &[f64], snp2: &[f64]) -> f64 {
    let pairs: Vec<(f64, f64)> = snp1
        .iter()
        .zip(snp2.iter())
        .filter(|(a, b)| !a.is_nan() && !b.is_nan())
        .map(|(a, b)| (*a, *b))
        .collect();

    if pairs.len() < 10 {
        return 0.0;
    }

    let n = pairs.len() as f64;
    let p_a: f64 = pairs.iter().map(|(a, _)| a).sum::<f64>() / n;
    let p_b: f64 = pairs.iter().map(|(_, b)| b).sum::<f64>() / n;
    let p_ab: f64 = pairs.iter().map(|(a, b)| a * b).sum::<f64>() / n;

    let d = p_ab - (p_a * p_b);

    let d_max = if d >= 0.0 {
        f64::min(p_a * (1.0 - p_b), (1.0 - p_a) * p_b)
    } else {
        f64::min(p_a * p_b, (1.0 - p_a) * (1.0 - p_b))
    };

    if d_max == 0.0 {
        0.0
    } else {
        d.abs() / d_max
    }
}

pub fn build_ld_haploblocks(
    hap_matrix: &Array2<f64>,
    snp_positions: Option<&[i64]>,
    d_prime_threshold: f64,
    min_block_snps: usize,
    verbose: bool,
) -> Vec<Haploblock> {
    let n_snps = hap_matrix.ncols();

    if verbose {
        println!("    Building haploblocks with D' threshold {}", d_prime_threshold);
    }

    let mut d_primes = Vec::new();
    let mut prev_col = hap_matrix.column(0).to_vec();
    for i in 0..n_snps - 1 {
        let next_col = hap_matrix.column(i + 1).to_vec();
        let d_prime = calculate_d_prime(&prev_col, &next_col);
        d_primes.push(d_prime);
        prev_col = next_col;  // reuse
    }

    if verbose {
        let mean_d = d_primes.iter().sum::<f64>() / d_primes.len() as f64;
        let above_thresh = d_primes.iter().filter(|&&d| d >= d_prime_threshold).count();
        println!("    Mean consecutive D': {:.3}", mean_d);
        println!("    Pairs ≥ threshold: {}/{}", above_thresh, d_primes.len());
    }

    let mut haploblocks = Vec::new();
    let mut current_block = vec![0];

    for (i, &d_prime) in d_primes.iter().enumerate() {
        if d_prime >= d_prime_threshold {
            current_block.push(i + 1);
        } else {
            if current_block.len() >= min_block_snps {
                let block_start = current_block[0];
                let block_end = *current_block.last().unwrap();

                let physical_span = if let Some(pos) = snp_positions {
                    Some(pos[block_end] - pos[block_start])
                } else {
                    None
                };

                haploblocks.push(Haploblock {
                    start_idx: block_start,
                    end_idx: block_end,
                    snp_indices: current_block.clone(),
                    n_snps: current_block.len(),
                    physical_span,
                });
            }
            current_block = vec![i + 1];
        }
    }

    if current_block.len() >= min_block_snps {
        let block_start = current_block[0];
        let block_end = *current_block.last().unwrap();

        let physical_span = if let Some(pos) = snp_positions {
            Some(pos[block_end] - pos[block_start])
        } else {
            None
        };

        haploblocks.push(Haploblock {
            start_idx: block_start,
            end_idx: block_end,
            snp_indices: current_block.clone(),
            n_snps: current_block.len(),
            physical_span,
        });
    }

    if verbose {
        let sizes: Vec<usize> = haploblocks.iter().map(|b| b.n_snps).collect();
        println!("    Built {} haploblocks", haploblocks.len());
        if !sizes.is_empty() {
            let mean = sizes.iter().sum::<usize>() as f64 / sizes.len() as f64;
            let mut sorted = sizes.clone();
            sorted.sort();
            let median = sorted[sorted.len() / 2];
            println!("    Block size - mean: {:.1}, median: {}, range: [{}, {}]",
                     mean, median, sizes.iter().min().unwrap(), sizes.iter().max().unwrap());

            if let Some(block) = haploblocks.first() {
                if let Some(_) = block.physical_span {
                    let spans: Vec<i64> = haploblocks.iter()
                        .filter_map(|b| b.physical_span)
                        .collect();
                    let mean_span = spans.iter().sum::<i64>() as f64 / spans.len() as f64;
                    let mut sorted_spans = spans.clone();
                    sorted_spans.sort();
                    let median_span = sorted_spans[sorted_spans.len() / 2];
                    println!("    Physical span (bp) - mean: {:.0}, median: {:.0}, range: [{:.0}, {:.0}]",
                             mean_span, median_span, 
                             spans.iter().min().unwrap(), 
                             spans.iter().max().unwrap());
                }
            }
        }
    }

    haploblocks
}

/// Calculate r² between two SNP slices
fn compute_r2_consecutive(a: &[f64], b: &[f64]) -> f64 {
    let pairs: Vec<(f64, f64)> = a.iter().zip(b.iter())
        .filter(|(x, y)| !x.is_nan() && !y.is_nan())
        .map(|(x, y)| (*x, *y))
        .collect();
    if pairs.len() < 10 { return 0.0; }
    let n = pairs.len() as f64;
    let pa = pairs.iter().map(|(x, _)| x).sum::<f64>() / n;
    let pb = pairs.iter().map(|(_, y)| y).sum::<f64>() / n;
    let pab = pairs.iter().map(|(x, y)| x * y).sum::<f64>() / n;
    let d = pab - pa * pb;
    let denom = pa * (1.0 - pa) * pb * (1.0 - pb);
    if denom <= 0.0 { 0.0 } else { (d * d / denom).min(1.0) }
}

/// Adaptively estimate r² threshold and max_snps from data.
///
/// Algorithm:
/// 1. Calculate r² for all consecutive pairs (i, i+1) → population r² distribution
/// 2. Find bimodal valley as threshold; fallback to 75th percentile if not bimodal
/// 3. For each SNP i, calculate "reach": how many SNPs to the right still have r² ≥ threshold
/// 4. max_snps = median(reach) per chromosome, capped between min_snps and 20
///
/// Return: (r2_threshold, max_snps)
pub fn estimate_adaptive_max_snps(
    hap_matrix: &Array2<f64>,
    min_snps: usize,
    verbose: bool,
) -> (f64, usize) {
    let n_snps = hap_matrix.ncols();

    // Step 1: r² konsekutif
    let mut r2_consec: Vec<f64> = Vec::with_capacity(n_snps - 1);
    let cols: Vec<Vec<f64>> = (0..n_snps)
        .map(|i| hap_matrix.column(i).to_vec())
        .collect();

    for i in 0..n_snps - 1 {
        r2_consec.push(compute_r2_consecutive(&cols[i], &cols[i + 1]));
    }

    // Step 2: threshold adaptif dari distribusi r² konsekutif
    let mut sorted_r2 = r2_consec.clone();
    sorted_r2.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let r2_threshold = {
        // Cari valley bimodal: scan KDE sederhana dengan 20 bins
        let n_bins = 20usize;
        let mut bins = vec![0usize; n_bins];
        for &v in &sorted_r2 {
            let idx = ((v * n_bins as f64) as usize).min(n_bins - 1);
            bins[idx] += 1;
        }
        // Cari lembah antara dua puncak
        let peak1 = bins[..n_bins/2].iter().enumerate().max_by_key(|(_, &v)| v).map(|(i, _)| i);
        let peak2 = bins[n_bins/2..].iter().enumerate().max_by_key(|(_, &v)| v).map(|(i, _)| i + n_bins/2);
        let valley = if let (Some(p1), Some(p2)) = (peak1, peak2) {
            if p1 < p2 {
                let valley_idx = (p1..p2).min_by_key(|&i| bins[i]).unwrap_or(p1);
                Some((valley_idx as f64 + 0.5) / n_bins as f64)
            } else { None }
        } else { None };

        // Fallback ke persentil-75 jika tidak bimodal
        valley.unwrap_or_else(|| sorted_r2[sorted_r2.len() * 3 / 4])
    };

    // Step 3: hitung reach tiap SNP
    let mut reaches: Vec<usize> = Vec::new();
    let scan_window = 50usize;
    for i in 0..n_snps {
        let limit = (i + scan_window).min(n_snps);
        let mut reach = 0usize;
        for j in (i + 1)..limit {
            if compute_r2_consecutive(&cols[i], &cols[j]) >= r2_threshold {
                reach = j - i;
            } else {
                break;
            }
        }
        if reach > 0 { reaches.push(reach); }
    }

    let max_snps = if reaches.is_empty() {
        min_snps.max(2)
    } else {
        let mut sorted_reach = reaches.clone();
        sorted_reach.sort();
        let median = sorted_reach[sorted_reach.len() / 2];
        median.max(min_snps).min(20)
    };

    if verbose {
        println!("    Adaptive r² threshold: {:.3}", r2_threshold);
        println!("    Adaptive max_snps: {}", max_snps);
    }

    (r2_threshold, max_snps)
}