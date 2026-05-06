// src/modules/types.rs
use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub struct SnpInfo {
    pub snpid: String,
    pub chr: i32,
    pub position: i64,
}

#[derive(Debug, Clone)]
pub struct Block {
    pub chr: i32,
    pub start_pos: i64,
    pub end_pos: i64,
    pub start_idx: usize,
    pub end_idx: usize,
    pub n_snps: usize,
    pub mean_ld_r2: Option<f64>,
    pub criterion_b_score: Option<f64>,
    pub selection_type: Option<String>,
    pub original_block_size: Option<usize>,
    pub physical_span: Option<i64>,
    pub split_type: Option<String>,
    pub rare_allele_count: Option<usize>,
    pub pic: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct Haploblock {
    pub start_idx: usize,
    pub end_idx: usize,
    pub snp_indices: Vec<usize>,
    pub n_snps: usize,
    pub physical_span: Option<i64>,
}

#[derive(Debug, Clone)]
pub struct SelectedHaplotype {
    pub snp_indices: Vec<usize>,
    pub n_snps: usize,
    pub criterion_b_score: f64,
    pub original_block_size: usize,
    pub selection_type: String,
}

#[derive(Debug, Clone)]
pub struct Microhaplotype {
    pub snp_indices: Vec<usize>,
    pub n_snps: usize,
    pub criterion_b_score: f64,
    pub physical_span: i64,
    pub split_type: String,
}

#[derive(Debug, Clone)]
pub enum Method {
    LdHaploblock,
    SnpCountSimple,
    FixedKb,
}

#[derive(Debug, Clone)]
pub enum HaplotypeType {
    Pure,
    Micro,
}

/// Configuration for block definition
#[derive(Debug, Clone)]
pub struct BlockDefinitionConfig {
    pub method: Method,
    pub window_bp: i32,
    pub haplotype_type: HaplotypeType,
    pub window_snps: usize,
    pub min_snps: usize,
    pub max_snps: usize,
    pub d_prime_threshold: f64,
    pub aft: f64,
    pub md: f64,
    pub min_ld: Option<f64>,
    pub max_criterion_b: Option<f64>,
    pub adaptive_max_snps: bool,
    pub noheader: bool,
    pub verbose: bool,
}

impl Default for BlockDefinitionConfig {
    fn default() -> Self {
        Self {
            method: Method::LdHaploblock,
            window_bp: 125,
            haplotype_type: HaplotypeType::Pure,
            window_snps: 4,
            min_snps: 2,
            max_snps: 4,
            d_prime_threshold: 0.45,
            aft: 0.08,
            md: 0.10,
            min_ld: None,
            max_criterion_b: None,
            adaptive_max_snps: false,
            noheader: false,
            verbose: false,
        }
    }
}

/// COnfigure LD Pruning
#[derive(Debug, Clone)]
pub struct LdPruneConfig {
    pub r2_threshold: f64,
    pub window_bp: i64,
    pub prioritize_rare: bool,
    pub rare_freq_threshold: f64,
}

impl Default for LdPruneConfig {
    fn default() -> Self {
        Self {
            r2_threshold: 0.8,
            window_bp: 500_000,
            prioritize_rare: true,
            rare_freq_threshold: 0.05,
        }
    }
}

/// Coverage statistics
#[derive(Debug, Clone)]
pub struct CoverageStats {
    pub total_genome_length_bp: i64,
    pub covered_length_bp: i64,
    pub coverage_pct: f64,
    pub total_blocks: usize,
    pub chromosomes: BTreeMap<i32, ChromosomeCoverage>,
}

#[derive(Debug, Clone)]
pub struct ChromosomeCoverage {
    pub chr: i32,
    pub n_blocks: usize,
    pub length_bp: i64,
    pub covered_bp: i64,
    pub coverage_pct: f64,
}

#[derive(Debug, Clone)]
pub struct AlleleMetrics {
    pub n_alleles: usize,
    pub n_predictable: usize,
    pub min_freq: f64,
    pub max_freq: f64,
    pub mean_freq: f64,
    pub freq_entropy: f64,
    pub expected_heterozygosity: f64,
    pub observed_heterozygosity: f64,
    pub pic: f64,
    pub effective_alleles: f64,
    pub rare_alleles_prop: f64,
}