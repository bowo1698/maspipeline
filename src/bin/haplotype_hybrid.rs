// src/bin/haplotype_hybrid.rs
use clap::{Parser, Subcommand};
use std::path::Path;
use anyhow::Result;
use microhaplotype_preprocessing::modules::*;

// =============================================================================
// GLOBAL ARGS
// =============================================================================

#[derive(Parser, Debug)]
#[command(name = "haplotype-hybrid")]
#[command(about = "Haplotype block definition for PHASED data")]
#[command(override_usage = "haplotype-hybrid -i <INPUT>... -m <MAP> -o <OUTPUT> [--generate-genotypes <DIR>] [OPTIONS] <COMMAND>")]
#[command(long_about = "
Haplotype block definition using hybrid approach for PHASED data.

Three methods available:
  ld-haploblock   LD-based blocks
  fixed-kb        Fixed physical distance blocks
  snp-count       Fixed SNP count per block

Examples:
  haplotype-hybrid -i hap_chr{1..18} -m map.txt -o out_micro --generate-genotypes out_micro  ld-haploblock micro --window-bp 125
  haplotype-hybrid -i hap_chr{1..18} -m map.txt -o out_pure  --generate-genotypes out_pure   ld-haploblock pure  --window 4
  haplotype-hybrid -i hap_chr{1..18} -m map.txt -o out_fixkb --generate-genotypes out_fixkb  fixed-kb --window-bp 100000
  haplotype-hybrid -i hap_chr{1..18} -m map.txt -o out_snp   --generate-genotypes out_snp    snp-count --window 4

  # With missing data handling:
  haplotype-hybrid -i hap_chr{1..18} -m map.txt -o out_micro --generate-genotypes out_micro --missing N --missing-value -9999  ld-haploblock micro --window-bp 125

  # Reuse block definitions from previous generation (e.g. G1 → G2):
  haplotype-hybrid -i hap_chrG2{1..18} -m map.txt -o out_G2 --generate-genotypes out_G2 --reuse-blocks out_G1/mh_info_ld_haploblock --save-allele-map out_G2/allele_maps  ld-haploblock micro --window-bp 125
")]
struct Cli {
    // --- General args ---
    #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
    #[arg(help = "Input haplotype files (supports {1..18} expansion)")]
    input: Vec<String>,

    #[arg(short = 'm', long = "map", required = true)]
    #[arg(help = "SNP map file (SNPID, chr, position)")]
    map: String,

    #[arg(long = "ncores", default_value_t = 1)]
    #[arg(help = "Number of parallel threads")]
    ncores: usize,

    #[arg(short = 'o', long = "output", default_value = "hap_info_microhap")]
    #[arg(help = "Output directory")]
    output: String,

    #[arg(long = "noheader", default_value_t = true)]
    #[arg(help = "Input files have no header row")]
    noheader: bool,

    #[arg(long = "generate-genotypes")]
    #[arg(help = "Generate haplotype genotype files in this subdirectory")]
    generate_genotypes: Option<String>,

    #[arg(long = "save-allele-map")]
    #[arg(help = "Save allele mapping per chromosome to this directory")]
    save_allele_map: Option<String>,

    #[arg(long = "reuse-blocks")]
    #[arg(help = "Reuse block definitions from this directory (hap_block_chr* files), skipping block discovery")]
    reuse_blocks: Option<String>,

    #[arg(long = "missing")]
    #[arg(help = "Missing allele code in haplotype files")]
    missing: Option<String>,

    #[arg(long = "missing-value", default_value = "-9999")]
    #[arg(help = "Output code for missing haplotypes")]
    missing_value: String,

    #[arg(long = "no-dedup")]
    #[arg(help = "Skip deduplication of overlapping blocks")]
    no_dedup: bool,

    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    // --- Post-hoc filtering (available for all methods) ---
    #[arg(long = "ld-prune", help_heading = "Advanced Filtering")]
    #[arg(help = "Apply LD pruning between blocks [default: off]")]
    ld_prune: bool,

    #[arg(long = "ld-prune-r2", default_value_t = 0.8, help_heading = "Advanced Filtering")]
    #[arg(help = "r² threshold for LD pruning between blocks [default: 0.8]")]
    ld_prune_r2: f64,

    #[arg(long = "ld-prune-window-bp", default_value_t = 500_000, help_heading = "Advanced Filtering")]
    #[arg(help = "Genomic window (bp) for LD pruning [default: 500000]")]
    ld_prune_window_bp: i64,

    #[arg(long = "no-rare-priority", help_heading = "Advanced Filtering")]
    #[arg(help = "Disable rare allele priority during LD pruning [default: off, rare priority active]")]
    no_rare_priority: bool,

    #[arg(long = "min-pic", help_heading = "Advanced Filtering")]
    #[arg(help = "Minimum PIC score to retain a block [default: none, no filter]")]
    min_pic: Option<f64>,

    #[arg(long = "top-k", help_heading = "Advanced Filtering")]
    #[arg(help = "Retain only top-K blocks ranked by PIC [default: none, keep all]")]
    top_k: Option<usize>,

    #[arg(long = "min-effective-alleles", help_heading = "Advanced Filtering")]
    #[arg(help = "Minimum effective number of alleles per block [default: none, no filter]")]
    min_effective_alleles: Option<f64>,

    #[arg(long = "max-rare-alleles-prop", help_heading = "Advanced Filtering")]
    #[arg(help = "Maximum proportion of rare alleles (freq < 0.05) per block [default: none, no filter]")]
    max_rare_alleles_prop: Option<f64>,

    #[command(subcommand)]
    command: Commands,
}

// TOP-LEVEL SUBCOMMANDS

#[derive(Subcommand, Debug)]
enum Commands {
    /// LD-based haplotype blocks using D' threshold
    #[command(name = "ld-haploblock")]
    LdHaploblock {
        #[arg(long = "d-prime-threshold", default_value_t = 0.45)]
        #[arg(help = "D' threshold for LD haploblock definition")]
        d_prime_threshold: f64,

        #[arg(long = "min-snps", default_value_t = 2)]
        #[arg(help = "Minimum SNPs per block")]
        min_snps: usize,

        #[arg(long = "aft", default_value_t = 0.08)]
        #[arg(help = "Allele frequency threshold for Criterion-B")]
        aft: f64,

        #[arg(long = "md", default_value_t = 0.10)]
        #[arg(help = "Minimum divergence parameter for Criterion-B")]
        md: f64,

        #[arg(long = "max-criterion-b")]
        #[arg(help = "Filter blocks with Criterion-B score above this threshold")]
        max_criterion_b: Option<f64>,

        #[command(subcommand)]
        haplotype_type: HaplotypeTypeCmd,
    },

    /// Fixed physical distance blocks
    #[command(name = "fixed-kb")]
    FixedKb {
        #[arg(long = "window-bp", required = true)]
        #[arg(help = "Block size in base pairs (e.g. 100000 = 100 kb)")]
        window_bp: i32,

        #[arg(long = "min-snps", default_value_t = 2)]
        #[arg(help = "Minimum SNPs required to form a block")]
        min_snps: usize,
    },

    /// Fixed SNP count per block
    #[command(name = "snp-count")]
    SnpCount {
        #[arg(short = 'w', long = "window", default_value_t = 4)]
        #[arg(help = "Number of SNPs per block")]
        window: usize,

        #[arg(long = "min-snps", default_value_t = 2)]
        #[arg(help = "Minimum SNPs required to retain a block")]
        min_snps: usize,
    },
}

#[derive(Subcommand, Debug)]
enum HaplotypeTypeCmd {
    /// Split LD blocks into fixed physical windows (microhaplotypes)
    #[command(name = "micro")]
    Micro {
        #[arg(long = "window-bp", default_value_t = 125)]
        #[arg(help = "Physical window size in bp for splitting LD blocks")]
        window_bp: i32,

        #[arg(long = "max-snps", default_value_t = usize::MAX)]
        #[arg(help = "Maximum SNPs per microhaplotype block [default: no limit]")]
        max_snps: usize,

        #[arg(long = "adaptive-max-snps", default_value_t = false)]
        #[arg(help = "Estimate max-snps adaptively from LD structure of the data [default: off]")]
        adaptive_max_snps: bool,
    },

    /// Select best fixed-size haplotype window from each LD block
    #[command(name = "pure")]
    Pure {
        #[arg(short = 'w', long = "window", default_value_t = 4)]
        #[arg(help = "Number of SNPs per haplotype (best window selected by Criterion-B)")]
        window: usize,
    },
}

// =============================================================================
// MAIN
// =============================================================================

fn main() -> Result<()> {
    let cli = Cli::parse();

    let config = build_config(&cli);

    // --- Print header ---
    println!("\n{}", "=".repeat(70));
    println!("MICROHAPLOTYPE BLOCK DEFINITION - HYBRID (PHASED DATA)");
    println!("{}", "=".repeat(70));
    println!("Input files:       {} files", cli.input.len());
    println!("Map file:          {}", cli.map);
    print_method_summary(&cli);
    println!("Output directory:  {}", cli.output);
    println!("{}", "=".repeat(70));

    // --- Setup rayon ---
    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.ncores)
        .build_global()
        .unwrap();

    if cli.verbose {
        println!("Threads:           {}", cli.ncores);
    }

    // --- Read inputs ---
    let snp_map = read_map_file(&cli.map)?;
    let hap_files = expand_and_sort_files(&cli.input)?;

    // --- Define blocks ---
    let t_define = std::time::Instant::now();
    let all_blocks = if let Some(ref block_dir) = cli.reuse_blocks {
        println!("Reusing block definitions from: {}", block_dir);
        load_block_definitions(Path::new(block_dir), &snp_map, cli.verbose)?
    } else {
        define_microhaplotype_blocks(&hap_files, &snp_map, &config)?
    };
    println!("[timing] define_blocks:     {:.1}s", t_define.elapsed().as_secs_f64());

    // --- Deduplication ---
    let all_blocks = if !cli.no_dedup {
        deduplicate_blocks(all_blocks, cli.verbose)
    } else {
        all_blocks
    };

    // --- LD pruning ---
    let all_blocks = if cli.ld_prune {
        if cli.verbose {
            println!("\nLD pruning blocks (r²>{}, window={}bp)...", cli.ld_prune_r2, cli.ld_prune_window_bp);
        }
        let prune_config = LdPruneConfig {
            r2_threshold: cli.ld_prune_r2,
            window_bp: cli.ld_prune_window_bp,
            prioritize_rare: !cli.no_rare_priority,
            rare_freq_threshold: 0.05,
        };
        ld_prune_blocks(all_blocks, &hap_files, &prune_config, cli.noheader, cli.verbose)
    } else {
        all_blocks
    };

     // --- Allele quality filter ---
    let all_blocks = if cli.min_effective_alleles.is_some() || cli.max_rare_alleles_prop.is_some() {
        filter_blocks_by_allele_quality(
            all_blocks, &hap_files,
            cli.min_effective_alleles, cli.max_rare_alleles_prop,
            cli.noheader, cli.verbose,
        )
    } else {
        all_blocks
    };

    // --- PIC ranking ---
    let all_blocks = if cli.min_pic.is_some() || cli.top_k.is_some() {
        rank_and_filter_blocks(
            all_blocks, &hap_files,
            cli.min_pic, cli.top_k,
            cli.noheader, cli.verbose,
        )
    } else {
        all_blocks
    };

    // --- Outputs ---
    let method_label = match config.method {
        Method::LdHaploblock => match config.haplotype_type {
            HaplotypeType::Pure  => "ld_haploblock_pure",
            HaplotypeType::Micro => "ld_haploblock_micro",
        },
        Method::SnpCountSimple => "snp_count_simple",
        Method::FixedKb        => "fixed_kb",
    };

    let t_csv = std::time::Instant::now();
    write_csv_outputs(&all_blocks, &hap_files, &snp_map, &cli.output,
                      method_label, cli.noheader, cli.verbose)?;
    println!("[timing] write_csv_outputs: {:.1}s", t_csv.elapsed().as_secs_f64());

    let t_cb = std::time::Instant::now();
    write_criterion_b_summary(&all_blocks, &cli.output, cli.verbose)?;
    println!("[timing] criterion_b_summary: {:.1}s", t_cb.elapsed().as_secs_f64());

    let t_snp = std::time::Instant::now();
    write_snp_selection_file(&all_blocks, &snp_map, &cli.output, cli.verbose)?;
    println!("[timing] snp_selection:     {:.1}s", t_snp.elapsed().as_secs_f64());

    println!("\nWriting output files:");
    let total_blocks = write_output_files(&all_blocks, &cli.output, cli.verbose)?;

    let t_geno = std::time::Instant::now();
    if let Some(ref geno_prefix) = cli.generate_genotypes {
        let parent_dir = Path::new(&cli.output).parent().unwrap_or(Path::new("."));
        let geno_output_dir = if geno_prefix == Path::new(&cli.output)
            .file_name().unwrap().to_str().unwrap()
        {
            parent_dir.join(format!("{}_geno", geno_prefix))
        } else {
            parent_dir.join(geno_prefix)
        };
        generate_haplotype_genotypes(
            &all_blocks, &hap_files, &geno_output_dir,
            cli.missing.as_deref(), &cli.missing_value,
            cli.noheader, cli.verbose,
            cli.save_allele_map.as_deref().map(Path::new),
        )?;
    }
    println!("[timing] generate_genotypes: {:.1}s", t_geno.elapsed().as_secs_f64());

    let t_cov = std::time::Instant::now();
    let coverage_stats = calculate_genome_coverage(&all_blocks, &snp_map, cli.verbose)?;
    save_coverage_stats(&coverage_stats, &cli.output, method_label, cli.verbose)?;
    println!("[timing] coverage:          {:.1}s", t_cov.elapsed().as_secs_f64());

    println!("\n{}", "=".repeat(70));
    println!("SUMMARY");
    println!("{}", "=".repeat(70));
    println!("Chromosomes processed:  {}", all_blocks.len());
    println!("Total blocks defined:   {}", total_blocks);
    println!("Output directory:       {}/", cli.output);
    println!("{}", "=".repeat(70));
    println!("\n✓ Block definition completed!");

    Ok(())
}

// =============================================================================
// HELPERS
// =============================================================================

fn build_config(cli: &Cli) -> BlockDefinitionConfig {
    match &cli.command {
        Commands::LdHaploblock {
            d_prime_threshold, min_snps, aft, md, max_criterion_b,
            haplotype_type,
        } => {
            let (method_ht, window_bp, window_snps, max_snps, adaptive_max_snps) = match haplotype_type {
                HaplotypeTypeCmd::Micro { window_bp, max_snps, adaptive_max_snps } =>
                    (HaplotypeType::Micro, *window_bp, 4, *max_snps, *adaptive_max_snps),
                HaplotypeTypeCmd::Pure { window } =>
                    (HaplotypeType::Pure, 125, *window, usize::MAX, false),
            };
            BlockDefinitionConfig {
                method: Method::LdHaploblock,
                haplotype_type: method_ht,
                window_bp,
                window_snps,
                min_snps: *min_snps,
                max_snps,
                adaptive_max_snps,
                d_prime_threshold: *d_prime_threshold,
                aft: *aft,
                md: *md,
                min_ld: None,
                max_criterion_b: *max_criterion_b,
                noheader: cli.noheader,
                verbose: cli.verbose,
            }
        }
        Commands::FixedKb { window_bp, min_snps } => {
            BlockDefinitionConfig {
                method: Method::FixedKb,
                haplotype_type: HaplotypeType::Micro,
                window_bp: *window_bp,
                window_snps: 4,
                min_snps: *min_snps,
                max_snps: usize::MAX,
                adaptive_max_snps: false,
                d_prime_threshold: 0.0,
                aft: 0.0,
                md: 0.0,
                min_ld: None,
                max_criterion_b: None,
                noheader: cli.noheader,
                verbose: cli.verbose,
            }
        }
        Commands::SnpCount { window, min_snps } => {
            BlockDefinitionConfig {
                method: Method::SnpCountSimple,
                haplotype_type: HaplotypeType::Micro,
                window_bp: 0,
                window_snps: *window,
                min_snps: *min_snps,
                max_snps: usize::MAX,
                adaptive_max_snps: false,
                d_prime_threshold: 0.0,
                aft: 0.0,
                md: 0.0,
                min_ld: None,
                max_criterion_b: None,
                noheader: cli.noheader,
                verbose: cli.verbose,
            }
        }
    }
}

fn print_method_summary(cli: &Cli) {
    match &cli.command {
        Commands::LdHaploblock {
            d_prime_threshold, aft, md, max_criterion_b,
            haplotype_type, ..
        } => {
            println!("Method:            LD-haploblock");
            println!("D' threshold:      {}", d_prime_threshold);
            println!("Criterion-B AFT:   {}  MD: {}", aft, md);
            match haplotype_type {
                HaplotypeTypeCmd::Micro { window_bp, max_snps, adaptive_max_snps } => {
                    println!("Haplotype type:    micro ({}bp windows)", window_bp);
                    if *adaptive_max_snps {
                        println!("Max SNPs/block:    adaptive (estimated from LD structure)");
                    } else if *max_snps != usize::MAX {
                        println!("Max SNPs/block:    {}", max_snps);
                    }
                }
                HaplotypeTypeCmd::Pure { window } => {
                    println!("Haplotype type:    pure ({} SNPs, best Criterion-B)", window);
                }
            }
            if let Some(b) = max_criterion_b { println!("Max Criterion-B:   {}", b); }
        }
        Commands::FixedKb { window_bp, min_snps } => {
            println!("Method:            fixed-kb");
            println!("Window size:       {} bp ({} kb)", window_bp, window_bp / 1000);
            println!("Min SNPs/block:    {}", min_snps);
        }
        Commands::SnpCount { window, min_snps } => {
            println!("Method:            snp-count");
            println!("Window size:       {} SNPs", window);
            println!("Min SNPs/block:    {}", min_snps);
        }
    }
    // Print global post-hoc filters if active
    if cli.ld_prune {
        println!("LD pruning:        r²>{}, window={}bp", cli.ld_prune_r2, cli.ld_prune_window_bp);
    }
    if let Some(p) = cli.min_pic             { println!("Min PIC:           {}", p); }
    if let Some(k) = cli.top_k               { println!("Top-K blocks:      {}", k); }
    if let Some(v) = cli.min_effective_alleles { println!("Min eff. alleles:  {}", v); }
    if let Some(v) = cli.max_rare_alleles_prop { println!("Max rare prop:     {}", v); }
    if let Some(ref d) = cli.reuse_blocks      { println!("Reuse blocks:      {}", d); }
    if let Some(ref d) = cli.save_allele_map   { println!("Save allele map:   {}", d); }
}