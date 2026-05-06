# Preprocessing pipeline for genomic prediction using microhaplotype markers

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Rust](https://img.shields.io/badge/Rust-1.9+-orange.svg)](https://www.rust-lang.org/)

## General information

We adapted and optimised the GVCHAP pipeline from [Prakapenka et al. (2020)](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00282/full) to handle large genomic data and utilise microhaplotype markers in genomic prediction. The pipeline consists of four key stages: (1) phasing, (2) converting phased data to haplotypes, (3) defining microhaplotype blocks, and finally, (4) encoding them into microhaplotype genotypes. The pipeline was optimised with Rust, so it has more efficient memory and enhanced scalability due to parallelisation. For phasing, we only used Beagle as the pipeline has not been developed for other phasing methods yet (e.g., FindHap, FImpute, etc.). In addition, we encourage the use of parallelisation techniques when phasing with Beagle. 

While haplotypes refer to long-range combinations of alleles along a chromosome, microhaplotypes can be defined as several SNPs covering a short DNA segment, typically 125 - 150 bp. So in this pipeline, conducted by `haplotype-hybrid`, we provide fixed window (SNP counts) and LD-based methods to define fully haplotype blocks and microhaplotype segments, respectively. When the `--method snp-count-simple` argument is used, it simply defines haplotype blocks depending on how many SNPs are set in each block. In contrast, when the `--method ld-haploblock` along with `--haplotype-type micro` arguments are used, it will discover microhaplotypes based on LD and window constraints. However, both methods will result in similar data outputs.

## Installation

We provide three tools for the microhaplotype discovery and genotyping preprocessing pipeline:

1.  **`convert-to-vcf`**: Converts genotype data from CSV format to standard VCF (Variant Call Format) file
    - Input: CSV file with SNP genotypes (0/1/2 coding)
    - Output: VCF file compatible with phasing software (e.g., Beagle)

2.  **`convert-from-vcf`**: Extracts phased haplotype data from VCF files
    - Input: Phased VCF file (after running phasing software)
    - Output: Separated haplotype files per chromosome for downstream analysis
    - Generates SNP map file for position information

3.  **`haplotype-hybrid`**: Discovers microhaplotype segments using LD-based haploblock identification and performs genotyping
    - Input: Phased haplotype files and SNP map
    - Methods: LD-based haploblock discovery (Jonas et al. 2017) or fixed-window segmentation
    - Output: Microhaplotype block definitions and numerical genotype matrix
    - Implements Criterion-B scoring for optimal haplotype selection

We provide two installation options:

1. **Build from source** (Recommended): Compile using Rust's Cargo package manager
2. **Pre-compiled binaries** (Quick start): Download platform-specific executables from the [Releases](https://github.com/bowo1698/MicrohapsSel/releases/tag/v1.0) page

> We strongly encourage building from source using Cargo, as different operating systems and hardware architectures may require specific optimisations.

### 1. Build from source

**Prerequisites for building from source**

- **Rust toolchain** (version 1.70 or later recommended): Install from [rustup.rs](https://rustup.rs)

**Build instructions**

```bash
# 1. Clone the repository
git clone https://github.com/bowo1698/MicrohapsSel.git
cd MicrohapsSel

# 2. Navigate to preprocessing directory
cd preprocessing

# 3. Build all tools in release mode (optimized)
cargo build --release

# 4. Compiled binaries will be available in:
# ./target/release/convert-to-vcf
# ./target/release/convert-from-vcf
# ./target/release/haplotype-hybrid

# 5. (Optional) Add to system PATH for global access
# For Linux/macOS:
export PATH="$PWD/target/release:$PATH"
# Add to ~/.bashrc or ~/.zshrc to make permanent

# For Windows (PowerShell):
$env:Path += ";$PWD\target\release"
# Add to system environment variables to make permanent

# 6. Verify installation
./convert-to-vcf --help
./convert-from-vcf --help
./haplotype-hybrid --help
```

### 2. Pre-compiled binaries

Download platform-specific binaries from the links below:

- **Linux (x86_64)**: [microhaplotype-tools-x64-linux.tar.gz](https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-x64-linux.tar.gz)
- **macOS (Intel)**: [microhaplotype-tools-x64-macos.tar.gz](https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-x64-macos.tar.gz)
- **macOS (Apple Silicon)**: [microhaplotype-tools-arm64-macos.tar.gz](https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-arm64-macos.tar.gz)
- **Windows (x86_64)**: [microhaplotype-tools-x64-windows.zip](https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-x64-windows.zip)

**Linux:**
```bash
# Download and extract
wget https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-x64-linux.tar.gz
tar -xzf microhaplotype-tools-x64-linux.tar.gz
cd microhaplotype-tools-x64-linux

# Make executable
chmod +x convert-to-vcf convert-from-vcf haplotype-hybrid

# Test installation
./convert-to-vcf --help
./convert-from-vcf --help
./haplotype-hybrid --help

# (Optional) Move to system PATH for global access
sudo mv convert-to-vcf convert-from-vcf haplotype-hybrid /usr/local/bin/
```

**macOS (Intel x86_64):**
```bash
# Download and extract
curl -L -O https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-x64-macos.tar.gz
tar -xzf microhaplotype-tools-x64-macos.tar.gz
cd microhaplotype-tools-x64-macos

# Remove macOS quarantine attribute and make executable
chmod +x convert-to-vcf convert-from-vcf haplotype-hybrid
xattr -d com.apple.quarantine convert-to-vcf convert-from-vcf haplotype-hybrid

# Test installation
./convert-to-vcf --help
./convert-from-vcf --help
./haplotype-hybrid --help

# (Optional) Move to system PATH for global access
sudo mv convert-to-vcf convert-from-vcf haplotype-hybrid /usr/local/bin/
```

**macOS (Apple silicon ARM64):**
```bash
# Download and extract
curl -L -O https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-arm64-macos.tar.gz
tar -xzf microhaplotype-tools-arm64-macos.tar.gz
cd microhaplotype-tools-arm64-macos

# Remove macOS quarantine attribute and make executable
chmod +x convert-to-vcf convert-from-vcf haplotype-hybrid
xattr -d com.apple.quarantine convert-to-vcf convert-from-vcf haplotype-hybrid

# Test installation
./convert-to-vcf --help
./convert-from-vcf --help
./haplotype-hybrid --help

# (Optional) Move to system PATH for global access
sudo mv convert-to-vcf convert-from-vcf haplotype-hybrid /usr/local/bin/
```

**Windows:**
```powershell
# Download (using PowerShell)
Invoke-WebRequest -Uri "https://github.com/bowo1698/MicrohapsSel/releases/download/v1.0/microhaplotype-tools-x64-windows.zip" -OutFile "microhaplotype-tools-x64-windows.zip"

# Extract
Expand-Archive -Path microhaplotype-tools-x64-windows.zip -DestinationPath microhaplotype-tools-x64-windows
cd microhaplotype-tools-x64-windows

# Test installation (no chmod needed on Windows)
.\convert-to-vcf.exe --help
.\convert-from-vcf.exe --help
.\haplotype-hybrid.exe --help

# (Optional) Add to system PATH:
# 1. Copy full path of current directory: Get-Location
# 2. Search "Environment Variables" in Windows Start menu
# 3. Edit "Path" variable → Add new entry with the copied path
```

**Note:** On macOS, if you still get security warnings after using `xattr`, go to System Preferences → Security & Privacy → General, and click "Allow Anyway" for each blocked binary.

## Data preparation

Raw genotype data is expected in a .VCF format, as it contains genetic map positions, reference/alternative alleles, and unphased genotype calls (e.g., 0/1), which serve as the essential input for the subsequent phasing step. However, if the raw geneotype data is in PLINK formats (provided as .bed, .bim, .fam), recoding to .VCF format can be done by PLINK as well, for example:

```bash
plink --bfile GENO_DATA --chr-set 95 no-xy --recode vcf --out GENO_DATA
```

This will generate unphased genotype data in a .VCF file, for example:

```
##fileformat=VCFv4.2
##fileDate=20250503
##source=PLINKv1.9
##contig=<ID=1,length=80480704>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Ind_1	Ind_2
1	16977	snp1	C	A	.	.	PR	GT	0/0	0/1
1	33723	snp2	G	T	.	.	PR	GT	0/0	0/0
```

Alternatively, if the genotype data is in form of .CSV:

```
ID,snp_1,snp_2,snp_3,...snp_n
ind_1,0,1,2,...
ind_2,1,2,0,...
ind_3,2,0,1,...
```

Along with a genetic map (in form .txt):

```
SNPID	Chr	Position
snp_1	1	0
snp_2	1	16746
snp_3	1	217334
```

We provide `convert-to-vcf` to convert such genotype data to a VCF file. This tool can be run as follows:

```
./convert-to-vcf \
      -i genotypes.csv \
      -m map.txt \
      -o genotypes.vcf \
      --split-by-chr \
      --parallel \
      --ncores 64 \
      --merge-after-split
```

Note that the conversion is done by splitting the genotype files based on the number of chromosomes processed in parallel, then merging them into a single VCF file.

## Phasing

Phasing is the first critical step in our pipeline, where ambiguous diploid genotypes are separated into their constituent haploid sequences to allow for the accurate identification of microhaplotype blocks. We encourage to do phasing using Beagle in parallel process, and to do that, we use `bcftools` to split the VCF file per chromosome and `tabix` for indexing a compressed VCF, as an example:

```bash
for chr in {1..29}; do
    bcftools view -r ${chr} /data/genotypes.vcf.gz -O z -o /data/tmp/chr${chr}.vcf.gz
    tabix -p vcf /data/tmp/chr${chr}.vcf.gz
done

for chr in {1..29}; do
    (
        java -Xmx28g -jar ${EBROOTBEAGLE}/beagle.jar \
            gt=/data/tmp/chr${chr}.vcf.gz \
            nthreads=4 \
            out=/data/tmp/chr${chr}_phased
    ) &
    
    if (( chr % 8 == 0 )); then
        wait
    fi
done
wait

for chr in {1..29}; do
    tabix -p vcf /data/tmp/chr${chr}_phased.vcf.gz
    bcftools reheader -h <(bcftools view -h /data/genotypes.vcf.gz) /data/tmp/chr${chr}_phased.vcf.gz | \
        bcftools view -O z -o /data/tmp/chr${chr}_phased_corrected.vcf.gz
    tabix -p vcf /data/tmp/chr${chr}_phased_corrected.vcf.gz
done

bcftools concat -O z -o /data/genotypes_phased.vcf.gz \
    /data/tmp/chr{1..29}_phased_corrected.vcf.gz

tabix -p vcf /data/genotypes_phased.vcf.gz

/usr/bin/rm -rf /data/tmp
gunzip /data/genotypes_phased.vcf.gz
```

This process results in a phased VCF file, for example:

```
##fileformat=VCFv4.2
##filedate=20260108
##source="beagle.22Jul22.46e.jar"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + 2*P(AA)]">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Ind_1	Ind_2
1	16977	snp1	C	A	.	.	PR	GT	0|0	1|1	
1	33723	snp2	G	T	.	.	PR	GT	1|0	0|0
```

## Convertion from phased genotypes to haplotypes

After phasing, we want to separate the phased haplotypes into individual haploid sequences for haploblock detection. Here, we provide `convert-from-vcf` for conversion, which accomplishes two key transformations:

**1. Haplotype separation**

Each phased genotype (e.g., `0|1`) is split into two separate haplotype columns:
- Individual_A with genotype `0|1` becomes:
  - Haplotype 1: `1` (recoded from `0`)
  - Haplotype 2: `2` (recoded from `1`)

This separation is essential because LD-based haploblock detection requires calculating pairwise $D'$ between SNPs, which depends on counting the some possible haplotype combinations (e.g. 0|0, 0|1, 1|0, 1|1).

**2. Allele recoding**

VCF format uses `0` (reference) and `1` (alternate) encoding. The conversion recodes these to:

- `0` → `1` (major allele)
- `1` → `2` (minor allele)

For genotypes (diploid coding):

- `0|0` → `0` (homozygous major)
- `0|1` or `1|0` → `1` (heterozygous)
- `1|1` → `2` (homozygous minor)

To handle large datasets, individuals are processed in chunks (batches) rather than all at once. The chunk size is determined by the `--interval` parameter:

- `--interval 1`: Process all individuals in one pass (high memory)
- `--interval 10`: Process in 10 batches (lower memory, slightly slower)

Each chunk reads through the VCF file sequentially, extracts the relevant individuals, and writes output incrementally. Data is automatically split by chromosome, creating separate output files (haplotypes and genotypes), for example:

**Haplotype Files (`hap/chr*.hap`)**

Format: Each individual contributes **two columns** (one per haplotype)
```
ID                snp1_h1  snp1_h2  snp2_h1  snp2_h2  ...
Ind_1             1        2        1        1        ...
Ind_2             2        2        1        2        ...
```

This provides input for LD-based haploblock detection

**Genotype Files (`geno/chr*.geno`)**

Format: Traditional genotype matrix with header
```
ID              snp1  snp2  snp3  ...
Ind_1           1     0     2     ...
Ind_2           2     1     1     ...
```

This provides input for GBLUP or other genotype-based prediction models that utilises biallelic SNPs.

### Usage Example

```bash
./convert-from-vcf \
      -i /data/genotypes_phased.vcf.gz \
      --map map.txt \
      --genofolder geno \
      --hapfolder hap \
      --nosort
```

### Key Parameters

- `-i, --input`: VCF file(s) to convert (required)
- `--interval`: Number of chunks for processing (default: 1)
- `--hapfolder`: Output directory for haplotype files (default: `hap`)
- `--genofolder`: Output directory for genotype files (default: `geno`)
- `--map`: Output path for SNP map file (default: `map_new.txt`)
- `--missing`: Code for missing genotypes (default: `-9999`)
- `-v, --verbose`: Display detailed progress information
- `-h, --help`: Print help information

The separated haplotype files (`hap/chr*`) then serve as direct input for the next step: defining microhaplotype segments using LD-based haploblock detection.

## Discovering microhaplotype segments and genotyping

Because our objective is to discover microhaplotype segments, we only use outputs from `hap/chr*` and perform arguments `--method ld-haploblock` and `--haplotype-type micro`. The `haplotype-hybrid` tool works with finding all potential haplotype candidates that meet an LD criteria (e.g., $D' > 0.45$) as the argument of `--d-prime-threshold 0.45` is applied. After that, align with the microhaplotype definition, they are selected from the haplotype block candidates, where in each segment should contain consecutive SNPs within 125 or 150 bp which then they are evaluated using Criterion-B following [Jónás et al. (2017)](https://www.journalofdairyscience.org/article/S0022-0302(16)30076-5/fulltext) to ensure balance between allele frequency and microhaplotype diversity, following the calculation:

$$CriterionB_{m h_i}=\sum_{k=1}^{N_i}\left(f_i-\frac{1}{H S}\right)^2-w N_i$$

In the first term, $m h_i$ denotes microhaplotype $i$, $f_i$ is the frequency of microhaplotype allele, and $HS=2^n$ is the theoretical maximum number of alleles for n SNPs. In the second term, $w$ is calculated as $(MD.N_i)/(HS.(N_i-1) )$, where $MD$ is the scaling parameter of maximum deviation to control the magnitude of microhaplotype diversity. We set the $MD$ as $0.1$ following (Jónás et al., 2017). The $N_i$ represents the number of predictable alleles for microhaplotype $i$ that have a frequency above the allele frequency threshold ($AFT≥0.08$). 

The final genotype data will be:

```
ID       hap_1_1  hap_1_1  hap_1_2  hap_1_2  hap_1_3  hap_1_3
ind_1        2        1        1        2        3        3
ind_2        1        3        2        2        3        2
ind_3        3        1        2        2        2        3
```

This genotype data consists of numerical microhaplotype alleles in diploid format, where each number represents a unique combination of SNP sequences in a specific genomic block. Each locus is presented in two side-by-side columns to show the pair of haplotypes inherited from each individual's parents. 

For example, at locus `hap_1_1`, individual `ind_1` has alleles `2` and `1`, meaning they inherited microhaplotype variant 2 from one parent and variant 1 from the other, where each variant represents a distinct multi-SNP sequence pattern (e.g., variant 1 might be "AACG" while variant 2 is "ATCG" for a 4-SNP microhaplotype).

> It is important to note that, as we apply multiple criteria, generated microhaplotype segments will be significantly lower than fully haplotype blocks. The main advantage of this system is that we effectively use lower computational resources than haplotype-based genotypes or even biallelic SNPs. Therefore, a specialised coding matrix and model are strongly required to handle such genotype data, particularly to prevent multicollinearity while effectively utilising the multiallelic power.

### Usage Example

```bash
# Microhaplotype discovery and genotyping (LD-based, physical window)
./haplotype-hybrid \
      -i hap/chr* \
      -m map.txt \
      -o mh_info_ld_micro \
      --generate-genotypes mh_genotypes \
      -v \
      ld-haploblock micro --window-bp 125

# Haplotype discovery and genotyping (LD-based, best SNP selection)
./haplotype-hybrid \
      -i hap/chr* \
      -m map.txt \
      -o mh_info_ld_pure \
      --generate-genotypes hap_genotypes \
      -v \
      ld-haploblock pure --window 4

# Fixed physical distance blocks (100 kb)
./haplotype-hybrid \
      -i hap/chr* \
      -m map.txt \
      -o mh_info_fixedkb \
      --generate-genotypes fixedkb_genotypes \
      -v \
      fixed-kb --window-bp 100000

# Fixed SNP count per block
./haplotype-hybrid \
      -i hap/chr* \
      -m map.txt \
      -o mh_info_snpcount \
      --generate-genotypes snpcount_genotypes \
      -v \
      snp-count --window 4
```

### Key Parameters

- `-i, --input <INPUT>...`: Phased haplotype files, one per chromosome (required, e.g., `hap/chr*`)
- `-m, --map <MAP>`: SNP map file with columns: SNPID, Chr, Position (required)
- `--method <METHOD>`: Block definition method (default: `ld_haploblock`, possible values: `ld-haploblock`, `snp-count-simple`)
- `--haplotype-type <HAPLOTYPE_TYPE>`: For LD method - `pure` selects best N-SNP haplotype per block, `micro` splits blocks into physical windows (default: `pure`, possible values: `pure`, `micro`)
- `--d-prime-threshold <D_PRIME_THRESHOLD>`: D' threshold for consecutive SNP LD in haploblock discovery (default: `0.45`, Jonas 2017)
- `--window-bp <WINDOW_BP>`: Physical window size in base pairs for microhaplotype segments (default: `125`)
- `-w, --window <WINDOW>`: Haplotype size in SNPs for `pure` mode or window size for `snp-count-simple` (default: `4`)
- `--min-snps <MIN_SNPS>`: Minimum SNPs required per block (default: `2`)
- `--max-snps <MAX_SNPS>`: Maximum SNPs allowed per block (default: `4`)
- `--aft <AFT>`: Allele Frequency Threshold for Criterion-B score calculation (default: `0.08`)
- `--md <MD>`: Maximum Deviation parameter for Criterion-B weighting (default: `0.1`)
- `--min-ld <MIN_LD>`: Minimum mean LD (r²) threshold for block filtering (optional, e.g., `0.3`)
- `--no-dedup`: Skip deduplication of blocks with identical SNP ranges
- `--noheader`: Specify if haplotype input files have no header row
- `-o, --output <OUTPUT>`: Output directory for block definition files (default: `hap_info_microhap`)
- `--generate-genotypes <GENERATE_GENOTYPES>`: Generate haplotype genotype files with specified prefix (e.g., `mh_geno_ld`)
- `--missing <MISSING>`: Code for missing alleles in genotype output (e.g., `9`)
- `--missing-value <MISSING_VALUE>`: Replacement value for blocks containing missing alleles (default: `-9999`)
- `-v, --verbose`: Display detailed progress and diagnostic information
- `-h, --help`: Print help information

## Want to help us?

Contributions are welcome and very beneficial!
You can email me to improve the Rust implementation, add a new model, documentation, benchmarks, or bug reporting. I will appreciate, really

## License

GPL-3 License - see [LICENSE](LICENSE) file

Copyright (c) 2025 Agus Wibowo

## Contact

- **Email**: aguswibowo1698@gmail.com

## References

- Da, Y. Multi-allelic haplotype model based on genetic partition for genomic prediction and variance component estimation using SNP markers. [BMC Genet. 16, 144(2015)](https://doi.org/10.1186/s12863-015-0301-1)

- Jonas, D. et.al. Alternative haplotype construction methods for genomic evaluation. [Journal of Dairy Science. 99, 6(2016)](https://www.journalofdairyscience.org/article/S0022-0302(16)30076-5/fulltext)

- Prakapenka, D et.al. GVCHAP: A computing pipeline for genomic prediction and variance component estimation using haplotypes and SNP markers. [Front Genet. 11, 282(2020)](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00282/full)

## Development Team

**Lead Developer:** Agus Wibowo  
James Cook University

**Supervisors:**  
- Prof. Kyall Zenger
- Dr. Cecile Massault