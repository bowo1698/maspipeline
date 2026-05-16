<a id="readme-top"></a>

<h1 align="center">maspipeline</h1>

<p align="center"><em>Preprocessing pipeline for genomic prediction using microhaplotype markers</em></p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="GPL v3"></a>
  <a href="https://www.rust-lang.org/"><img src="https://img.shields.io/badge/Rust-000000?logo=rust&logoColor=white" alt="Rust"></a>
  <a href="https://bowo1698.github.io/masgenomics-docs/"><img src="https://img.shields.io/badge/docs-masgenomics--docs-success" alt="Docs"></a>
</p>

---

## Overview

`maspipeline` adapts and optimises the GVCHAP pipeline ([Prakapenka et al. 2020](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00282/full)) for microhaplotype-based genomic prediction. It runs in four stages — phasing → haplotype conversion → microhaplotype block discovery → genotyping. All stages are implemented in Rust for efficient memory use and parallelism across chromosomes.

---

## Features

- **`convert-to-vcf`** — convert CSV genotype data (0/1/2 coding) into a VCF file compatible with Beagle
- **`convert-from-vcf`** — extract per-chromosome phased haplotype files and a SNP map from a phased VCF
- **`haplotype-hybrid`** — discover microhaplotype segments (LD-based haploblock or fixed-window) and emit the final genotype matrix
- Phasing via Beagle (parallel by chromosome)
- LD-based haploblock detection with Criterion-B scoring ([Jonas et al. 2016](https://www.journalofdairyscience.org/article/S0022-0302(16)30076-5/fulltext))
- Parallelised Rust backend for large datasets

---

## Part of the masgenomics suite

`maspipeline` is one of three packages for end-to-end genomic prediction:

- **maspipeline** *(this repo)* — preprocessing (phasing → haploblock discovery → microhaplotype genotyping)
- **[masreml](https://github.com/bowo1698/masreml)** — REML-BLUP, GWAS (EMMAX), GWABLUP
- **[masbayes](https://github.com/bowo1698/masbayes)** — Bayesian genomic prediction (BayesA, BayesR)

Full documentation, tutorials, theory, and reference: **<https://bowo1698.github.io/masgenomics-docs/>**

---

## Installation

### Build from source (recommended)

Install Rust via [rustup](https://rustup.rs/), then:

```bash
git clone https://github.com/bowo1698/maspipeline.git
cd maspipeline
cargo build --release
# Binaries: ./target/release/{convert-to-vcf, convert-from-vcf, haplotype-hybrid}
```

### Pre-compiled binaries

Download platform-specific binaries (Linux x86_64, macOS Intel / Apple Silicon, Windows x86_64) from the [Releases](https://github.com/bowo1698/maspipeline/releases/tag/v1) page.

Detailed installation, PATH setup, and platform-specific notes:
[masgenomics-docs › maspipeline](https://bowo1698.github.io/masgenomics-docs/reference/maspipeline/).

---

## Quick start

The pipeline runs in four stages. A minimal end-to-end sketch:

```bash
# 1. (Optional) Convert CSV genotypes to VCF
./convert-to-vcf -i genotypes.csv -m map.txt -o genotypes.vcf \
    --split-by-chr --parallel --ncores 4 --merge-after-split

# 2. Phase per chromosome with Beagle (see docs for parallel script)

# 3. Convert phased VCF into per-chromosome haploid files
./convert-from-vcf -i genotypes_phased.vcf.gz --map map.txt \
    --genofolder geno --hapfolder hap

# 4. Discover microhaplotype segments and genotype
./haplotype-hybrid --ncores 4 -i hap/chr* -m map.txt -o mh_info_ld_micro \
    --generate-genotypes mh_genotypes -v \
    ld-haploblock micro --window-bp 125
```

The microhaplotype genotype matrix produced by Stage 4 is consumed by `masreml::build_G_mh()` and `masbayes::construct_wah_matrix()`.

Full stage-by-stage walkthroughs (data preparation, phasing, conversion, microhaplotype discovery) with parameter references:
[masgenomics-docs › maspipeline](https://bowo1698.github.io/masgenomics-docs/reference/maspipeline/).

---

## Citation

```bibtex
@software{maspipeline,
  author = {Wibowo, Agus},
  title  = {maspipeline: Preprocessing pipeline for genomic prediction using microhaplotype markers},
  url    = {https://github.com/bowo1698/maspipeline},
  year   = {2026}
}
```

---

## Development Team

**Lead Developer**

- Agus Wibowo — James Cook University

**Supervisors**

- Prof. Kyall Zenger
- Dr. Cecile Massault
- Dr. Dave Jones

---

## License

[GPL-3](LICENSE) © 2025 Agus Wibowo · Contact: aguswibowo1698@gmail.com

<p align="right"><a href="#readme-top">↑ back to top</a></p>
