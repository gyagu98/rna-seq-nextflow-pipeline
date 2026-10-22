# Portable RNA-Seq Analysis Pipeline (Nextflow DSL2)

[![CI Pipeline Test](https://github.com/gyagu98/rna-seq-nextflow-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/gyagu98/rna-seq-nextflow-pipeline/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-ready-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A reproducible, portable RNA-seq processing and quantification workflow written in **Nextflow DSL2**. The pipeline performs quality control, transcript quantification via Salmon, differential expression aggregation, and comprehensive quality reporting via MultiQC.

---

## Workflow Overview

```
                      +-----------------------------+
                      |   samplesheet.csv (FASTQs)  |
                      +--------------+--------------+
                                     |
               +---------------------+---------------------+
               |                                           |
               v                                           v
     +-------------------+                       +-------------------+
     |   FastQC (QC)     |                       |  Salmon Indexing  |
     +---------+---------+                       +---------+---------+
               |                                           |
               |                                           v
               |                                 +-------------------+
               |                                 |   Salmon Quant    |
               |                                 +---------+---------+
               |                                           |
               |                       +-------------------+-------------------+
               |                       |                                       |
               v                       v                                       v
     +-----------------------------------+                   +-------------------+
     |        MultiQC Aggregate          |                   |  DESeq2 Analysis  |
     +-----------------------------------+                   +-------------------+
```

---

## Features

- **Nextflow DSL2 Modular Design:** Processes are decoupled into reusable components (`modules/fastqc`, `modules/salmon`, `modules/multiqc`, `modules/deseq2`).
- **Containerization Support:** Out-of-the-box support for **Docker**, **Singularity / Apptainer**, and **Conda** via BioContainers.
- **Automated CI/CD Testing:** GitHub Actions workflow executes dry-run validation and profile testing on every commit.
- **Zero-Download Test Fixtures:** Includes self-contained synthetic paired-end FASTQ reads and a miniature transcriptome for instant testing without multi-gigabyte downloads.
- **HPC & Cloud Ready:** Dynamic resource allocation with error-retry strategies in `conf/base.config`.

---

## Quickstart

### Prerequisites
- [Nextflow](https://www.nextflow.io/) `>= 22.10.1`
- Java 11+
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/) or [Conda](https://docs.conda.io/)

### 1. Run the Built-In Test Profile
Test the entire pipeline in seconds using the lightweight fixture dataset:

```bash
# Using Docker
nextflow run gyagu98/rna-seq-nextflow-pipeline -profile test,docker

# Or using Conda
nextflow run gyagu98/rna-seq-nextflow-pipeline -profile test,conda
```

### 2. Run on Your Own RNA-Seq Data

Prepare a `samplesheet.csv` with paired-end FASTQ paths:

```csv
sample,fastq_1,fastq_2,condition
control_rep1,/path/to/ctrl1_R1.fastq.gz,/path/to/ctrl1_R2.fastq.gz,control
control_rep2,/path/to/ctrl2_R1.fastq.gz,/path/to/ctrl2_R2.fastq.gz,control
treated_rep1,/path/to/treat1_R1.fastq.gz,/path/to/treat1_R2.fastq.gz,treated
treated_rep2,/path/to/treat2_R1.fastq.gz,/path/to/treat2_R2.fastq.gz,treated
```

Execute the pipeline:

```bash
nextflow run gyagu98/rna-seq-nextflow-pipeline \
    --input samplesheet.csv \
    --transcriptome /path/to/transcripts.fa \
    --outdir my_results \
    -profile docker
```

---

## Pipeline Parameters

| Parameter | Description | Default |
| :--- | :--- | :--- |
| `--input` | Path to `samplesheet.csv` | `null` |
| `--transcriptome` | Path to reference transcriptome FASTA | `null` |
| `--outdir` | Output results directory | `'results'` |
| `--salmon_kmer` | K-mer length for Salmon indexing | `21` |
| `-profile` | Execution profile (`docker`, `singularity`, `conda`, `test`) | `'standard'` |

---

## Output Structure

```
results/
├── fastqc/
│   ├── sample1_1_fastqc.html
│   └── sample1_1_fastqc.zip
├── salmon_index/
│   └── (Index binary files)
├── salmon_quant/
│   ├── sample1/quant.sf
│   └── sample2/quant.sf
├── deseq2/
│   └── differential_expression_summary.csv
└── multiqc/
    └── multiqc_report.html
```

---

## Repository Structure

```
├── .github/workflows/
│   └── ci.yml               # GitHub Actions CI workflow
├── bin/
│   └── deseq2_analysis.R    # Differential expression script
├── conf/
│   ├── base.config          # Computational resource allocations
│   └── test.config          # Minimal test fixture configuration
├── modules/
│   ├── deseq2/main.nf       # DESeq2 module
│   ├── fastqc/main.nf       # FastQC module
│   ├── multiqc/main.nf      # MultiQC module
│   └── salmon/
│       ├── index.nf         # Salmon index module
│       └── quant.nf         # Salmon quantification module
├── tests/
│   └── data/                # Lightweight test FASTQs and transcriptome
├── Dockerfile               # Pipeline Docker container definition
├── environment.yml          # Conda environment definition
├── main.nf                  # Main workflow execution script
└── nextflow.config          # Master Nextflow configuration
```

---

## Author & Citation

Developed by **Gretta Y. Rozenberg**  
Ph.D. Researcher in Computational Biology & Bioinformatics  
- **GitHub:** [@gyagu98](https://github.com/gyagu98)  
- **LinkedIn:** [Gretta Rozenberg](https://www.linkedin.com/in/gyagudayeva/)

If you find this pipeline helpful in your research, feel free to star the repository!
