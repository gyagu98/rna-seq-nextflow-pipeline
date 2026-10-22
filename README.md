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

## Setup & Installation

### 1. Prerequisites
Ensure you have the following installed on your machine or HPC environment:
- **Java 11 or later**:
  ```bash
  java -version
  ```
- **Nextflow** (>= 22.10.1):
  ```bash
  curl -s https://get.nextflow.io | bash
  chmod +x nextflow
  sudo mv nextflow /usr/local/bin/
  # Verify installation
  nextflow info
  ```
- **Container Engine or Package Manager** (at least one):
  - [Docker Desktop](https://www.docker.com/) (recommended for local development / macOS / Windows)
  - [Singularity / Apptainer](https://apptainer.org/) (standard for shared HPC clusters)
  - [Conda / Mamba](https://github.com/conda-forge/miniforge) (lightweight environments)

### 2. Clone the Repository
```bash
git clone https://github.com/gyagu98/rna-seq-nextflow-pipeline.git
cd rna-seq-nextflow-pipeline
```
*(Alternatively, Nextflow can pull and run the repository directly from GitHub without cloning).*

---

## Test-Profile Usage

The pipeline includes a self-contained test profile with synthetic paired-end FASTQ reads and a miniature transcriptome in `tests/data/`. This enables zero-download validation of all pipeline components.

### A. Preview / Syntax Validation
Verify pipeline syntax and channel resolution without launching tasks:
```bash
nextflow run main.nf -profile test -preview
```

### B. Stub Run (Instant DAG & Channel Verification)
Validate task dependencies, command stubs, and file staging in seconds without container downloads or heavy computation:
```bash
nextflow run main.nf -profile test -stub-run
```

### C. End-to-End Test Execution
Execute all processes (FastQC, Salmon Index, Salmon Quant, MultiQC, and DESeq2) with real containers:

```bash
# With Docker
nextflow run main.nf -profile test,docker

# With Singularity / Apptainer (HPC)
nextflow run main.nf -profile test,singularity

# With Conda / Mamba
nextflow run main.nf -profile test,conda

# Run directly from GitHub without local cloning
nextflow run gyagu98/rna-seq-nextflow-pipeline -profile test,docker
```

---

## Running on Your Own RNA-Seq Data

### 1. Prepare your `samplesheet.csv`
Create a CSV file referencing your paired-end FASTQ reads and experimental conditions:

```csv
sample,fastq_1,fastq_2,condition
control_rep1,/data/fastqs/ctrl_1_R1.fastq.gz,/data/fastqs/ctrl_1_R2.fastq.gz,control
control_rep2,/data/fastqs/ctrl_2_R1.fastq.gz,/data/fastqs/ctrl_2_R2.fastq.gz,control
treated_rep1,/data/fastqs/treat_1_R1.fastq.gz,/data/fastqs/treat_1_R2.fastq.gz,treated
treated_rep2,/data/fastqs/treat_2_R1.fastq.gz,/data/fastqs/treat_2_R2.fastq.gz,treated
```

### 2. Standard Production Execution
Run the workflow with a reference transcriptome (FASTA) and output directory:

```bash
nextflow run main.nf \
    --input /path/to/samplesheet.csv \
    --transcriptome /path/to/reference_transcripts.fa \
    --outdir ./results \
    -profile docker
```

### 3. Resuming an Interrupted Run
Nextflow uses deterministic task hashes to cache completed steps. If a run is interrupted or inputs change downstream, resume without recomputing finished tasks:

```bash
nextflow run main.nf \
    --input /path/to/samplesheet.csv \
    --transcriptome /path/to/reference_transcripts.fa \
    --outdir ./results \
    -profile docker \
    -resume
```

### 4. HPC Cluster Submission (SLURM / SGE)
When deploying on high-performance computing clusters with Singularity/Apptainer:

```bash
nextflow run main.nf \
    --input /path/to/samplesheet.csv \
    --transcriptome /path/to/reference_transcripts.fa \
    --outdir /scratch/project/results \
    -profile singularity \
    -with-report execution_report.html \
    -with-timeline execution_timeline.html \
    -with-trace execution_trace.txt
```

### 5. Custom Parameter Overrides
You can tune index k-mer size or memory/CPU directives via the command line or custom config:

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --transcriptome transcripts.fa \
    --salmon_kmer 31 \
    --outdir run_k31_results \
    -profile docker
```

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
