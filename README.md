# Gene Taxonomic Prevalence Estimation Tool for Bacterial Taxonomic Groups

## Overview
This tool is designed for estimating the prevalence of specific genes in bacterial taxa, integrating NCBI genome retrieval, BLAST database construction, and automated query analysis.

## Features

### Genomic Data Acquisition

> ⚠️ **Note:** Users must choose between:
>
> - **Pre-built default database**: Use the included curated BLAST database of complete genomes from seven Enterobacteriaceae genus.
> - **Custom genome panel (`-d`)**: Create your own complete-genome BLAST database using any bacterial taxon.
>

* **Storage‑friendly metadata:** Because of the large storage requirements, genome sequence files used to build the BLAST database are **not** stored in this repository. Instead, `database/metadata` holds the assembly‑accession lists so users can re‑download any sequence on demand.  

* **Pre‑built complete‑genome database:** A BLAST database built from complete genomes of the default seven Enterobacteriaceae genus is provided for quick, high‑quality searches and hosted on USDA Ag Data Commons. The script `build_EB_complete_genomes_database.sh` is also provided for reproducing or updating the default database. 

* **Heavy mode (`‑H heavy` + `‑t <taxon_file>`):** Adds draft genomes to the analysis. Draft assemblies for each target taxon are downloaded with *ncbi‑genome‑download*, shuffled, and sampled in iterations (Cochran’s formula with finite‑population correction; ≤ 20 iterations) to capture diversity while keeping runtime reasonable.  

* **Custom genome panel (`‑d <download_file>`):** Enables users to work outside the default Enterobacteriaceae genus or combine targets from both within and outside this group. Provide a plain text file with one taxon per line (e.g., genus, species, or serotype), or pass a single genus/species/serotype directly. The pipeline will download the corresponding complete genomes, build a custom BLAST database in real time, and run either LIGHT or HEAVY mode against that database.

> **Note:** When using -d, you do not need to provide the -t flag for HEAVY mode—the custom panel already defines the target taxa.
> **System requirements** Building either the pre‑built archive or a large custom panel requires significant memory (>= 64 GB).

### Advanced Options

* **Output file naming (`-o`)**  
  By default, each run is written to a unique results directory based on the run timestamp and target gene coverage and identity (e.g. `20260125_110648_i95_c90`), and the summary file is named as `gene_summary.csv`. But user may optionally specify a custom output filename using -o: `-o salmonella_gene_summary.csv`.

* **Force rebuild of custom database (`-F true`)**  
  If a custom genome panel (`-d`) was previously built, the pipeline will skip rebuilding it unless forced. Use `-F true` to **delete and rebuild** all custom genome and BLAST database files from scratch.

* **Expand all species under genus (`--get-all-species`)**  
  When used with a genus in the `-d` file, this option will automatically retrieve **all species** under that genus using NCBI Taxonomy and classfiy these genomes into species level for each genus.

> These options are especially useful for power users working on broad taxonomic groups or managing reproducible BLAST database builds across runs.
---

### Genome files Organization for Pre-built Default Database

* Creates a structured directory hierarchy per genus → species → (for *Salmonella enterica*) subspecies → serotype. Unclassified genomes are placed in an `unclassified/` subfolder for each taxon.  

* Complete‑genome downloads live under `complete_genome/`; draft‑genome iterations are stored separately under `draft_genome/`.  
The script to download complete genomes and construct their corresponding BLAST DB is build_EB_complete_genomes_db.sh
<pre>
│   ├── Escherichia/
│   │   ├── unclassified/
│   │   ├── Escherichia_coli/
│   │   ├── Escherichia_fergusonii/
│   │   ├── Escherichia_albertii/
│   │   ├── ...
│   ├── Salmonella/
│   │   ├── unclassified/
│   │   ├── Salmonella_enterica/
│   │   │   ├── unclassified/
│   │   │   ├── enterica/
│   │   │   │   ├── unclassified/
│   │   │   │   ├── Typhimurium/
│   │   │   │   ├── Infantis/
│   │   │   │   ├── Newport/
│   │   │   │   ├── Heidelberg/
│   │   │   │   ├── ...
│   │   │   ├── salamae/
│   │   │   ├── arizonae/
│   │   │   ├── diarizonae/
│   │   │   ├── houtenae/
│   │   │   ├── indica/
│   │   ├── Salmonella_bongori/
│   ├── Shigella/
│   │   ├── unclassified/
│   │   ├── Shigella_flexneri/
│   │   ├── Shigella_sonnei/
│   │   ├── Shigella_boydii/
│   │   ├── ...
│   ├── Klebsiella/
│   │   ├── unclassified/
│   │   ├── Klebsiella_pneumoniae/
│   │   ├── Klebsiella_oxytoca/
│   │   ├── ...
│   ├── Enterobacter/
│   │   ├── unclassified/
│   │   ├── Enterobacter_cloacae/
│   │   ├── Enterobacter_hormaechei/
│   │   ├── ...
│   ├── Citrobacter/
│   │   ├── unclassified/
│   │   ├── Citrobacter_freundii/
│   │   ├── Citrobacter_koseri/
│   │   ├── ...
│   ├── Cronobacter/
│   │   ├── unclassified/
│   │   ├── Cronobacter_sakazakii/
│   │   ├── Cronobacter_malonaticus/
│   │   ├── ...
</pre>
---

### BLAST Query & Analysis

* **Light mode (default):** Runs BLAST only against the complete‑genome database (pre‑built or custom `‑d`). Fastest, highest assembly quality.  

* **Heavy mode:** Adds a draft‑genome database, dynamically constructed during runtime to improve prevalence estimates in taxonomic groups with sparse complete genomes. Requires both `‑H heavy` and a target taxon file via `‑t`. If custom target feature is initiated (`-d`), no `-t` target taxon file needs to be provided.  

* **Query genes:** Provide a multi‑FASTA file with one or many genes via `‑g`.  

* **Filtering:** Results are filtered by user‑defined minimum identity (`‑i`) and coverage (`‑c`) thresholds.

---

### Gene prevalence calculation

* **Light mode:** Prevalence is calculated from hits in complete genomes only.  

* **Heavy mode:** Prevalence combines complete‑genome hits and the averaged result of up to 20 draft‑genome iterations, improving robustness to sampling bias.

---

### SLURM Integration

GeTPrev is designed to run efficiently on high‑performance computing (HPC) systems managed by **SLURM** workload managers.

Built‑in SLURM support enables:

* **Automatic job submission**: GeTPrev can generate SLURM batch scripts internally when queues (`‑q`), resources (`‑m`, `‑r`), and node constraints are specified.

* **Flexible resource requests**: Users can control memory, CPU cores, wall‑time, and partition selection through command‑line options.

* **Robust parallelism**: Multiple `blastn` processes and download threads are automatically parallelized across available CPUs.

* **Queue compatibility**: Tested on partitions like `ceres`, `short`, and `long` on SCINet clusters (Ceres, Atlas).

> **Important:** It is recommended to run database‑building steps (`build_EB_complete_genomes_database.sh`, custom complete genome panels via `‑d`) on HPC nodes with ≥ 250 GB free storage.

#### Example SLURM submission (automatic through GeTPrev)

```bash
# Run default light mode using custom database on SCINet (partition ceres, 8 CPUs, 16 GB RAM, 4 hours max)
bash getprev.sh -g genes.fasta -d download_taxon.txt -q ceres -C 8 -m 16G -r 04:00:00
```

### SLURM Parameters Exposed in GeTPrev

| Flag | Meaning |
|:-----|:--------|
| `‑q` | Partition (queue) name (e.g., `ceres`, `short`, `long`) |
| `‑r` | Requested wall‑time (e.g., `04:00:00`) |
| `‑m` | Requested memory (e.g., `16G`) |
| `‑C` | Number of CPU cores |

#### Monitoring Jobs

Monitor your jobs using the following SLURM commands:

```bash
# View running and queued jobs
squeue -u $USER

# Check detailed information on a specific job
scontrol show job <jobID>
```
## Getting Started

This guide walks you through setting up everything needed to run the **GeTPrev** pipeline, including Conda setup, cloning the repository, and installing dependencies. If you already have Conda installed, you can skip to the [Installation](#installation) section.

---

## Conda Setup (If Not Installed)

If you do not already have Conda or Miniconda installed, follow the steps below to install **Miniconda** (recommended for a lightweight setup):

### 1. Download Miniconda

Go to the [official Miniconda installation page](https://docs.conda.io/en/latest/miniconda.html) and download the appropriate installer for your operating system.

Or use the terminal (for Linux/macOS):

```sh
# For Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# For macOS
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```

### 2. Run the Installer

Once downloaded, run the installer script:

```sh
bash Miniconda3-latest-Linux-x86_64.sh   # or the macOS version if applicable
```

### 3. Restart Shell & Initialize Conda

After installation, restart your terminal or manually initialize Conda by sourcing your shell configuration file:

```sh
# For Bash users
source ~/.bashrc

# For Zsh users
source ~/.zshrc
```

To verify installation:

```sh
conda --version
```

## Installation

To run this pipeline, set up a Conda environment with the required dependencies.
1. Clone the Repository
```sh
git clone --branch main --depth 1 https://github.com/Weifanwu66/GeTPrev.git getprev
cd getprev
```
2. Create and Activate the Conda Environment
The pipeline requires a Conda environment with all necessary dependencies. To create and activate it, run:
```sh
conda env create -f environment.yml
conda activate getprev
```
To verify the installation, check if all tools are installed:
```sh
conda list | grep -E "blast|ncbi-genome-download|entrez-direct"
```
If any package is missing, please install it manually:
```sh
conda install -c bioconda <package_name>
```

---

## Dependencies
1. ncbi-genome-download: Blin, K. (2023). ncbi-genome-download (0.3.3). Zenodo. https://doi.org/10.5281/zenodo.8192486
2. NCBI BLAST: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
3. entrez-direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/

-----
## Download the pre-built default database from Zenodo
The pre-built BLAST database archive is named blastdb_all.tar.gz. Download it from Zenodo, extract it, and place the contents directly under database/complete_blast_db/.

### HPC (computer node)
On many HPC systems, large external downloads are throttled on login nodes. For best speed, start an interactive session and download from a compute node, ideally using aria2c for parallel connections and resume support.
```bash
# Link to the database
ZENODO_URL="https://zenodo.org/records/17346156/files/blastdb_all.tar.gz?download=1"

# Create target folder
mkdir -p database/complete_blast_db

# Download
aria2c -x 8 -s 8 -k 10M -o blastdb_all.tar.gz "$ZENODO_URL"

# Extract into the expected folder
tar -xzf blastdb_all.tar.gz -C database/complete_blast_db

# Clean up
rm -f blastdb_all.tar.gz
```
If aria2c is not available, try:
```bash
curl -L -C - -o blastdb_all.tar.gz "$ZENODO_URL"
```

### Local machine (Linux/macOS)
```bash
ZENODO_URL="https://zenodo.org/records/17346156/files/blastdb_all.tar.gz?download=1"

# create target folder
mkdir -p database/complete_blast_db

# download
curl -L -o blastdb_all.tar.gz "$ZENODO_URL"

# extract directly into the expected folder
tar -xzf blastdb_all.tar.gz -C database/complete_blast_db

# clean up
rm -f blastdb_all.tar.gz
```

### Local machine (Windows)
Running the pipeline directly on Windows is not recommended, as many bioinformatics tools and dependencies are designed for Unix-like environments and may not be fully compatible or supported on Windows.
```powerShell
$ZENODO_URL = "https://zenodo.org/records/17346156/files/blastdb_all.tar.gz?download=1"

# create target folder
New-Item -ItemType Directory -Force -Path "database/complete_blast_db" | Out-Null

# download
Invoke-WebRequest -Uri $ZENODO_URL -OutFile "blastdb_all.tar.gz"

# extract (Windows 10+ has tar)
tar -xzf .\blastdb_all.tar.gz -C .\database\complete_blast_db

# clean up
Remove-Item .\blastdb_all.tar.gz
```

## Build default database locally using script `build_EB_complete_genomes_db.sh`

The database-building script can be executed **either locally (no job scheduler required)** or on an **HPC system with a job scheduler (e.g., SLURM)**.  
Scheduler-related parameters are **optional** and are **only used on systems with a job scheduler**.

---

### Option 1: Local execution (NO job scheduler required)

If you are running this script on a system **without a job scheduler** (e.g., a workstation or server) and your machine has sufficient available storage (≥128 GB), **no scheduler-related options are required**:

```bash
bash build_EB_complete_genomes_db.sh
```

### Option 2: Local execution (job scheduler required)

To download and build the default Enterobacteriaceae database on HPC systems with job schedulers, run:

```bash
ACCOUNT=my_account QUEUE=your_queue bash build_EB_complete_genomes_db.sh
```

You may override the defaults by setting environment variables when running the script. The default resources are: RAM=128 GB, CPUS=16, RUNTIME=12:00:00.

```bash
ACCOUNT=my_account QUEUE=your_queue RAM=256GB CPUS=24 RUNTIME=08:00:00 bash build_EB_complete_genomes_db.sh
```

## Example commands for analysis:

> ⚠️ **Note:** if running on HPC, please replace ceres with a valid partition in your environment. Please download the pre-built database from Ag Data Commons or download the database locally when working with the default database.

### 1. Run default light mode with 95% of identity and 90% of coverage (no target (-t) is defined, so the pipeline will loop through all taxonomic group available in pre-built database)

```bash
bash getprev.sh -g test/test_gene.fasta -i 95 -c 90 -q ceres -r 04:00:00 -m 16G -C 8
```

### 2. Run with a single taxon target in light mode

```bash
bash getprev.sh -g test/test_gene.fasta -t "Salmonella" -i 95 -c 90 -q ceres -r 04:00:00 -m 16G -C 8
```

### 3. Run heavy mode with multiple targets listed in a text file

```bash
bash getprev.sh -g test/test_gene.fasta -t test/test_taxon.txt -q ceres -r 08:00:00 -m 64G -C 16 -H heavy
```

### 4. Custom genome panel — light mode

```bash
bash getprev.sh -g test/test_gene.fasta -d test/download_taxon.txt -q ceres -r 08:00:00 -m 16G -C 8
```

### 5. Custom genome panel — heavy mode

```bash
bash getprev.sh -g test/test_gene.fasta -d test/download_taxon.txt -q ceres -r 12:00:00 -m 64G -C 16 -H heavy
```

### 6. Customize your output file name

```bash
bash getprev.sh -g test/test_gene.fasta -d test/download_taxon.txt -q ceres -r 08:00:00 -m 64G -C 16 -F true -o my_gene_summary.csv
```

### 7. Expand genus into species

```bash
bash getprev.sh -g test/test_gene.fasta -d test/download_taxon.txt -q ceres -r 08:00:00 -m 64G -C 16 --get-all-species
```

---

*Note:* 
- `-g` (gene FASTA) is always required.
- `-t` (species target file) must be provided in heavy mode for default Enterobacteriaceae database.
- `-d` (custom genus panel file) is used to download genera outside the default Enterobacteriaceae set.
- Ensure requested memory (`-m`) and runtime (`-r`) are appropriate for your HPC environment.
- Building custom databases with `-d` should be performed on nodes with ≥250 GB available disk space.
- Test input files (`test_gene.fasta`, `test_taxon.txt`, and `download_taxon.txt`) can be found in the `test/` directory.
---

## Output
Examples of the final output file:

**Light mode**
| Organism                            | Gene_ID                      | Min_percentage_of_coverage | Min_percentage_of_identity | Total_draft_genomes | Total_complete_genomes | Complete_genomes_with_target_genes | Percentage_with_target_genes_complete_genomes |
| ----------------------------------- | ---------------------------- | -------------------------- | -------------------------- | ------------------- | ---------------------- | ---------------------------------- | --------------------------------------------- |
| Salmonella enterica subsp. enterica | NC_003197.2:c3010816-3009905 | 90                         | 95                         | 235792              | 2307                   | 1924                               | 83.00%                                        |
| Salmonella enterica subsp. enterica | NC_003197.2:3019846-3021524  | 90                         | 95                         | 235792              | 2307                   | 2270                               | 98.00%                                        |
| Salmonella enterica subsp. enterica | NC_003197.2:c2925778-2923593 | 90                         | 95                         | 235792              | 2307                   | 2276                               | 98.00%                                        |

**Heavy mode**
| Organism          | Gene_ID                     | Min_percentage_of_coverage | Min_percentage_of_identity | Total_draft_genomes | Total_complete_genomes | Draft_genomes_sample_size | Number_of_iterations | Complete_genomes_with_target_genes | Draft_genomes_with_target_genes | Percentage_with_target_genes_complete_genomes | Percentage_with_target_genes_draft_genomes |
| ----------------- | --------------------------- | -------------------------- | -------------------------- | ------------------- | ---------------------- | ------------------------- | -------------------- | ---------------------------------- | ------------------------------- | --------------------------------------------- | ------------------------------------------ |
| Salmonella Uganda | NC_003197.2:1707344-1707789 | 80                         | 90                         | 1173                | 12                     | 290                       | 1                    | 12                                 | 290                             | 100.00%                                       | 100.00%                                    |
