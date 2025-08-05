# Proteomics Quantification and Quality Control Pipeline

This **Snakemake** pipeline enables proteomics researchers to **quantify large-cohort DIA mass spectrometry data** and **automatically generate a quality control (QC) report in PDF format**. It supports both **DIA-NN** and **Spectronaut** quantification strategies, and is designed to be modular, robust, and user-friendly—even for those new to Snakemake.

---

## Table of Contents

- [Installation](#installation)
- [Configuration](#configuration)
- [Input Data Preparation](#input-data-preparation)
- [Usage](#usage)
- [Pipeline Structure & Modes](#pipeline-structure--modes)
- [Output](#output)
- [Features](#features)
- [Troubleshooting](#troubleshooting)
- [Overview](#overview)

---

## Installation

Clone the repository into your server's directory:

```bash
git clone --recursive https://github.com/nigelkurgan/snakepro.git
```

Install Snakemake (version 9.9.0 or compatible):

```bash
conda create -n snakemake -c bioconda -c conda-forge snakemake
conda activate snakemake
```

---

## Configuration

Edit the **config/config.yml** file with the appropriate paths and options.  
**Example:**

```yaml
# Tool selection
quant_tool: diann          # "diann" or "spectronaut"
run_mode: cohort           # "cohort", "workflow", or "all"
run_qc: true               # enable QC report generation
sn_batch_correct: false    # Spectronaut-specific flag

# Paths
project_root: '/path/to/pipeline/'
metadata: 'data_input/ms_raw_files.csv'
fasta_file: 'data_input/uniprot_reference.fa'

# QC settings
filtering_option: all       # "all" or "group"
group: null                 # required if filtering_option: group
filtering_percentage: 0.20
missing_values_grouping: 'condition'
pca_factors: ["plate", "condition", "time"]
contamination_panel: true
contamination_panel_file: 'data_input/contamination_panel.xlsx'
exclude_samples: []         # e.g. ["sample_001", "sample_042"]
batch_correction: false
batch_effect_column: "plate"
```

**Key parameters:**
- `quant_tool`: Choose `"diann"` or `"spectronaut"` for quantification.
- `run_mode`: Controls how quantification is performed:
    - `"cohort"`: Per workflow and cohort (default).
    - `"workflow"`: Per workflow (all cohorts pooled).
    - `"all"`: All samples pooled together.
- `run_qc`: If `true`, generates QC reports.
- `exclude_samples`: List of samples to exclude from analysis.
- `batch_correction`: Enable/disable batch correction.

---

## Input Data Preparation

1. **Raw data files**  
   Place your raw files (e.g., Thermo `.raw`, Bruker `.d` folders) in:
   ```
   data_input/raw_files/
   ```
   Example upload command:
   ```bash
   scp -r {directory}/raw_files user@server:/path/to/pipeline/data_input/raw_files/
   ```

2. **FASTA file and metadata**  
   Place your reference FASTA and metadata CSV in:
   ```
   data_input/
   ```
   - Metadata file: `ms_raw_files.csv` (must include columns: `file_name`, `workflow`, `cohort`)
   - FASTA file: e.g., `uniprot_reference.fa`

3. **(Optional) Contamination panel**  
   If using contamination analysis, place the panel file in `data_input/`.

---

## Usage

We recommend running the pipeline in a **tmux** session for stability:

```bash
tmux
```

To run the full workflow and generate a QC report:

```bash
bash run_slurm.sh all
```

This will perform quantification, generate QC plots, and output a full PDF report.

### Running Individual Pipeline Steps

You can also run specific steps:

```bash
bash run_slurm.sh convert_all            # Data conversion only
bash run_slurm.sh A_all                  # Full quantification (DIANN or Spectronaut)
bash run_slurm.sh summarize_qc           # QC reports only
```

### Post-Quality Control Actions

**Apply batch correction:**  
Edit `config/config.yml` to enable batch correction and specify the batch column, then run:
```bash
bash run_slurm.sh summarize_qc
```

**Remove identified outlier samples:**  
Add sample IDs to `exclude_samples` in `config/config.yml`, then run:
```bash
bash run_slurm.sh remove_all_samples 
bash run_slurm.sh summarize_qc
```

---

## Pipeline Structure & Modes

The pipeline is **modular** and adapts to your configuration:

- **Quantification tool:**  
  Select between DIA-NN and Spectronaut via `quant_tool` in the config.

- **Run modes:**  
  - `"cohort"`: Quantifies each workflow/cohort combination separately.
  - `"workflow"`: Quantifies all cohorts pooled per workflow.
  - `"all"`: Quantifies all samples pooled together.

- **Automatic file handling:**  
  - Raw file types (e.g., `.raw`, `.d`) are inferred from your metadata.
  - Output files are always in the correct format for the selected tool.

- **Final marker files:**  
  - `data_output/A_all_diann_complete.marker` (for DIA-NN)
  - `data_output/S_all_spectronaut_complete.marker` (for Spectronaut)
  These indicate successful completion of the pipeline.

---

## Output

The pipeline generates:
- **Quantification results** (per run mode and tool)
- **Quality Control visualizations** (PNG)
- **Comprehensive individual PDF report**
- **Summary PDF report**

To download results:
```bash
scp -r user@server:/path/to/pipeline/data_output/ ~/local_directory/
```

---

## Features

**Preprocessing**
   - Log2 transformation
   - Removal of contaminants
   - Imputation
   - Normalization
   - Filtering

**Sample Distribution Plot**
   - Intensity distribution of all samples
   - Outlier highlighting

**Missing Values Heatmap**
   - Visualizes missing data across proteins and samples

**Missing Values Per Group**
   - Compares missing data per group

**Missing Proteins Distribution by Group**
   - Shows missing protein counts by group

**Sample Correlation Heatmap**
   - Pearson correlation matrix between samples

**Hierarchical Clustering Dendrogram**
   - Clusters samples using Ward algorithm

**PCA Plots**
   - Colored by cluster and metadata factors (e.g., plate, condition, time)
   - Outlier detection

**Enrichment of Principal Components**
   - Visualizes enrichment of principal components

**Boxplot of Intra and Inter Plate Distances**
   - Within-plate vs between-plate variation

**Density Plots of Sample Intensities**
   - Density distribution per sample, outlier identification

**Mean Intensity Density Plot**
   - Mean density with outliers overlaid

**Intensity Histogram (Imputed vs Non-Imputed)**
   - Compare before/after imputation

**Albumin Concentration Plot**
   - Albumin per sample

**Protein Rank-Abundance Plot**
   - Proteins ranked by abundance (log-scale)

**Intra/Inter-individual CV Scatter Plots**
   - For all samples and by condition

**Contamination Panel**
   - Barplot of contaminant proteins, outlier highlighting

**Outlier Table**
   - Overview of potential outliers and detection methods

---

## Troubleshooting

- **Missing required columns in metadata:**  
  Ensure your `ms_raw_files.csv` contains at least `file_name`, `workflow`, and `cohort`.

- **No FASTA file found:**  
  Place a single FASTA file in `data_input/`.

- **Pipeline fails on log file creation:**  
  All rules create log directories as needed; if you see errors, check directory permissions.

- **Not sure which run mode to use?**  
  - Use `"cohort"` for most studies (default).
  - Use `"workflow"` if you want to pool all cohorts per workflow.
  - Use `"all"` to pool all samples together.

---

## Overview

This pipeline provides an automated solution for **quantitative proteomics analysis** based on **Data-Independent Acquisition (DIA) mass spectrometry**. The workflow is designed to handle large cohorts and generate both **individual** and **summary quality control reports**.

<img width="848" alt="Pipeline Overview" src="https://github.com/user-attachments/assets/b0bc5efe-4f48-4738-87e5-2d450beba3fe" />

---

**For further help, open an issue on the [GitHub repository](https://github.com/nigelkurgan/snakepro/issues) or contact the pipeline maintainer.**