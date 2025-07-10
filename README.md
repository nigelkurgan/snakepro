# Proteomics Quantification and Quality Control Pipeline

This **Snakemake** pipeline allows proteomics researchers to **quantify large-cohort DIA mass spectrometry data** and **automatically generate a quality control (QC) report in PDF format**. It integrates preprocessing, filtering, normalization, PCA, outlier detection, and contamination analysis into a reproducible and scalable workflow.

---

## Installation

Clone the repository into your server's directory:

```bash
git clone --recursive https://github.com/mafaldacustodio/Proteomics-QC-Pipeline.git
```

After cloning, configure your pipeline by editing the configuration file.

---

## Configuration

Edit the **config.yml** file with the appropriate paths and options:

```yaml
# Path to the root directory of the pipeline
project_root: '.../pipeline/'

# Path to the metadata CSV file 
metadata: ".../pipeline/data_input/metadata_file.csv"

# Path to the FASTA file containing the reference protein sequence database (used for spectral library generation)
fasta_file: '.../pipeline/data_input/fasta_file'

# Filtering settings for Quality Control
# Options:
#"all" (all samples),
#"group" (at least one group)
filtering_option: all
# Specify group for filtering when using filtering_option: group
group: None

#Filtering percentage - default 20%
filtering_percentage: 0.20

# Metadata column for grouping in missing value analysis
missing_values_grouping: 'condition'

# Metadata columns used for PCA visualization
pca_factors: ["plate", "condition", "time"]

# Enable or disable the contamination panel analysis (True or False)
contamination_panel: True

# Path to the contamination panel Excel file
contamination_panel_file: "/.../pipeline/data_input/contamination_panel.xlsx"

#Indicate the samples to exclude, in the format of a list ["sxxx", "sxxx"]
exclude_samples: []

# Batch correction settings
batch_correction: False
# Metadata column indicating the batch factor
batch_effect_column: "plate"
```

### Add files

Place the following files in the specified directories:

- **Raw data files** 
  ```bash
  pipeline/data_input/raw/
  ```
Use the following command:  
  ```bash
  scp -r {directory}/raw wqs176@esrumhead01fl.unicph.domain:/projects/cbmr_fpm_soup-AUDIT/data/pipeline/data_input/raw/
  ```
- **FASTA file and metadata**
  ```bash
  pipeline/data_input/
  ```



---

## Usage
We recommend the creation of a **tmux session** to run these commands.
```bash
tmux
```

To run the full workflow and generate a QC report, execute:

```bash
bash run_slurm.sh convert_all
bash run_slurm.sh A_all
bash run_slurm.sh summarize_qc
```
This will perform quantification, generate QC plots, and output a full PDF report.

---

### Post-Quality Control

After analysing the report, you may want to:

#### Apply batch correction
Edit the  **config.yml** file in the batch correction settings and run the following command:
```bash
bash run_slurm.sh summarize_qc
```

#### Remove identified outlier samples
Edit the  **config.yml** file by adding the samples you want to remove in the "exclude_samples" list and run the following command:
```bash
bash run_slurm.sh remove_all_samples 
bash run_slurm.sh summarize_qc
```

---

## Output

The pipeline generates:
- **Quality Control visualizations in png** 
- A **comprehensive individual PDF report**
- A **summary PDF report** 

To download the files use:
```bash
scp -r kuxxx@esrumhead01fl.unicph.domain:/{directory}/pipeline/data_output/{QC_remove_samples_and_batch_effect/}/ ~/{directory}/
```
---

## Overview

This pipeline provides an automated solution for **quantitative proteomics analysis** based on **Data-Independent Acquisition (DIA) mass spectrometry**. The workflow is designed to handle large cohorts and generate both **individual** and **summary quality control reports**.


<img width="848" alt="Screenshot 2025-04-20 at 21 30 15" src="https://github.com/user-attachments/assets/b0bc5efe-4f48-4738-87e5-2d450beba3fe" />



The process begins with two main inputs:
- A **FASTA file** containing the reference protein sequence database.
- **Raw DIA-MS files**, which are converted into the open-format `.mzML` files using a raw file converter (MSConvert).

These inputs are then used in **DIA-NN**, which performs:
1. **Spectral library generation**, based on the FASTA file.
2. **Quantification**, where the spectral library and `.mzML` files are used to identify and quantify peptides and proteins.

Once quantification is complete, the pipeline proceeds to the **Quality Control (QC)** step.

The QC step generates two types of reports:
- An **Individual report** for detailed insights of each workflow.
- A **Summary report** providing a summary of all workflows.

Following QC analysis, users are given two options:
- **Batch Correction**: If batch effects are identified, users can correct the batch effects by updating the configuration file.
- **Sample Removal**: If outlier samples are detected, users can remove them by adding their sample identifiers to the configuration file.

If no corrections or removals are needed, the workflow is finalized.

---

## Features

**Preprocessing**
   - Log 2 transformation
   - Removal of contaminants
   - Imputation
   - Normalization
   - Filtering

**Sample Distribution Plot**
   - Displays the intensity distribution of all samples
   - Highlights potential outliers

**Missing Values Heatmap**
   - Visualizes missing data across proteins and samples

**Missing Values Per Group**
   - Compares the amount of missing data per group

**Missing Proteins Distribution by Group**
   - Shows how missing protein counts vary across experimental groups

**Sample Correlation Heatmap**
   - Correlation matrix between samples based on protein intensities (Pearson Correlation)

**Hierarchical Clustering Dendrogram**
   - Clusters samples using a distance-based method (Ward algorithm)

**PCA Plots**
   - PCA colored by cluster
   - Shaped by by metadata factors (e.g., plate, condition, time)
   - Pontential Outliers highlighthed by Euclidian and Mahnolobis distances

**Enrichment of Principal Components**
   - Visualizes enrichment of Principal Components
    
**Boxplot of Intra and Inter Plate Distances**
   - Visualizes within-plate vs between-plate variation

**Density Plots of Sample Intensities**
   - Density distribution of intensities per sample
   - Potential outliers are identified

**Mean Intensity Density Plot**
   - Mean density of all samples with potential outliers samples overlaid for better analysis

**Intensity Histogram (Imputed vs Non-Imputed)**
   - Comparison of the intensity distributions before and after imputation

**Albumin Concentration Plot**
   - Visualizes albumin concentration per sample

**Protein Rank-Abundance Plot**
   - Displays proteins ranked by abundance (log-scale)

**Intra-individual CV Scatter Plots**
   - For all samples and by condition

**Inter-individual CV Scatter Plots**
   - For all samples and by condition

**Contamination Panel**
   - Barplot of contaminant proteins
   - Potential outlier samples highlighted 

**Outlier table**
   - An overview table of potential outliers and the different analysis they are highlighted


## Limitations

