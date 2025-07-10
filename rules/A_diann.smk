# ---------------------------------------------------------------------
# DIANN pipeline rules
# ---------------------------------------------------------------------

conda: "envs/proteomics-diann.yml"

import pandas as pd
from pathlib import Path
import glob, sys

# ---------------------------------------------------------------------
# 1. Load metadata & Define global wildcard lists
# ---------------------------------------------------------------------
metadata = pd.read_csv(f"{config['project_root']}data_input/ms_raw_files.csv", sep=";")
metadata.columns = metadata.columns.str.strip()
if not {"file_name", "workflow", "cohort"}.issubset(metadata.columns):
    sys.exit("[metadata error] Missing required columns in ms_raw_files.csv")

metadata["file_base"] = metadata["file_name"].str.replace(r"\.raw$", "", regex=True)

WORKFLOWS = sorted(metadata["workflow"].unique())
COHORTS   = sorted(metadata["cohort"].unique())
SAMPLES   = sorted(metadata["file_base"].unique())

# ---------------------------------------------------------------------
# 2. Generate predicted spectral library
# ---------------------------------------------------------------------
rule generate_spectral_library:
    input:
        fasta = f"{config['project_root']}data_input/uniprotkb_reviewed_true_AND_model_organ_2025_02_10.fasta"
    output:
        library = f"{config['project_root']}data_output/library.predicted.speclib"
    threads: 16
    shell:
        """
        module load diann/2.0.1 gcc
        diann \
          --fasta {input.fasta} \
          --predictor \
          --fasta-search \
          --gen-spec-lib \
          --out-lib {output.library} \
          --threads {threads}
        """

# ---------------------------------------------------------------------
# 3. Convert raw files → mzML
# ---------------------------------------------------------------------
rule convert_raw_to_mzml:
    input:
        raw=lambda wc: (
            f"{config['project_root']}data_input/raw_files/" +
            metadata.loc[
                (metadata.file_base == wc.sample) &
                (metadata.workflow == wc.workflow) &
                (metadata.cohort == wc.cohort),
                "file_name"
            ].iloc[0]
        )
    output:
        mzml = (
            f"{config['project_root']}data_input/converted_files/" +
            "{workflow}/{cohort}/{sample}.mzML"
        )
    threads: 16
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p $(dirname {output.mzml})
        module load msconvert/20250218_1
        msconvert \
          --64 \
          --zlib \
          --filter "peakPicking" \
          --filter "zeroSamples removeExtra 1-" \
          --outdir $(dirname {output.mzml}) \
          {input.raw}
        """

# ---------------------------------------------------------------------
# 4. Marker rule — ensure all mzML conversion is complete
# ---------------------------------------------------------------------
rule convert_all:
    input:
        expand(
            f"{config['project_root']}data_input/converted_files/{{workflow}}/{{cohort}}/{{sample}}.mzML",
            workflow=WORKFLOWS,
            cohort=COHORTS,
            sample=SAMPLES
        )
    output:
        marker = f"{config['project_root']}data_output/A_all_converted.marker"
    shell:
        "touch {output.marker}"

# ---------------------------------------------------------------------
# 5. DIANN analysis per workflow (all cohorts pooled)
# ---------------------------------------------------------------------
rule diann_analysis_workflows:
    input:
        library = f"{config['project_root']}data_output/library.predicted.speclib",
        raw_data_dir = lambda wc: f"{config['project_root']}data_input/converted_files/{wc.workflow}"
    output:
        stats = f"{config['project_root']}data_output/{{workflow}}/{{workflow}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p {config['project_root']}data_output/{wildcards.workflow}
        module load diann/2.0.1 gcc
        diann \
          --lib {input.library} \
          --dir-all {input.raw_data_dir} \
          --threads {threads} \
          --out {output.stats} \
          --qvalue 0.01 \
          --matrix-qvalue 0.01 \
          --mass-acc-ms1 10 \
          --mass-acc 10 \
          --matrices \
          --missed-cleavages 1
        """

# ---------------------------------------------------------------------
# 6. DIANN analysis per workflow + cohort
# ---------------------------------------------------------------------
rule diann_analysis_cohorts:
    input:
        library = f"{config['project_root']}data_output/library.predicted.speclib",
        raw_data_dir = lambda wc: (
            f"{config['project_root']}data_input/converted_files/{wc.workflow}/{wc.cohort}"
        )
    output:
        stats = f"{config['project_root']}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p {config['project_root']}data_output/{wildcards.workflow}_{wildcards.cohort}
        module load diann/2.0.1 gcc
        diann \
          --lib {input.library} \
          --dir {input.raw_data_dir} \
          --threads {threads} \
          --out {output.stats} \
          --qvalue 0.01 \
          --matrix-qvalue 0.01 \
          --mass-acc-ms1 10 \
          --mass-acc 10 \
          --matrices \
          --missed-cleavages 1
        """

# ---------------------------------------------------------------------
# 7. Final marker — all DIANN results complete
# ---------------------------------------------------------------------
rule A_all:
    input:
        expand(
            f"{config['project_root']}data_output/{{workflow}}/{{workflow}}.stats.tsv",
            workflow=WORKFLOWS
        ),
        expand(
            f"{config['project_root']}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",
            workflow=WORKFLOWS, cohort=COHORTS
        )
    output:
        marker = f"{config['project_root']}data_output/A_all_diann_complete.marker"
    shell:
        "touch {output.marker}"
