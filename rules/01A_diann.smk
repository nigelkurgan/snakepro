# ---------------------------------------------------------------------
# DIANN pipeline rules
# ---------------------------------------------------------------------

conda: "envs/proteomics-diann.yml"

import pandas as pd
from pathlib import Path
import glob, sys
import os
import re

# ---------------------------------------------------------------------
# 0. Load metadata & Define global wildcard lists
# ---------------------------------------------------------------------
# Load metadata
metadata = pd.read_csv(f"{config['project_root']}data_input/ms_raw_files.csv", sep=";")
metadata.columns = metadata.columns.str.strip()

# Ensure required columns
if not {"file_name", "workflow", "cohort"}.issubset(metadata.columns):
    sys.exit("[metadata error] Missing required columns in ms_raw_files.csv")

# Strip known raw file extensions to get base name
def strip_extension(f):
    return re.sub(r'\.(raw|d)$', '', f, flags=re.IGNORECASE)

metadata["file_base"] = metadata["file_name"].apply(strip_extension)

# global wildcards
WORKFLOWS = sorted(metadata["workflow"].unique())
COHORTS   = sorted(metadata["cohort"].unique())
SAMPLES   = sorted(metadata["file_base"].unique())

# ---------------------------------------------------------------------
# 1. Fetch FASTA file: defined in 00_fetch_fasta.smk
# ---------------------------------------------------------------------
FASTA = Path(FASTA_OUT)

# ---------------------------------------------------------------------
# 2. Generate predicted spectral library
# ---------------------------------------------------------------------
rule generate_spectral_library:
    input:
        fasta=lambda wc: FASTA,
    output:
        library = f"{config['project_root']}data_output/library.predicted.speclib"
    threads: 16
    log:
        f"{config['project_root']}logs/step2/generate_spectral_library.log"
    shell:
        """
        mkdir -p $(dirname {log})
        module load diann/2.0.1 gcc
        diann \
          --fasta {input.fasta} \
          --predictor \
          --fasta-search \
          --gen-spec-lib \
          --out-lib {output.library} \
          --threads {threads} \
          &>> {log}
        """

# ---------------------------------------------------------------------
# 3. Convert raw files → mzML
# ---------------------------------------------------------------------
rule convert_raw_to_mzml:
    input:
        raw=lambda wc: str(
            Path(config['project_root']) / "data_input" / "raw_files" /
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
    log:
        lambda wc: f"{config['project_root']}logs/step3/convert_raw_to_mzml_{wc.workflow}_{wc.cohort}_{wc.sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.mzml})
        mkdir -p $(dirname {log})
        module load msconvert/20250218_1

        msconvert \
          --64 \
          --zlib \
          --filter "peakPicking" \
          --filter "zeroSamples removeExtra 1-" \
          --outdir $(dirname {output.mzml}) \
          {input.raw} \
          &>> {log}
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
    log:
        f"{config['project_root']}logs/step4/convert_all.log"
    shell:
        """
        mkdir -p $(dirname {log})
        touch {output.marker} &>> {log}
        """

# ---------------------------------------------------------------------
# 5. DIANN analysis per workflow (all cohorts pooled) 
### Use this when you have different sample types = workflow but want different groups/cohorts together
# ---------------------------------------------------------------------
rule diann_analysis_workflows:
    input:
        library = f"{config['project_root']}data_output/library.predicted.speclib",
        raw_data_dir = lambda wc: f"{config['project_root']}data_input/converted_files/{wc.workflow}"
    output:
        stats = f"{config['project_root']}data_output/{{workflow}}/{{workflow}}.stats.tsv"
    threads: 56
    resources:
        mem_mb = 262144
    log:
        lambda wc: f"{config['project_root']}logs/step5/diann_analysis_workflows_{wc.workflow}.log"
    shell:
        """
        mkdir -p {config['project_root']}data_output/{wildcards.workflow}
        mkdir -p $(dirname {log})
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
          --missed-cleavages 1 \
          &>> {log}
        """

# ---------------------------------------------------------------------
# 6. DIANN analysis per worklow and per cohort
# ---------------------------------------------------------------------
rule diann_analysis_cohorts:
    input:
        library = f"{config['project_root']}data_output/library.predicted.speclib",
        raw_data_dir = lambda wc: (
            f"{config['project_root']}data_input/converted_files/{wc.workflow}/{wc.cohort}"
        )
    output:
        stats = f"{config['project_root']}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv"
    log:
        lambda wc: f"{config['project_root']}logs/step6/diann_analysis_cohorts_{wc.workflow}_{wc.cohort}.log"
    threads: 56
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p {config['project_root']}data_output/{wildcards.workflow}_{wildcards.cohort}
        mkdir -p $(dirname {log})
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
          --missed-cleavages 1 \
          &>> {log}
        """

# ---------------------------------------------------------------------
# 7. DIANN analysis — ALL samples together (no cohort/workflow split)
# ---------------------------------------------------------------------
rule diann_analysis_all:
    input:
        library = f"{config['project_root']}data_output/library.predicted.speclib",
        raw_data_dir = f"{config['project_root']}data_input/converted_files/"
    output:
        stats = f"{config['project_root']}data_output/ALL/ALL.stats.tsv"
    log:
        f"{config['project_root']}logs/step6/diann_analysis_all.log"
    threads: 56
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p {config['project_root']}data_output/ALL
        mkdir -p $(dirname {log})
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
          --missed-cleavages 1 \
          &>> {log}
        """

# ---------------------------------------------------------------------
# 8. Final marker — based on run_mode in config
# ---------------------------------------------------------------------
RUN_MODE = config.get("run_mode", "all")

if RUN_MODE == "cohort":
    rule diann_all:
        input:
            expand(
                f"{config['project_root']}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",
                workflow=WORKFLOWS,
                cohort=COHORTS
            )
        output:
            marker = f"{config['project_root']}data_output/A_all_diann_complete.marker"
        log:
            f"{config['project_root']}logs/step7/A_all_cohort.log"
        shell:
            "mkdir -p $(dirname {log}) && touch {output.marker} &>> {log}"

elif RUN_MODE == "workflow":
    rule diann_all:
        input:
            expand(
                f"{config['project_root']}data_output/{{workflow}}/{{workflow}}.stats.tsv",
                workflow=WORKFLOWS
            )
        output:
            marker = f"{config['project_root']}data_output/A_all_diann_complete.marker"
        log:
            f"{config['project_root']}logs/step7/A_all_workflow.log"
        shell:
            "mkdir -p $(dirname {log}) && touch {output.marker} &>> {log}"

elif RUN_MODE == "all":
    rule diann_all:
        input:
            f"{config['project_root']}data_output/ALL/ALL.stats.tsv"
        output:
            marker = f"{config['project_root']}data_output/A_all_diann_complete.marker"
        log:
            f"{config['project_root']}logs/step7/A_all_all.log"
        shell:
            "mkdir -p $(dirname {log}) && touch {output.marker} &>> {log}"
