# ---------------------------------------------------------------------
# SPECTRONAUT pipeline rules 
# ---------------------------------------------------------------------

conda: "envs/proteomics-spectronaut.yml"

import pandas as pd
from pathlib import Path
import glob, sys
import os

# ---------------------------------------------------------------------
# 0. Load metadata & Define global wildcard lists
# ---------------------------------------------------------------------
metadata = pd.read_csv(f"{config['project_root']}data_input/ms_raw_files.csv", sep=";")
metadata.columns = metadata.columns.str.strip()
if not {"file_name", "workflow", "cohort"}.issubset(metadata.columns):
    sys.exit("[metadata error] Missing required columns in ms_raw_files.csv")

metadata["file_base"] = metadata["file_name"].apply(lambda f: os.path.splitext(f)[0])

WORKFLOWS = sorted(metadata["workflow"].unique())
COHORTS   = sorted(metadata["cohort"].unique())
SAMPLES   = sorted(metadata["file_base"].unique())

# ---------------------------------------------------------------------
# Load activation key and schemas from config
# ---------------------------------------------------------------------
ACTIVATION_KEY = config.get("spectronaut_activation_key")
DIRECT_SCHEMA = config.get("spectronaut_direct_schema", "")
REPORT_SCHEMA = config.get("spectronaut_report_schema", "")      

if not ACTIVATION_KEY:
    sys.exit("[config error] Missing 'spectronaut_activation_key' in config.yml")

# ---------------------------------------------------------------------
# 1. Fetch FASTA file: defined in 00_fetch_fasta.smk
# ---------------------------------------------------------------------
FASTA = Path(FASTA)

# ---------------------------------------------------------------------
# CLI flag helpers
# ---------------------------------------------------------------------
def direct_schema_flag():
    return f"-s \"{DIRECT_SCHEMA}\"" if DIRECT_SCHEMA else ""

def report_schema_flag():
    return f"-rs {REPORT_SCHEMA}" if REPORT_SCHEMA else ""

def fasta_flag():
    return f"-fasta {FASTA}"

# ---------------------------------------------------------------------
# 2. Convert raw files → htrms (Spectronaut format)
# ---------------------------------------------------------------------
rule convert_raw_to_htrms:
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
        htrms=(
            f"{config['project_root']}data_input/converted_files/"
            + "{workflow}/{cohort}/{sample}.htrms"
        )
    log:
        f"{config['project_root']}logs/step3/convert_raw_to_htrms_{{workflow}}_{{cohort}}_{{sample}}.log"
    shell:
        f"""
        mkdir -p $(dirname {{output.htrms}})
        mkdir -p $(dirname {{log}})
        module load singularity
        module load spectronaut/20.0.250602.92449
        spectronaut -activate {ACTIVATION_KEY}
        spectronaut convert \
          -i {{input.raw}} \
          -o {{output.htrms}} \
          &>> {{log}}
        spectronaut -deactivate
        """


# ---------------------------------------------------------------------
# 3. Marker rule — ensure all conversion is complete
# ---------------------------------------------------------------------
rule convert_all:
    input:
        expand(
            f"{config['project_root']}data_input/converted_files/{{workflow}}/{{cohort}}/{{sample}}.htrms",
            workflow=WORKFLOWS,
            cohort=COHORTS,
            sample=SAMPLES
        )
    output:
        marker = f"{config['project_root']}data_output/S_all_converted.marker"
    log:
        f"{config['project_root']}logs/step4/convert_all.log"
    shell:
        """
        mkdir -p $(dirname {log})
        touch {output.marker} &>> {log}
        """

# ---------------------------------------------------------------------
# 4. Spectronaut analysis per workflow (all cohorts pooled)
# ---------------------------------------------------------------------
rule spectronaut_analysis_workflows:
    input:
        spectra_dir = lambda wc: f"{config['project_root']}data_input/converted_files/{wc.workflow}"
    output:
        report = f"{config['project_root']}data_output/{{workflow}}/{{workflow}}_Spectronaut_Report.tsv"
    log:
        f"{config['project_root']}logs/step5/spectronaut_analysis_workflows_{{workflow}}.log"
    shell:
        f"""
        mkdir -p {config['project_root']}data_output/{{wildcards.workflow}}
        mkdir -p $(dirname {{log}})
        module load spectronaut
        spectronaut -activate {ACTIVATION_KEY}
        spectronaut direct \
          -n {config['project_name']}
          {direct_schema_flag()} \
          -o {{output.report}} \
          {fasta_flag()} \
          -d {{input.spectra_dir}} \
          {report_schema_flag()} &>> {{log}}
        spectronaut -deactivate
        """

# ---------------------------------------------------------------------
# 5. Spectronaut analysis per workflow + cohort
# ---------------------------------------------------------------------
rule spectronaut_analysis_cohorts:
    input:
        spectra_dir = lambda wc: (
            f"{config['project_root']}data_input/converted_files/{wc.workflow}/{wc.cohort}"
        )
    output:
        report = f"{config['project_root']}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}_Spectronaut_Report.tsv"
    log:
        lambda wc: f"{config['project_root']}logs/step6/spectronaut_analysis_cohorts_{wc.workflow}_{wc.cohort}.log"
    threads: 16
    resources:
        mem_mb = 262144
    shell:
        f"""
        mkdir -p {config['project_root']}data_output/{{wildcards.workflow}}_{{wildcards.cohort}}
        mkdir -p $(dirname {{log}})
        module load spectronaut
        spectronaut -activate {ACTIVATION_KEY}
        spectronaut direct \
          -n {config['project_name']}
          {direct_schema_flag()} \
          -o {{output.report}} \
          {fasta_flag()} \
          -d {{input.spectra_dir}} \
          {report_schema_flag()} &>> {{log}}
        spectronaut -deactivate
        """

# ---------------------------------------------------------------------
# 6. Spectronaut analysis — ALL samples together
# ---------------------------------------------------------------------
rule spectronaut_analysis_all:
    input:
        spectra_dir = f"{config['project_root']}data_input/converted_files/"
    output:
        report = f"{config['project_root']}data_output/ALL/ALL_Spectronaut_Report.tsv"
    log:
        f"{config['project_root']}logs/step6/spectronaut_analysis_all.log"
    threads: 112
    resources:
        mem_mb = 262144
    shell:
        f"""
        mkdir -p {config['project_root']}data_output/ALL
        mkdir -p $(dirname {{log}})
        module load spectronaut
        spectronaut -activate {ACTIVATION_KEY}
        spectronaut direct \
          -n {config['project_name']}
          {direct_schema_flag()} \
          -o {{output.report}} \
          {fasta_flag()} \
          -d {{input.spectra_dir}} \
          {report_schema_flag()} &>> {{log}}
        spectronaut -deactivate
        """
# ---------------------------------------------------------------------
# 7. Final marker — based on run_mode in config
# ---------------------------------------------------------------------
RUN_MODE = config.get("run_mode", "all")

if RUN_MODE == "cohort":
    rule S_all:
        input:
            expand(
                f"{config['project_root']}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}_Spectronaut_Report.tsv",
                workflow=WORKFLOWS,
                cohort=COHORTS
            )
        output:
            marker = f"{config['project_root']}data_output/S_all_spectronaut_complete.marker"
        log:
            f"{config['project_root']}logs/step7/S_all_cohort.log"
        shell:
            "mkdir -p $(dirname {log}) && touch {output.marker} &>> {log}"

elif RUN_MODE == "workflow":
    rule S_all:
        input:
            expand(
                f"{config['project_root']}data_output/{{workflow}}/{{workflow}}_Spectronaut_Report.tsv",
                workflow=WORKFLOWS
            )
        output:
            marker = f"{config['project_root']}data_output/S_all_spectronaut_complete.marker"
        log:
            f"{config['project_root']}logs/step7/S_all_workflow.log"
        shell:
            "mkdir -p $(dirname {log}) && touch {output.marker} &>> {log}"

elif RUN_MODE == "all":
    rule S_all:
        input:
            f"{config['project_root']}data_output/ALL/ALL_Spectronaut_Report.tsv"
        output:
            marker = f"{config['project_root']}data_output/S_all_spectronaut_complete.marker"
        log:
            f"{config['project_root']}logs/step7/S_all_all.log"
        shell:
            "mkdir -p $(dirname {log}) && touch {output.marker} &>> {log}"