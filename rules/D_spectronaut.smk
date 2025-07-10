# Spectronaut pipeline
conda: "envs/proteomics-spectronaut.yml"

import pandas as pd
from pathlib import Path
import glob, sys

# ---------------------------------------------------------------------
# 1. Load and validate metadata
# ---------------------------------------------------------------------
metadata = pd.read_csv("data_input/ms_raw_files.csv", sep=";")
metadata.columns = metadata.columns.str.strip()
required = {"file_name", "workflow", "cohort"}
if not required.issubset(metadata.columns):
    sys.exit(f"[metadata error] Missing columns: {required}")

# Split into base name and extension
metadata[["file_base", "ext"]] = metadata["file_name"].str.rsplit(".", 1, expand=True)

# ---------------------------------------------------------------------
# 2. Define global lists for workflows and cohorts
# ---------------------------------------------------------------------
WORKFLOWS = sorted(metadata["workflow"].unique())
COHORTS = sorted(metadata["cohort"].unique())
SAMPLES = sorted(metadata["file_base"].unique())
EXTS = sorted(metadata["ext"].unique())

# ---------------------------------------------------------------------
# 3. Locate FASTA file
# ---------------------------------------------------------------------
fasta_files = glob.glob("data_input/*.fa*")
if len(fasta_files) != 1:
    sys.exit("[FASTA error] Expect exactly one FASTA in data_input/")
FASTA = Path(fasta_files[0]).name  # used in shell commands

# ---------------------------------------------------------------------
# 4. Convert raw → .htrms
# ---------------------------------------------------------------------
rule convert_to_htrms:
    input:
        raw=lambda wc: metadata.loc[
            (metadata.file_base == wc.sample) &
            (metadata.workflow == wc.workflow) &
            (metadata.cohort == wc.cohort),
            "file_name"
        ].iloc[0]
    output:
        htrms="data_input/converted_files/{workflow}/{cohort}/{sample}.htrms"
    threads: 4
    resources:
        mem_mb=64000
    shell:
        r"""
        mkdir -p data_input/converted_files/{wildcards.workflow}/{wildcards.cohort}
        module load --auto spectronaut/20.0.250602.92449
        spectronaut -convert \
            -i {input.raw} \
            -o {output.htrms} \
            -nogui
        """

# ---------------------------------------------------------------------
# 5. Spectronaut DirectDIA per workflow (pooled cohorts)
# ---------------------------------------------------------------------
rule spectronaut_workflow:
    input:
        json="data_input/spectronaut_schema.json",
        fasta=lambda wc: FASTA,
        htrms=expand(
            "data_input/converted_files/{workflow}/{cohort}/{sample}.htrms",
            workflow="{workflow}",
            cohort=COHORTS,
            sample=SAMPLES
        )
    output:
        report="data_output/{workflow}/{workflow}_Spectronaut_Report.tsv"
    threads: 16
    resources:
        mem_mb=262144
    shell:
        r"""
        mkdir -p data_output/{wildcards.workflow}
        module load --auto spectronaut/20.0.250602.92449
        spectronaut direct \
            -s {input.json} \
            -n {wildcards.workflow} \
            -o {output.report:r} \
            -fasta data_input/{FASTA} \
            -d data_input/converted_files/{wildcards.workflow}
        """

# ---------------------------------------------------------------------
# 6. Spectronaut DirectDIA per workflow + cohort
# ---------------------------------------------------------------------
rule spectronaut_cohort:
    input:
        json="data_input/spectronaut_schema.json",
        fasta=lambda wc: FASTA,
        htrms=expand(
            "data_input/converted_files/{workflow}/{cohort}/{sample}.htrms",
            workflow="{workflow}",
            cohort="{cohort}",
            sample=SAMPLES
        )
    output:
        report="data_output/{workflow}_{cohort}/{workflow}_{cohort}_Spectronaut_Report.tsv"
    threads: 16
    resources:
        mem_mb=262144
    shell:
        r"""
        mkdir -p data_output/{wildcards.workflow}_{wildcards.cohort}
        module load --auto spectronaut/20.0.250602.92449
        spectronaut direct \
            -s {input.json} \
            -n {wildcards.workflow}_{wildcards.cohort} \
            -o {output.report:r} \
            -fasta data_input/{FASTA} \
            -d data_input/converted_files/{wildcards.workflow}/{wildcards.cohort}
        """

# ---------------------------------------------------------------------
# 7. Spectronaut DirectDIA across all data
# ---------------------------------------------------------------------
rule spectronaut_all:
    input:
        json="data_input/spectronaut_schema.json",
        fasta=lambda wc: FASTA,
        htrms=lambda wc: glob.glob("data_input/converted_files/**/*.htrms", recursive=True)
    output:
        report="data_output/all_workflows/all_Spectronaut_Report.tsv"
    threads: 16
    resources:
        mem_mb=262144
    shell:
        r"""
        mkdir -p data_output/all_workflows
        module load --auto spectronaut/20.0.250602.92449
        spectronaut direct \
            -s {input.json} \
            -n all_workflows \
            -o {output.report:r} \
            -fasta data_input/{FASTA} \
            -d data_input/converted_files
        """

# ---------------------------------------------------------------------
# 8. Final marker rule
# ---------------------------------------------------------------------
rule spectronaut_complete:
    input:
        expand("data_output/{wf}/{wf}_Spectronaut_Report.tsv", wf=WORKFLOWS),
        expand("data_output/{wf}_{ct}/{wf}_{ct}_Spectronaut_Report.tsv", wf=WORKFLOWS, ct=COHORTS),
        "data_output/all_workflows/all_Spectronaut_Report.tsv"
    output:
        marker="data_output/S_all_spectronaut_complete.marker"
    shell:
        "touch {output.marker}"
