# ---------------------------------------------------------------------
# DIANN sample-removal + re-analysis pipeline
# ---------------------------------------------------------------------

import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------
# 1. Load metadata & global wildcard lists
# ---------------------------------------------------------------------
directory_root = config["project_root"]

metadata = pd.read_csv(f"{directory_root}data_input/ms_raw_files.csv", delimiter=";")
metadata.columns = metadata.columns.str.strip()
WORKFLOWS = sorted(metadata["workflow"].unique())
COHORTS = sorted(metadata["cohort"].unique())

# ---------------------------------------------------------------------
# 2. Remove excluded samples from converted data
# ---------------------------------------------------------------------
rule remove_excluded_samples:
    input:
        indir = f"{directory_root}data_input/converted_files"
    output:
        marker = f"{directory_root}data_input/converted_files/.removed_excluded"
    params:
        exclude = ",".join(config["exclude_samples"])
    shell:
        """
        python {directory_root}scripts/remove_samples.py \
          --indir {input.indir} \
          --exclude {params.exclude}
        touch {output.marker}
        """

# ---------------------------------------------------------------------
# 3. DIANN analysis on workflows (after sample removal)
# ---------------------------------------------------------------------
rule diann_analysis_workflows_removed:
    input:
        library = f"{directory_root}data_output/library.predicted.speclib",
        raw_data_dir = lambda wc: f"{directory_root}data_input/converted_files/{wc.workflow}"
    output:
        stats = f"{directory_root}data_output/removed_samples_data/{{workflow}}/{{workflow}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p {directory_root}data_output/removed_samples_data/{wildcards.workflow}
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
# 4. DIANN analysis on workflows + cohorts (sample removal)
# ---------------------------------------------------------------------
rule diann_analysis_cohorts_removed:
    input:
        library = f"{directory_root}data_output/library.predicted.speclib",
        raw_data_dir = lambda wc: f"{directory_root}data_input/converted_files/{wc.workflow}/{wc.cohort}"
    output:
        stats = f"{directory_root}data_output/removed_samples_data/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144
    shell:
        """
        mkdir -p {directory_root}data_output/removed_samples_data/{wildcards.workflow}_{wildcards.cohort}
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
# 5. Marker rule — all removed-sample DIANN results done
# ---------------------------------------------------------------------
rule A_quantify_all:
    input:
        # workflow level
        expand(
            f"{directory_root}data_output/removed_samples_data/{{workflow}}/{{workflow}}.stats.tsv",
            workflow=WORKFLOWS
        ),
        # cohort level
        expand(
            f"{directory_root}data_output/removed_samples_data/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",
            workflow=WORKFLOWS, cohort=COHORTS
        )
    output:
        marker = f"{directory_root}data_output/A_all_diann_removed_complete.marker"
    shell:
        "touch {output.marker}"

# ---------------------------------------------------------------------
# 6. Cleanup rule — optional, no shell commands
# ---------------------------------------------------------------------
rule remove_all_samples:
    input:
        f"{directory_root}data_input/converted_files/.removed_excluded",
        expand(
            f"{directory_root}data_output/removed_samples_data/{{workflow}}/{{workflow}}.stats.tsv",
            workflow=WORKFLOWS
        ),
        expand(
            f"{directory_root}data_output/removed_samples_data/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",
            workflow=WORKFLOWS, cohort=COHORTS
        ),
        f"{directory_root}data_output/A_all_diann_removed_complete.marker"
    # This rule has no shell—just serves as a DAG node if you want to add cleanup actions.
