import pandas as pd

directory_root = config["project_root"]

# --- Data Preparation ---

# Load file mapping raw files to workflows and cohorts
metadata = pd.read_csv(f"{directory_root}data_input/ms_raw_files.csv", delimiter=";")
metadata.columns = metadata.columns.str.strip()

# Remove '.raw' from raw file names
metadata['file_name'] = metadata['file_name'].str.replace('.raw', '', regex=False)

# List of unique workflows and cohorts (if needed)
workflows = metadata["workflow"].unique()
cohorts = sorted(metadata["cohort"].unique())

# --- Rules ---

rule remove_excluded_samples:
    input:
        indir = f"{directory_root}data_input/converted_files"
    output:
        marker = f"{directory_root}data_input/converted_files/.removed_excluded"
    params:
        exclude = ",".join(config["exclude_samples"])
    shell:
        """
        python {directory_root}scripts/remove_samples.py --indir {input.indir} --exclude {params.exclude}
        touch {output.marker}
        """



# Rule: DIANN Analysis for Entire Workflows
# This runs DIANN on all converted files for a workflow, regardless of cohort (using --dir-all).
rule diann_analysis_workflows_removed:
    input:
        library = f"{directory_root}data_output/library.predicted.speclib",
        raw_data_dir = f"{directory_root}data_input/converted_files/{{workflow}}"
    output:
        results = f"{directory_root}data_output/removed_samples_data/{{workflow}}/{{workflow}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144,
        slurm_partition = "standardqueue"
    shell:
        """
        mkdir -p {directory_root}data_output/removed_samples_data/{wildcards.workflow}
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto diann/2.0.1
        module load gcc/13.2.0
        diann --lib {input.library} \
              --dir-all {input.raw_data_dir} \
              --threads {threads} \
              --out {directory_root}data_output/removed_samples_data/{wildcards.workflow}/{wildcards.workflow} \
              --qvalue 0.01 \
              --matrix-qvalue 0.01 \
              --mass-acc-ms1 10 \
              --mass-acc 10 \
              --matrices \
              --missed-cleavages 1
        """

# Rule: DIANN Analysis for Each Cohort Separately
# This runs DIANN on each workflow/cohort folder individually.

rule diann_analysis_cohorts_removed:
    input:
        library = f"{directory_root}data_output/library.predicted.speclib",
        raw_data_dir = f"{directory_root}data_input/converted_files/{{workflow}}/{{cohort}}"
    output:
        results = f"{directory_root}data_output/removed_samples_data/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144,
        slurm_partition = "standardqueue"
    shell:
        """
        mkdir -p {directory_root}data_output/removed_samples_data/{wildcards.workflow}_{wildcards.cohort}
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto diann/2.0.1
        module load gcc/13.2.0
        diann --lib {input.library} \
              --dir {input.raw_data_dir} \
              --threads {threads} \
              --out {directory_root}data_output/removed_samples_data/{wildcards.workflow}_{wildcards.cohort}/{wildcards.workflow}_{wildcards.cohort} \
              --qvalue 0.01 \
              --matrix-qvalue 0.01 \
              --mass-acc-ms1 10 \
              --mass-acc 10 \
              --matrices \
              --missed-cleavages 1
        """

rule A_quantify_all:
    input:
        # Workflow-level DIANN results
        expand(f"{directory_root}data_output/removed_samples_data/{{workflow}}/{{workflow}}.stats.tsv", workflow=workflows),
        # Cohort-level DIANN results: for each row in metadata, construct a target path for that workflow/cohort
        expand(f"{directory_root}data_output/removed_samples_data/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",
               workflow=workflows, cohort=sorted(metadata["cohort"].unique()))
    output:
        marker = f"{directory_root}data_output/A_all_diann_removed_complete.marker"
    shell:
        "touch {output.marker}"


rule remove_all_samples:
    input:
        f"{directory_root}data_input/converted_files/.removed_excluded",
        f"{directory_root}data_input/converted_files",
        expand(f"{directory_root}data_output/removed_samples_data/{{workflow}}/{{workflow}}.stats.tsv",workflow=workflows),
        expand(f"{directory_root}data_output/removed_samples_data/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",workflow=workflows, cohort=sorted(metadata["cohort"].unique())),
        f"{directory_root}data_output/A_all_diann_removed_complete.marker"

