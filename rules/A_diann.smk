import pandas as pd

directory_root = config["project_root"]
print()

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

# Rule: Generate Predicted Spectral Library
rule generate_spectral_library:
    input:
        fasta = f"{directory_root}data_input/uniprotkb_reviewed_true_AND_model_organ_2025_02_10.fasta"
    output:
        library = f"{directory_root}data_output/library.predicted.speclib"
    threads: 16
    shell:
        """
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto diann/2.0.1
        module load gcc/13.2.0
        diann --fasta {input.fasta} \
              --predictor \
              --fasta-search \
              --gen-spec-lib \
              --out-lib {output.library} \
              --threads {threads}
        """

# Rule: Convert Raw Files to mzML
rule convert_raw_files:
    input:
        raw = f"{directory_root}data_input/raw_files/{{sample}}.raw"
    output:
        mzml = f"{directory_root}data_input/converted_files/{{workflow}}/{{cohort}}/{{sample}}.mzML"
    params:
        workflow = lambda wc: metadata[metadata["file_name"] == wc.sample]["workflow"].values[0],
        cohort = lambda wc: metadata[metadata["file_name"] == wc.sample]["cohort"].values[0],
        outdir = lambda wc: f"{directory_root}data_input/converted_files/{metadata[metadata['file_name'] == wc.sample]['workflow'].values[0]}/{metadata[metadata['file_name'] == wc.sample]['cohort'].values[0]}"
    threads: 16
    resources:
        mem_mb = 262144,
        slurm_partition = "standardqueue"
    shell:
        """
        mkdir -p {params.outdir}
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto msconvert/20250218_1
        msconvert --64 --zlib --filter "peakPicking" --filter "zeroSamples removeExtra 1-" --outdir {params.outdir} {input.raw}
        """

rule convert_all:
    input:
        [f"{directory_root}data_input/converted_files/{row['workflow']}/{row['cohort']}/{row['file_name']}.mzML"
         for idx, row in metadata.iterrows()]
    output:
        marker = f"{directory_root}data_output/A_all_converted.marker"
    shell:
        "touch {output.marker}"


# Rule: DIANN Analysis for Entire Workflows
# This runs DIANN on all converted files for a workflow, regardless of cohort (using --dir-all).
rule diann_analysis_workflows:
    input:
        library = f"{directory_root}data_output/library.predicted.predicted.speclib",
        raw_data_dir = f"{directory_root}data_input/converted_files/{{workflow}}"
    output:
        results = f"{directory_root}data_output/{{workflow}}/{{workflow}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144,
        slurm_partition = "standardqueue"
    shell:
        """
        mkdir -p {directory_root}data_output/{wildcards.workflow}
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto diann/2.0.1
        module load gcc/13.2.0
        diann --lib {input.library} \
              --dir-all {input.raw_data_dir} \
              --threads {threads} \
              --out {directory_root}data_output/{wildcards.workflow}/{wildcards.workflow} \
              --qvalue 0.01 \
              --matrix-qvalue 0.01 \
              --mass-acc-ms1 10 \
              --mass-acc 10 \
              --matrices \
              --missed-cleavages 1
        """

# Rule: DIANN Analysis for Each Cohort Separately
# This runs DIANN on each workflow/cohort folder individually.

rule diann_analysis_cohorts:
    input:
        library = f"{directory_root}data_output/library.predicted.predicted.speclib",
        raw_data_dir = f"{directory_root}data_input/converted_files/{{workflow}}/{{cohort}}"
    output:
        results = f"{directory_root}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv"
    threads: 16
    resources:
        mem_mb = 262144,
        slurm_partition = "standardqueue"
    shell:
        """
        mkdir -p {directory_root}data_output/{wildcards.workflow}_{wildcards.cohort}
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto diann/2.0.1
        module load gcc/13.2.0
        diann --lib {input.library} \
              --dir {input.raw_data_dir} \
              --threads {threads} \
              --out {directory_root}data_output/{wildcards.workflow}_{wildcards.cohort}/{wildcards.workflow}_{wildcards.cohort} \
              --qvalue 0.01 \
              --matrix-qvalue 0.01 \
              --mass-acc-ms1 10 \
              --mass-acc 10 \
              --matrices \
              --missed-cleavages 1
        """

rule A_all:
    input:
        # Workflow-level DIANN results
        expand(f"{directory_root}data_output/{{workflow}}/{{workflow}}.stats.tsv", workflow=workflows),
        # Cohort-level DIANN results: for each row in metadata, construct a target path for that workflow/cohort
        expand(f"{directory_root}data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.stats.tsv",
               workflow=workflows, cohort=sorted(metadata["cohort"].unique()))
    output:
        marker = f"{directory_root}data_output/A_all_diann_complete.marker"
    shell:
        "touch {output.marker}"
