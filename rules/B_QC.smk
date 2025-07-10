import pandas as pd

directory_root = config["project_root"]
data_dir = f"{directory_root}data_output/"
# Define workflows and cohorts before using them
metadata = pd.read_csv(f"{directory_root}data_input/ms_raw_files.csv", delimiter=";")
metadata.columns = metadata.columns.str.strip()

# List of unique workflows and cohorts (if needed)
workflows = metadata["workflow"].unique()
cohorts = sorted(metadata["cohort"].unique())


if config["batch_correction"] == False and config["exclude_samples"] == []:
    qc_version = "QC"
elif config["batch_correction"] == True and config["exclude_samples"] == []:
    qc_version = "QC_batch_correction"
elif config["batch_correction"] == True and config["exclude_samples"] != []:
    qc_version = "QC_remove_samples_and_batch_effect"
    data_dir = f"{directory_root}data_output/removed_samples_data/"
else:
    qc_version = "QC_remove_samples"
    data_dir = f"{directory_root}data_output/removed_samples_data/"

ruleorder: run_alphastats_analysis_cohort > run_alphastats_analysis_workflow

rule setup_env:
    output:
        env_flag=f"{directory_root} .env_setup"
    shell:
        """
        curl -LsSf https://astral.sh/uv/install.sh | sh
        module purge
        uv venv --python 3.9
        uv pip install venn
        uv pip install alphastats
        uv pip install adjustText
        uv pip install snakemake
        touch {output.env_flag}
        """
rule run_alphastats_analysis_workflow:
    input:
        matrix= f"{data_dir}{{workflow}}/{{workflow}}.pg_matrix.tsv"
    output:
        new_matrix=f"{data_dir}{{workflow}}/{{workflow}}_corrected.pg_matrix.tsv",
        output_dir = directory(f"{directory_root}data_output/{qc_version}/{{workflow}}"),
        report = f"{directory_root}data_output/{qc_version}/{{workflow}}/{{workflow}}_report.pdf"
    threads: 16
    resources:
        mem_mb = 262144,
        slurm_partition = "standardqueue"
    shell:
        """
        uv run {directory_root}scripts/B_QC.py --matrix {input.matrix} \
      --metadata  {config[metadata]} \
      --contamination {config[contamination_panel_file]} \
      --outputMatrix {output.new_matrix} \
      --outputReport {output.report} \
      --output {output.output_dir} \
      --filtering_option {config[filtering_option]} \
      --group {config[group]} \
      --filtering_percentage {config[filtering_percentage]} \
      --pca_factors '{config[pca_factors]}' \
      --missing_values_group '{config[missing_values_grouping]}' \
      --contamination_panel_flag {config[contamination_panel]} \
      --exclude_samples '{config[exclude_samples]}' \
      --batch_column {config[batch_column]}  \
      --batch_correction {config[batch_correction]}  \
      --batch_correction_column {config[batch_effect_column]}
        """

rule run_alphastats_analysis_cohort:
    input:
        matrix=f"{data_dir}{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}.pg_matrix.tsv"
    output:
        new_matrix=f"{data_dir}{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}_corrected.pg_matrix.tsv",
        output_dir = directory(f"{directory_root}data_output/{qc_version}/{{workflow}}_{{cohort}}"),
        report = f"{directory_root}data_output/{qc_version}/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}_report.pdf"
    threads: 16
    resources:
        mem_mb=262144,
        slurm_partition="standardqueue"
    shell:
        """
        uv run {directory_root}scripts/B_QC.py --matrix {input.matrix} \
              --metadata  {config[metadata]} \
      --contamination {config[contamination_panel_file]} \
      --outputMatrix {output.new_matrix} \
      --outputReport {output.report} \
      --output {output.output_dir} \
      --filtering_option {config[filtering_option]} \
      --group {config[group]} \
      --filtering_percentage {config[filtering_percentage]} \
      --pca_factors '{config[pca_factors]}' \
      --missing_values_group '{config[missing_values_grouping]}' \
      --contamination_panel_flag {config[contamination_panel]} \
      --exclude_samples '{config[exclude_samples]}' \
      --batch_column {config[batch_column]}  \
      --batch_correction {config[batch_correction]} \
      --batch_correction_column {config[batch_effect_column]}
        """

rule run_all_QC:
    input:
        expand(f"{data_dir}{{workflow}}/{{workflow}}_corrected.pg_matrix.tsv",
               workflow=workflows),
        expand(f"{data_dir}{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}_corrected.pg_matrix.tsv",
               workflow=workflows, cohort=cohorts)
    output:
        marker = f"{directory_root}data_output/{qc_version}/B_QC_individual_complete.marker"
    shell:
        "touch {output.marker}"

rule summarize_qc:
    input:
        marker = f"{directory_root}data_output/{qc_version}/B_QC_individual_complete.marker",
        dir_input = f"{directory_root}data_output/{qc_version}"
    output:
        directory = directory(f"{directory_root}data_output/{qc_version}/summary_plots"),
        marker = f"{directory_root}data_output/{qc_version}/B_QC_summary_complete.marker"
    resources:
        mem_mb=262144,
        slurm_partition="standardqueue"
    shell:
        """
        uv run {directory_root}scripts/summary_report.py --input {input.dir_input} --output {output.directory}
        touch {output.marker}
        """




