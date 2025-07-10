# Spectronaut pipeline

conda: "envs/proteomics-spectronaut.yml"

# ---------------------------------------------------------------------
# 1. Setup is done in config.yml file and work dir is already set
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# 2. Load and validate metadata and FASTA file
# ---------------------------------------------------------------------
metadata = pd.read_csv(directory_root / "data_input/ms_raw_files.csv", delimiter=";")
metadata.columns = metadata.columns.str.strip()

# Ensure required columns exist
for col in ["file_name", "workflow", "cohort"]:
    if col not in metadata.columns:
        sys.exit(f"[metadata error] Missing column: '{col}'")

# Split extensions (e.g., '.raw', '.d', '.wiff') and base name
metadata["ext"] = metadata["file_name"].apply(lambda x: os.path.splitext(x)[1].lower())
metadata["file_name"] = metadata["file_name"].apply(lambda x: os.path.splitext(x)[0])

# Print unique extensions for reference
print("Unique file extensions found:", metadata["ext"].unique())

# Collect unique sample, workflow, cohort values
workflows = sorted(metadata["workflow"].unique())
cohorts   = sorted(metadata["cohort"].unique())
samples   = sorted(metadata["file_name"].unique())
exts       = sorted(metadata["ext"].unique())

# fetch fasta fil
# Accept common FASTA suffixes
fasta_candidates = glob.glob(str(directory_root / "data_input" / "*.fa*"))

# Ensure exactly ONE FASTA file is present
if len(fasta_candidates) == 0:
    sys.exit("[FASTA error] No FASTA file found in data_input/")
elif len(fasta_candidates) > 1:
    sys.exit(f"[FASTA error] Expected exactly one FASTA file in data_input/, "
             f"found: {', '.join(Path(p).name for p in fasta_candidates)}")

# Save as Path object for convenience
FASTA = Path(fasta_candidates[0]).resolve()

# Print to stdout so the user sees which FASTA is being used
print(f"Using reference FASTA: {FASTA.name}")

# ---------------------------------------------------------------------
# 3. Convert raw files (.raw, .d, .wiff, etc.) → .htrms
# ---------------------------------------------------------------------
rule convert_to_htrms:
    """
    Convert vendor format (.raw/.d/.wiff) to .htrms using Spectronaut CLI.
    """
    input:
        raw="data_input/raw_files/{sample}.{ext}"
    output:
         htrms = "data_input/converted_files/{{workflow}}/{{cohort}}/{{sample}}.htrms"
    params:
        workflow = lambda wc: metadata[metadata["file_name"] == wc.sample]["workflow"].values[0],
        cohort = lambda wc: metadata[metadata["file_name"] == wc.sample]["cohort"].values[0],
        outdir = lambda wc: f"{directory_root}data_input/converted_files/{metadata[metadata['file_name'] == wc.sample]['workflow'].values[0]}/{metadata[metadata['file_name'] == wc.sample]['cohort'].values[0]}"
    threads: 4
    resources:
        mem_mb=64000,
        slurm_partition="standardqueue"
    message:
        "Converting {input.raw} → {output.htrms}"
    shell:
        r"""
        mkdir -p {params.outdir}
        source /projects/cbmr_shared/apps/modules/activate.sh
        module load --auto spectronaut/20.0.250602.92449
        spectronaut -convert -i {input.raw} -o {output.htrms} -nogui
        """

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



# ---------------------------------------------------------------------
# 4. Spectronaut Pipeline – Per workflow (all cohorts pooled)
# ---------------------------------------------------------------------
rule spectronaut_workflow:
    """
    Run Spectronaut DirectDIA on all samples from a workflow (all cohorts pooled).
    """
    input:
        json  = directory_root / "data_input/spectronaut_schema.json",
        fasta = directory_root / "data_input/uniprotkb_reviewed_true_AND_model_organ_2025_02_10.fasta",
        hdir  = lambda wc: directory_root / "data_input/converted_files" / wc.workflow
    output:
        report = directory_root / "data_output/{workflow}/{workflow}_Spectronaut_Report.tsv"
    threads: 16
    resources:
        mem_mb=262144,
        slurm_partition="standardqueue"
    shell:
        r"""
        mkdir -p {directory_root}/data_output/{wildcards.workflow}
        module load --auto spectronaut/20.0.250602.92449
        spectronaut direct \
            -s {input.json} \
            -n {wildcards.workflow} \
            -o {directory_root}/data_output/{wildcards.workflow} \
            -fasta {input.fasta} \
            -d {input.hdir}
        """


# ---------------------------------------------------------------------
# 5. Spectronaut Pipeline – Per workflow + cohort
# ---------------------------------------------------------------------
rule spectronaut_cohort:
    """
    Run Spectronaut DirectDIA per workflow+cohort subset.
    """
    input:
        json  = directory_root / "data_input/spectronaut_schema.json",
        fasta = FASTA,
        hdir  = lambda wc: directory_root / "data_input/converted_files" / wc.workflow / wc.cohort
    output:
        report = directory_root / "data_output/{workflow}_{cohort}/{workflow}_{cohort}_Spectronaut_Report.tsv"
    threads: 16
    resources:
        mem_mb=262144,
        slurm_partition="standardqueue"
    shell:
        r"""
        mkdir -p {directory_root}/data_output/{wildcards.workflow}_{wildcards.cohort}
        module load --auto spectronaut/20.0.250602.92449
        spectronaut direct \
            -s {input.json} \
            -n {wildcards.workflow}_{wildcards.cohort} \
            -o {directory_root}/data_output/{wildcards.workflow}_{wildcards.cohort} \
            -fasta {input.fasta} \
            -d {input.hdir}
        """


# ---------------------------------------------------------------------
# 6. pectronaut Pipeline – All files combined
# ---------------------------------------------------------------------
rule spectronaut_all:
    """
    Run Spectronaut DirectDIA across all samples from all workflows/cohorts.
    """
    input:
        json  = directory_root / "data_input/spectronaut_schema.json",
        fasta = FASTA,
        hdir  = directory_root / "data_input/converted_files"
    output:
        report = directory_root / "data_output/all_workflows/all_Spectronaut_Report.tsv"
    threads: 16
    resources:
        mem_mb=262144,
        slurm_partition="standardqueue"
    shell:
        r"""
        mkdir -p {directory_root}/data_output/all_workflows
        module load --auto spectronaut/20.0.250602.92449
        spectronaut direct \
            -s {input.json} \
            -n all_workflows \
            -o {directory_root}/data_output/all_workflows \
            -fasta {input.fasta} \
            -d {input.hdir}
        """


# ---------------------------------------------------------------------
# 7. Ensure all report exist
# ---------------------------------------------------------------------
rule S_all:
    """
    Final marker rule – checks that all analysis reports exist.
    """
    input:
        # Workflow-level
        expand(
            f"{directory_root}/data_output/{{workflow}}/{{workflow}}_Spectronaut_Report.tsv",
            workflow=workflows
        ),
        # Cohort-level
        expand(
            f"{directory_root}/data_output/{{workflow}}_{{cohort}}/{{workflow}}_{{cohort}}_Spectronaut_Report.tsv",
            workflow=workflows, cohort=cohorts
        ),
        # All combined
        f"{directory_root}/data_output/all_workflows/all_Spectronaut_Report.tsv"
    output:
        marker = f"{directory_root}/data_output/S_all_spectronaut_complete.marker"
    shell:
        "touch {output.marker}"
