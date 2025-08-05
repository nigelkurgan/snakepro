# ---------------------------------------------------------------------
# Load configuration and set workdir
# ---------------------------------------------------------------------


configfile: "config/config.yml"  # Loads parameters from config file :contentReference[oaicite:1]{index=1}

workdir: config["project_root"]  # Set working directory to project root

# Pull key options from config
TOOL   = config["quant_tool"].lower()
RUN_QC = config.get("run_qc", True)

# Read workflows and cohorts from config, or else infer from metadata
WORKFLOWS = config.get("workflows", None)
COHORTS   = config.get("cohorts", None)

if TOOL == "diann" or TOOL == "spectronaut":
    import pandas as pd
    md = pd.read_csv("data_input/ms_raw_files.csv", sep=";")
    md.columns = md.columns.str.strip()
    if WORKFLOWS is None:
        WORKFLOWS = sorted(md["workflow"].unique())
    if COHORTS is None:
        COHORTS = sorted(md["cohort"].unique())

# ---------------------------------------------------------------------
# Define top-level default targets
# ---------------------------------------------------------------------

def final_outputs():
    if TOOL == "diann":
        return expand("data_output/{wf}/A_all_diann_complete.marker", wf=WORKFLOWS)
    elif TOOL == "spectronaut":
        return expand("data_output/S_all_spectronaut_complete.marker", wf=WORKFLOWS)
    else:
        raise ValueError(f"Unsupported quant_tool: {TOOL}")

rule all:
    input:
        final_outputs()

# ---------------------------------------------------------------------
# Include rule modules by tool and QC
# ---------------------------------------------------------------------

if TOOL == "diann":
    include: "rules/01A_diann.smk"
elif TOOL == "spectronaut":
    include: "rules/01B_spectronaut.smk"
else:
    raise ValueError(f"Unsupported quant_tool: {TOOL}")

if RUN_QC:
    include: "rules/02_QC.smk"

# Only run sample-removal branch if DIANN and exclude_samples is set
if TOOL == "diann" and config.get("exclude_samples"):
    include: "rules/03_diann_remove.smk"
