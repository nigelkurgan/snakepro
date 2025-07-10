############################################
# 1️⃣ Load configuration and set workdir
############################################

configfile: "config/config.yml"

# Optionally set the working directory from config:
workdir: config["project_root"]

# Extract key user options
TOOL   = config["quant_tool"].lower()
RUN_QC = config.get("run_qc", True)


############################################
# 2️⃣ Define the top-level default target
############################################

def final_outputs():
    if TOOL == "diann":
        return expand("data_output/{wf}/A_all_diann_complete.marker",
                      wf=config["workflows"])
    elif TOOL == "spectronaut":
        return expand("data_output/{wf}/A_all_spectronaut_complete.marker",
                      wf=config["workflows"])
    else:
        raise ValueError(f"Unsupported quant_tool: {TOOL}")

rule all:
    input:
        final_outputs()


############################################
# 3️⃣ Conditionally include rule modules
############################################

if TOOL == "diann":
    include: "rules/A_diann.smk"
elif TOOL == "spectronaut":
    include: "rules/D_spectronaut.smk"
else:
    raise ValueError(f"Unsupported quant_tool: {TOOL}")

# If QC is requested, include QC rules
if RUN_QC:
    include: "rules/B_QC.smk"

# Only relevant for DIA-NN
if TOOL == "diann":
    include: "rules/C_diann_remove.smk"


