# rules/00_fetch_fasta.smk
# Defines FASTA and (if needed) a rule to fetch the UniProt FASTA.
# Downstream rules can just use: input: fasta=FASTA

from pathlib import Path

# --------------------------
# Config
# --------------------------
FASTA_CFG   = config.get("fasta", "fetch_fasta")   # path or "fetch_fasta"
OUT_DIR     = config.get("fasta_out_dir", "data_input")
ORG_ID      = str(config.get("organism_id", 9606))
INCLUDE_ISO = bool(config.get("include_isoforms", False))
REVIEWED    = bool(config.get("reviewed", True))           # Swiss-Prot only if True
EXCL_FRAG   = bool(config.get("exclude_fragments", True))  # drop fragments if True
PREFIX      = config.get("filename_prefix", "uniprot")     # name only

# --------------------------
# Helpers
# --------------------------
def fetched_fasta_name():
    tags = [f"org{ORG_ID}"]
    if REVIEWED:    tags.append("reviewed")
    if INCLUDE_ISO: tags.append("isoforms")
    if EXCL_FRAG:   tags.append("noFrag")
    return str(Path(OUT_DIR) / (PREFIX + "_" + "_".join(tags) + ".fasta"))

# --------------------------
# Define FASTA + optional rule
# --------------------------
if FASTA_CFG == "fetch_fasta":
    FASTA = fetched_fasta_name()

    # Make sure the dir exists early (no-op if present)
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    rule fetch_fasta_uniprot:
        output:
            FASTA
        params:
            out_dir     = OUT_DIR,
            org_id      = ORG_ID,
            include_iso = "true" if INCLUDE_ISO else "false",
            reviewed    = "1" if REVIEWED else "0",
            excl_frag   = "1" if EXCL_FRAG else "0",
        log:
            f"{OUT_DIR}/download_uniprot_fasta.log"
        message:
            "Fetching UniProt FASTA (organism_id={params.org_id}, reviewed={params.reviewed}, isoforms={params.include_iso}, excl_frag={params.excl_frag})"
        shell:
            r"""
            set -euo pipefail
            mkdir -p "{params.out_dir}"

            QUERY="organism_id:{params.org_id}"
            if [ "{params.reviewed}" = "1" ]; then
              QUERY="$QUERY AND reviewed:true"
            fi
            if [ "{params.excl_frag}" = "1" ]; then
              QUERY="$QUERY AND fragment:false"
            fi

            tmp=$(mktemp)
            curl -sS -f --retry 3 --retry-delay 5 --max-time 1200 \
                 -L -G "https://rest.uniprot.org/uniprotkb/stream" \
                 --data-urlencode "query=${QUERY}" \
                 --data-urlencode "format=fasta" \
                 --data-urlencode "includeIsoform={params.include_iso}" \
                 --compressed \
            > "$tmp"

            mv "$tmp" "{output}" 2>> "{log}"
            """

else:
    # User provided a path; make it available to all rules as FASTA.
    FASTA = FASTA_CFG

    # Optional: fail early if the file is missing
    rule check_user_fasta_exists:
        input:
            FASTA
        output:
            touch(f"{OUT_DIR}/.user_fasta_ok")
        message:
            "Checking user-provided FASTA exists: {input}"
        shell:
            r"""
            test -s "{input}"
            touch "{output}"
            """

# --------------------------
# How to use FASTA downstream
# --------------------------
# Example:
# rule spectronaut_run:
#     input:
#         fasta=FASTA,
#         raw="data_input/raw_files/"
#     output:
#         "results/{sample}/report.tsv"
#     shell:
#         "spectronaut direct -fasta {input.fasta} -d {input.raw} -o {output}"
