#!/usr/bin/env python3
import urllib.request
import pathlib
import tempfile
import shutil

# --- Settings ---
organism_id = "9606"         # Human; use 10090 for mouse
include_isoforms = "false"    # "true" or "false"
reviewed = True              # Swiss-Prot only if True
exclude_fragments = True     # Drop fragments if True

# --- Build query ---
query = f"organism_id:{organism_id}"
if reviewed:
    query += " AND reviewed:true"
if exclude_fragments:
    query += " AND fragment:false"

params = (
    f"https://rest.uniprot.org/uniprotkb/stream"
    f"?query={urllib.request.quote(query)}"
    f"&format=fasta"
    f"&includeIsoform={include_isoforms}"
)

# --- Output path ---
downloads = pathlib.Path.home() / "Downloads"
output_file = downloads / "uniprot_reviewed_human.fasta"
output_file.parent.mkdir(parents=True, exist_ok=True)

print(f"[INFO] Downloading from UniProt to {output_file} ...")

# --- Download and save ---
with urllib.request.urlopen(params) as response, tempfile.NamedTemporaryFile(delete=False) as tmp:
    shutil.copyfileobj(response, tmp)
    tmp_path = tmp.name

shutil.move(tmp_path, output_file)
print(f"[OK] FASTA saved: {output_file}")


# after running this script, you can examine how many entries were downloaded with: grep -c "^>" ~/Downloads/uniprot_reviewed_human.fasta


