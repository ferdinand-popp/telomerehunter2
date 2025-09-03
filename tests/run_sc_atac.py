#!/usr/bin/env python3

import os
import subprocess

# === Path to your TH2 wrapper script ===
# Get the directory path of the current script or file
script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
# Combine the script directory and data folder using os.path.join
data_folder = os.path.join(script_dir, "data")

script_path = os.path.join(
    os.getcwd(), "..", "src", "telomerehunter2", "sc_barcode_splitter_run.py"
)

# === SET YOUR PARAMETERS HERE ===
input_bam = data_folder + "/atac_hgmm_1k_nextgem_possorted_bam_subsampled_10pct.bam"
cytoband_file = os.path.join(
    script_dir, "src", "telomerehunter2", "cytoband_files", "hg38_cytoBand.txt"
)
patient_id = "Patient5"
output_dir = script_dir + "/results" + f"/scATAC_{patient_id}"
max_parallel = os.cpu_count()

# === Construct command ===
cmd = [
    "python",
    str(script_path),
    "-ibt",
    input_bam,
    "-b",
    cytoband_file,
    "-p",
    patient_id,
    "-o",
    output_dir,
    "--max-parallel",
    str(max_parallel),
    "--max-barcodes",
    "50",
    "--keep-bams",
    "--min-reads-per-barcode",
    "10000",
    "--steps",
    "th2",
]

print("Running TelomereHunter2 SC wrapper with command:")
print(*cmd)

# Run the command and capture its output
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)

# Iterate over the output line by line and print it in real-time
for line in proc.stdout:
    print(line, end="")  # The end='' argument is used to prevent adding extra newlines

# Wait for the command to complete
proc.wait()
