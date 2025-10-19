#!/usr/bin/python

# Copyright 2024 Ferdinand Popp, Lina Sieverling, Philip Ginsbach, Lars Feuerbach

# This file is part of TelomereHunter2.

# TelomereHunter2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter2. If not, see <http://www.gnu.org/licenses/>.


import os
import pandas as pd
import pysam
from telomerehunter2 import filter_telomere_reads, get_repeat_threshold

def run_fast_mode(args):
    out_dir = os.path.join(args.parent_outdir, args.pid)
    os.makedirs(out_dir, exist_ok=True)
    band_file = getattr(args, 'banding_file', None)

    tumor_summary = None
    control_summary = None
    if args.tumor_bam:
        tumor_summary = process_fast_mode_sample(args, args.tumor_bam, "tumor", out_dir, band_file)
    if args.control_bam:
        control_summary = process_fast_mode_sample(args, args.control_bam, "control", out_dir, band_file)
    write_fast_mode_combined_summary(args, out_dir, tumor_summary, control_summary)


def process_fast_mode_sample(args, bam_path, sample_name, out_dir, band_file):
    sample_out_dir = os.path.join(out_dir, f"{sample_name}_TelomerCnt_{args.pid}")
    os.makedirs(sample_out_dir, exist_ok=True)
    temp_unmapped_bam = os.path.join(sample_out_dir, f"{args.pid}_{sample_name}_unmapped.bam")
    print(f"Extracting unmapped reads from {bam_path} to {temp_unmapped_bam}")
    pysam.view('-b', '-f', '4', bam_path, '-o', temp_unmapped_bam, catch_stdout=False)
    pysam.index(temp_unmapped_bam)

    with pysam.AlignmentFile(temp_unmapped_bam, "rb") as unmapped_bam:
        unmapped_count = unmapped_bam.count(until_eof=True)
    print(f"Number of unmapped reads: {unmapped_count}")

    with pysam.AlignmentFile(bam_path, "rb" if bam_path.endswith('.bam') else "rc") as bam_file:
        total_count = bam_file.count(until_eof=True)
    print(f"Total number of reads in input: {total_count}")

    (_, read_lengths, _, _, _, _, repeat_threshold) = get_repeat_threshold.get_read_lengths_and_repeat_thresholds(args, None, bam_path)

    filter_telomere_reads.parallel_filter_telomere_reads(
        bam_path=temp_unmapped_bam,
        out_dir=sample_out_dir,
        pid=args.pid,
        sample=sample_name,
        repeat_threshold_calc=repeat_threshold,
        mapq_threshold=args.mapq_threshold,
        repeats=args.repeats,
        consecutive_flag=args.consecutive,
        remove_duplicates=args.remove_duplicates,
        band_file=band_file,
        num_processes=args.cores,
        singlecell_mode=getattr(args, "singlecell_mode", False),
        fast_mode=True,
    )
    print(f"Fast mode filtering complete for {sample_name}.")

    summary = create_fast_mode_sample_summary(args, sample_out_dir, bam_path, sample_name, total_count, read_lengths, repeat_threshold)
    return summary


def create_fast_mode_sample_summary(args, sample_out_dir, bam_path, sample_name, total_count, read_lengths, repeat_threshold):
    readcount_path = os.path.join(sample_out_dir, f"{args.pid}_readcount.tsv")
    filtered_bam_path = os.path.join(sample_out_dir, f"{args.pid}_filtered.bam")
    trpm_threshold = repeat_threshold if isinstance(repeat_threshold, int) else str(repeat_threshold)
    read_lengths_str = str(read_lengths) if read_lengths is not None else "NA"

    unmapped_reads = 0
    if os.path.exists(readcount_path):
        try:
            readcount_df = pd.read_csv(readcount_path, sep='\t')
            unmapped_reads = readcount_df[readcount_df['chr'] == 'unmapped']['reads'].sum()
        except Exception as e:
            print(f"Error reading readcount.tsv: {e}")

    tel_read_count = 0
    if os.path.exists(filtered_bam_path):
        try:
            with pysam.AlignmentFile(filtered_bam_path, "rb") as bam:
                tel_read_count = bam.count(until_eof=True)
        except Exception as e:
            print(f"Error reading filtered.bam for tel_read_count: {e}")

    TRPM = (tel_read_count / total_count * 1e6) if total_count > 0 else 0.0

    summary = {
        "PID": args.pid,
        "sample": sample_name,
        "tel_reads_per_million_reads": f"{TRPM:.6f}",
        "total_reads": total_count,
        "read_lengths": read_lengths_str,
        "repeat_threshold_set": trpm_threshold,
        "repeat_threshold_used": trpm_threshold,
        "tel_read_count": tel_read_count,
        "unmapped_reads": unmapped_reads
    }

    summary_path = os.path.join(sample_out_dir, f"{args.pid}_{sample_name}_summary.tsv")
    header = ["PID", "sample", "tel_reads_per_million_reads", "total_reads", "read_lengths", "repeat_threshold_set", "repeat_threshold_used", "tel_read_count", "unmapped_reads"]
    with open(summary_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(str(summary[h]) for h in header) + "\n")
    print(f"Summary file written for {sample_name}: {summary_path}")
    return summary


def write_fast_mode_combined_summary(args, out_dir, tumor_summary, control_summary):
    header = ["PID", "sample", "tel_content", "total_reads", "read_lengths", "repeat_threshold_set", "repeat_threshold_used", "tel_read_count", "unmapped_reads"]
    rows = []
    tumor_trpm = None
    control_trpm = None

    def map_row(summary):
        return [str(summary.get(col, "NA")) for col in header]

    if tumor_summary:
        rows.append(map_row(tumor_summary))
        tumor_trpm = float(tumor_summary.get("TRPM", 0))
    if control_summary:
        rows.append(map_row(control_summary))
        control_trpm = float(control_summary.get("TRPM", 0))

    combined_summary_path = os.path.join(out_dir, f"{args.pid}_summary.tsv")
    with open(combined_summary_path, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"Combined summary file written: {combined_summary_path}")
