#!/usr/bin/env python3
"""
Single-cell TelomereHunter2: Filtering, Sorting, and Summary for Barcodes

This script performs telomere analysis on single-cell sequencing data:
1. Runs TH2 bulk analysis with single-cell mode enabled
2. Processes barcodes above a read count threshold
3. Calculates per-barcode telomere content, TVR occurrences, and singleton variants
4. Outputs a comprehensive summary table with normalized metrics

The analysis follows the same approach as bulk TH2 but applied to each barcode independently.
"""
import sys
import os
import pandas as pd
import subprocess
import pysam
import numpy as np
from telomerehunter2.telomerehunter2_main import parse_command_line_arguments
from telomerehunter2.TVR_screen import find_prioritized_telomeric_indices
from telomerehunter2.utils import get_reverse_complement

# Default TVR hexamers if not specified in arguments
DEFAULT_TVR_HEXAMERS = [
    "TCAGGG", "TGAGGG", "TTGGGG", "TTCGGG", "TTTGGG",
    "ATAGGG", "CATGGG", "CTAGGG", "GTAGGG", "TAAGGG"
]

# BAM file types for telomere analysis
BAM_TYPES = ["intratelomeric", "junctionspanning", "subtelomeric", "intrachromosomal"]

# Canonical telomere repeat
CANONICAL_REPEAT = "TTAGGG"

# Default quality threshold for base calling
DEFAULT_QUALITY_THRESHOLD = 20


def count_reads_per_barcode_in_bamfiles(bam_files, barcodes):
    """
    Count reads for each barcode in each BAM file (intratelomeric, junctionspanning, subtelomeric, intrachromosomal).

    Args:
        bam_files: Dictionary mapping BAM types to file paths
        barcodes: List of barcodes to count reads for

    Returns:
        Dictionary mapping barcodes to counts per BAM type
    """
    barcode_counts = {bc: {bam_type: 0 for bam_type in bam_files} for bc in barcodes}
    total_barcodes = len(barcodes)

    print(f"Counting reads for {total_barcodes} barcodes across {len(bam_files)} BAM files")
    for bam_type, bam_path in bam_files.items():
        if not os.path.exists(bam_path):
            print(f"BAM file not found: {bam_path} (type: {bam_type})")
            continue
        print(f"Processing {bam_type} BAM file: {bam_path}")
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bamfile:
                read_count = 0
                barcode_matched = 0
                for read in bamfile.fetch(until_eof=True):
                    read_count += 1
                    if read_count % 1000000 == 0:
                        print(f"Processed {read_count} reads, matched {barcode_matched} to barcodes")
                    bc = read.get_tag("CB") if read.has_tag("CB") else None
                    if bc in barcode_counts:
                        barcode_counts[bc][bam_type] += 1
                        barcode_matched += 1
                print(f"Completed {bam_type}: processed {read_count} reads, matched {barcode_matched} to barcodes")
        except Exception as e:
            print(f"Error processing {bam_type} BAM file: {str(e)}")
    return barcode_counts


def screen_bam_all_barcodes(intratel_bam_path, TVR_HEXAMERS, repeat_threshold_used=6, qual_threshold=DEFAULT_QUALITY_THRESHOLD):
    """
    Screen intratelomeric BAM for TVR hexamers and singletons for all barcodes in a single pass.

    Args:
        intratel_bam_path: Path to intratelomeric BAM file
        TVR_HEXAMERS: List of telomere variant repeat (TVR) hexamers to search for
        repeat_threshold_used: Minimum number of repeats required
        qual_threshold: Minimum base quality score required

    Returns:
        Tuple of (barcode_tvr_counts, barcode_singleton_counts, barcode_total_intratelomeric_bp)
    """
    barcode_tvr_counts = {}
    barcode_singleton_counts = {}
    barcode_total_intratelomeric_bp = {}

    if not os.path.exists(intratel_bam_path):
        print(f"Intratelomeric BAM not found: {intratel_bam_path}")
        return barcode_tvr_counts, barcode_singleton_counts, barcode_total_intratelomeric_bp

    print(f"Screening BAM file for TVRs and singletons: {intratel_bam_path}")
    print(f"Using repeat threshold: {repeat_threshold_used}, quality threshold: {qual_threshold}")
    with pysam.AlignmentFile(intratel_bam_path, "rb") as bamfile:
        read_count = 0
        for read in bamfile.fetch(until_eof=True):
            read_count += 1
            if read_count % 1000000 == 0:
                print(f"Processed {read_count} reads...")
            bc = read.get_tag("CB") if read.has_tag("CB") else None
            if not bc:
                continue

            seq = read.query_sequence
            qual = read.query_qualities

            # Initialize data structures for new barcodes
            if bc not in barcode_tvr_counts:
                barcode_tvr_counts[bc] = {tvr: 0 for tvr in TVR_HEXAMERS}
                barcode_singleton_counts[bc] = {tvr: 0 for tvr in TVR_HEXAMERS}
                barcode_total_intratelomeric_bp[bc] = 0

            barcode_total_intratelomeric_bp[bc] += len(seq)

            # Find telomeric pattern indices
            indices = find_prioritized_telomeric_indices(seq, telomere_pattern=r"(.{3})(GGG)")

            # Check reverse complement if not enough repeats found
            if len(indices) < repeat_threshold_used:
                seq = get_reverse_complement(seq)
                qual = read.query_qualities[::-1]  # Reverse quality scores too
                indices = find_prioritized_telomeric_indices(seq, telomere_pattern=r"(.{3})(GGG)")

            # Process each hexamer found
            for i in indices:
                if i + 6 > len(seq):  # Ensure we don't go past end of sequence
                    continue

                pattern_segment = seq[i:i+6]

                # Skip patterns with N or low quality
                if "N" in pattern_segment:
                    continue

                q_scores = qual[i:i+6]
                if any(q < qual_threshold for q in q_scores):
                    continue

                # Count pattern if it's in our list of TVRs
                if pattern_segment in barcode_tvr_counts[bc]:
                    barcode_tvr_counts[bc][pattern_segment] += 1

                    # Check for singleton (TVR flanked by canonical repeats)
                    if i >= 18 and i + 24 <= len(seq):
                        flank_left = seq[i-18:i]
                        flank_right = seq[i+6:i+24]
                        if flank_left == CANONICAL_REPEAT*3 and flank_right == CANONICAL_REPEAT*3:
                            barcode_singleton_counts[bc][pattern_segment] += 1

    print(f"Completed screening. Processed {read_count} reads across {len(barcode_tvr_counts)} barcodes")
    return barcode_tvr_counts, barcode_singleton_counts, barcode_total_intratelomeric_bp


def prepare_summary(barcodes_above, barcode_df, barcode_counts, barcode_tvr_counts, barcode_singleton_counts, barcode_total_intratelomeric_bp, args_th2, summary_df, TVR_HEXAMERS):
    """
    Prepare a summary dataframe for each barcode with telomere content, TVR occurrences, and singleton variants.

    Args:
        barcodes_above: List of barcodes above read count threshold
        barcode_df: DataFrame of barcode information
        barcode_counts: Dictionary mapping barcodes to read counts per BAM type
        barcode_tvr_counts: Dictionary mapping barcodes to TVR counts
        barcode_singleton_counts: Dictionary mapping barcodes to singleton counts
        barcode_total_intratelomeric_bp: Dictionary mapping barcodes to total intratelomeric base pair counts
        args_th2: Arguments from command line
        summary_df: Summary DataFrame from bulk TH2 run
        TVR_HEXAMERS: List of TVR hexamers to analyze

    Returns:
        List of dictionaries, each containing summary data for one barcode
    """
    total_reads_with_tel_gc = int(summary_df['total_reads_with_tel_gc'].values[0])
    total_reads = int(summary_df['total_reads'].values[0])
    pid = str(summary_df['PID'].values[0])
    read_length = str(summary_df['read_lengths'].values[0])
    repeat_threshold_set = str(summary_df['repeat_threshold_set'].values[0])
    repeat_threshold_used = str(summary_df['repeat_threshold_used'].values[0])
    gc_bins_for_correction = str(summary_df['gc_bins_for_correction'].values[0])

    print(f"Preparing summary for {len(barcodes_above)} barcodes")

    summary_rows = []
    for barcode in barcodes_above:
        # Get barcode-specific counts
        bc_total_reads = int(barcode_df[barcode_df['barcode'] == barcode]['read_count'].values[0])
        intratelomeric_reads = barcode_counts[barcode]["intratelomeric"]
        junctionspanning_reads = barcode_counts[barcode]["junctionspanning"]
        subtelomeric_reads = barcode_counts[barcode]["subtelomeric"]
        intrachromosomal_reads = barcode_counts[barcode]["intrachromosomal"]

        # Calculate telomere read count and content
        tel_read_count = intratelomeric_reads + junctionspanning_reads + subtelomeric_reads + intrachromosomal_reads
        bc_total_reads_with_tel_gc = (total_reads_with_tel_gc / total_reads) * bc_total_reads
        tel_content = (float(intratelomeric_reads) / float(bc_total_reads_with_tel_gc) * 1e6) if bc_total_reads_with_tel_gc > 0 else np.nan

        # Get TVR and singleton counts
        patterns = barcode_tvr_counts.get(barcode, {tvr: 0 for tvr in TVR_HEXAMERS})
        singleton_counts = barcode_singleton_counts.get(barcode, {tvr: 0 for tvr in TVR_HEXAMERS})
        summed_intratel_read_length = barcode_total_intratelomeric_bp.get(barcode, 0)

        # Calculate normalization values
        context_norm = {}
        context_norm_100bp = {}
        singleton_norm = {}

        for tvr in TVR_HEXAMERS:
            count = patterns.get(tvr, 0)
            context_norm[f"{tvr}_arbitrary_context_norm_by_intratel_reads"] = (count / intratelomeric_reads) if intratelomeric_reads > 0 else np.nan
            context_norm_100bp[f"{tvr}_arbitrary_context_per_100bp_intratel_read"] = (count * 100 / summed_intratel_read_length) if summed_intratel_read_length and summed_intratel_read_length > 0 else np.nan
            singleton_norm[f"{tvr}_singletons_norm_by_all_reads"] = (singleton_counts.get(tvr, 0) / bc_total_reads) if bc_total_reads > 0 else np.nan

        # Create summary row
        summary_row = {
            'PID': pid,
            'sample': barcode,
            'tel_content': tel_content,
            'total_reads': bc_total_reads,
            'read_lengths': read_length,
            'repeat_threshold_set': repeat_threshold_set,
            'repeat_threshold_used': repeat_threshold_used,
            'intratelomeric_reads': intratelomeric_reads,
            'junctionspanning_reads': junctionspanning_reads,
            'subtelomeric_reads': subtelomeric_reads,
            'intrachromosomal_reads': intrachromosomal_reads,
            'tel_read_count': tel_read_count,
            'gc_bins_for_correction': gc_bins_for_correction,
            'total_reads_with_tel_gc': bc_total_reads_with_tel_gc,
            **context_norm,
            **context_norm_100bp,
            **singleton_norm,
        }
        summary_rows.append(summary_row)

    print(f"Generated {len(summary_rows)} summary rows")
    return summary_rows


def main():
    """
    Main function to run single-cell TelomereHunter2 analysis.
    1. Run bulk TH2 analysis
    2. Process barcodes and count reads
    3. Screen for TVR hexamers and singletons
    4. Generate summary table
    """
    print("Running in TH2 single-cell mode: full file analysis...")

    # Parse command line arguments
    args_th2 = parse_command_line_arguments()

    # Ensure single-cell mode and disable plotting
    cli_args = list(sys.argv[1:])
    if '--plotNone' not in cli_args:
        cli_args.append('--plotNone')
    if '--singlecell_mode' not in cli_args:
        cli_args.append('--singlecell_mode')

    # Run bulk TH2 analysis first
    th2_cmd = [sys.executable, '-m', 'telomerehunter2.telomerehunter2_main', *cli_args]
    print(f"Running bulk TH2: {' '.join(map(str, th2_cmd))}")
    try:
        result = subprocess.run(th2_cmd, check=True)
        if result.returncode != 0:
            print("Bulk TH2 run failed. Exiting.")
            sys.exit(result.returncode)
    except subprocess.CalledProcessError as e:
        print(f"Bulk TH2 run failed with error: {str(e)}")
        sys.exit(1)

    # Load results from bulk run for SC analysis
    barcode_table = os.path.join(args_th2.outdir, f"{args_th2.pid}_barcode_counts.tsv")
    summary_file = os.path.join(args_th2.outdir, f"{args_th2.pid}_summary.tsv")

    if not os.path.exists(barcode_table):
        print(f"Barcode table not found: {barcode_table}")
        return
    if not os.path.exists(summary_file):
        print(f"Summary file not found: {summary_file}")
        return

    # Read result files
    try:
        barcode_df = pd.read_table(barcode_table)
        summary_df = pd.read_table(summary_file)
    except Exception as e:
        print(f"Error reading result files: {str(e)}")
        return

    # Check for required fields
    if 'total_reads_with_tel_gc' not in summary_df.columns:
        print("Required field 'total_reads_with_tel_gc' not found in summary file. Exiting.")
        sys.exit(1)

    # Filter barcodes by read count threshold
    barcodes_above = barcode_df[barcode_df['read_count'] >= args_th2.min_reads_per_barcode]['barcode'].tolist()
    print(f"Barcodes above threshold ({args_th2.min_reads_per_barcode} reads): {len(barcodes_above)}")
    if not barcodes_above:
        print("No barcodes above the threshold. Exiting.")
        return

    # Get BAM files and count reads per barcode
    bam_files_all = {
        bam_type: os.path.join(args_th2.outdir, f"tumor_TelomerCnt_{args_th2.pid}", f"{args_th2.pid}_filtered_{bam_type}.bam")
        for bam_type in BAM_TYPES
    }

    # Count reads per barcode in each BAM file
    barcode_counts = count_reads_per_barcode_in_bamfiles(bam_files=bam_files_all, barcodes=barcodes_above)

    # Use TVR hexamers from args or default list
    TVR_HEXAMERS = args_th2.TVRs_for_context if hasattr(args_th2, 'TVRs_for_context') else DEFAULT_TVR_HEXAMERS

    # Screen intratelomeric BAM for TVRs and singletons
    intratel_bam_path = os.path.join(
        args_th2.outdir,
        f"tumor_TelomerCnt_{args_th2.pid}",
        f"{args_th2.pid}_filtered_intratelomeric.bam"
    )

    barcode_tvr_counts, barcode_singleton_counts, barcode_total_intratelomeric_bp = screen_bam_all_barcodes(
        intratel_bam_path=intratel_bam_path,
        TVR_HEXAMERS=TVR_HEXAMERS,
        repeat_threshold_used=int(summary_df['repeat_threshold_used'].values[0]),
        qual_threshold=DEFAULT_QUALITY_THRESHOLD
    )

    # Prepare summary rows
    summary_rows = prepare_summary(
        barcodes_above,
        barcode_df,
        barcode_counts,
        barcode_tvr_counts,
        barcode_singleton_counts,
        barcode_total_intratelomeric_bp,
        args_th2,
        summary_df,
        TVR_HEXAMERS
    )

    print(f"Summary rows generated: {len(summary_rows)}")
    if summary_rows:
        # Get dynamic column names based on the TVR hexamers
        tvr_arbitrary = [f"{tvr}_arbitrary_context_norm_by_intratel_reads" for tvr in TVR_HEXAMERS]
        tvr_100bp = [f"{tvr}_arbitrary_context_per_100bp_intratel_read" for tvr in TVR_HEXAMERS]
        tvr_singletons = [f"{tvr}_singletons_norm_by_all_reads" for tvr in TVR_HEXAMERS]

        # Standard columns
        std_columns = [
            'PID', 'sample', 'tel_content', 'total_reads', 'read_lengths', 'repeat_threshold_set', 'repeat_threshold_used',
            'intratelomeric_reads', 'junctionspanning_reads', 'subtelomeric_reads', 'intrachromosomal_reads',
            'tel_read_count', 'gc_bins_for_correction', 'total_reads_with_tel_gc'
        ]

        # All columns
        columns = std_columns + tvr_arbitrary + tvr_100bp + tvr_singletons

        # Create DataFrame and select columns
        summary_df = pd.DataFrame(summary_rows)

        # Ensure all expected columns exist
        for col in columns:
            if col not in summary_df.columns:
                print(f"Expected column '{col}' not found in results. Adding empty column.")
                summary_df[col] = np.nan

        summary_final = summary_df[columns]

        # Write to file
        summary_path = os.path.join(args_th2.outdir, f"{args_th2.pid}_sc_summary.tsv")
        summary_final.to_csv(summary_path, sep="\t", index=False)
        print(f"SC summary file written to: {summary_path}")
    else:
        print("No barcode summaries generated.")
    print("Completed barcode-specific analysis and summary writing.")


if __name__ == "__main__":
    main()
