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

import multiprocessing as mp
import os
import re
import shutil
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

import pysam

from telomerehunter2.utils import (get_band_info, get_reverse_complement,
                                   measure_time)


def compile_patterns(repeats):
    patterns_regex_forward = "|".join(repeats)
    patterns_regex_reverse = "|".join(
        [get_reverse_complement(repeat) for repeat in repeats]
    )
    return re.compile(patterns_regex_forward), re.compile(patterns_regex_reverse)


def is_telomere_read(
        consecutive_flag,
        patterns_regex_forward,
        patterns_regex_reverse,
        sequence,
        repeat_threshold_calc,
):
    # Important filtering logic: check if read has the specified amount of patterns, else skip
    if consecutive_flag:
        # Check if consecutive repeats of forward or reverse patterns in the sequence meet the repeat threshold
        return patterns_regex_forward.search(
            f"({patterns_regex_forward.pattern}){{{repeat_threshold_calc},}}", sequence
        ) or patterns_regex_reverse.search(
            f"({patterns_regex_reverse.pattern}){{{repeat_threshold_calc},}}", sequence
        )
    else:
        # Check if the count of forward or reverse patterns in the sequence meets the repeat threshold
        return (
                len(patterns_regex_forward.findall(sequence)) >= repeat_threshold_calc
                or len(patterns_regex_reverse.findall(sequence)) >= repeat_threshold_calc
        )


def initialize_chromosome_and_band_data(bamfile, band_file):
    """Initialize chromosome and band data structures for filtering."""
    references = bamfile.references
    bam_chr_prefix = "chr" if references and references[0].startswith("chr") else ""

    # Use the utility method to get band info
    band_info, sorted_chromosomes = get_band_info(band_file)

    # Create a pickle-compatible structure for band lookups
    bands_dict = {}
    # Track band order for each chromosome
    band_order = {}

    # Process band info into our structure
    for _, row in band_info.iterrows():
        chrom = row["chr"]
        if chrom not in bands_dict:
            bands_dict[chrom] = {"bands": []}
            band_order[chrom] = []

        bands_dict[chrom]["bands"].append({"name": row["band_name"], "end": row["end"]})
        band_order[chrom].append(row["band_name"])

    # Add unmapped category
    bands_dict["unmapped"] = {"bands": [{"name": "unmapped", "end": 0}]}
    band_order["unmapped"] = ["unmapped"]

    # Ensure bands are sorted by end position while maintaining name reference
    for chrom in bands_dict:
        bands_dict[chrom]["bands"].sort(key=lambda x: x["end"])

    return {
        "bam_chr_prefix": bam_chr_prefix,
        "bands": bands_dict,
        "band_order": band_order,
        "chromosome_list": sorted_chromosomes,
        "references": references,
    }


def write_output(out_dir, pid, sample, gc_content_list, read_counts, band_info, barcode_counts=None):
    # Write read counts
    readcount_file_path = os.path.join(out_dir, f"{pid}_readcount.tsv")
    with open(readcount_file_path, "w") as readcount_file:
        readcount_file.write("chr\tband\treads\n")
        # Use the original chromosome order
        for ref in band_info["chromosome_list"]:
            if ref in read_counts:
                # Use the original band order for each chromosome
                for band in band_info["band_order"][ref]:
                    count = read_counts[ref].get(band, 0)
                    readcount_file.write(f"{ref}\t{band}\t{count}\n")

        # Write unmapped reads last if they exist
        if "unmapped" in read_counts:
            readcount_file.write(
                f"unmapped\tunmapped\t{read_counts['unmapped']['unmapped']}\n"
            )

    # Write GC content
    gc_content_file_path = os.path.join(out_dir, f"{pid}_{sample}_gc_content.tsv")
    with open(gc_content_file_path, "w") as gc_content_file:
        gc_content_file.write("gc_content_percent\tread_count\n")
        for gc_content in range(101):
            gc_content_file.write(
                f"{gc_content}\t{gc_content_list.get(gc_content, 0)}\n"
            )

    # Write barcode counts if present
    if barcode_counts is not None:
        barcode_file_path = os.path.join(os.path.dirname(out_dir), f"{pid}_barcode_counts.tsv")
        print(f"Writing barcode table to: {barcode_file_path}")
        complete_entries = sum(1 for bc in barcode_counts if bc)
        print(f"Number of complete barcode entries: {complete_entries}")
        with open(barcode_file_path, "w") as barcode_file:
            barcode_file.write("barcode\tread_count\n")
            for bc, count in barcode_counts.items():
                barcode_file.write(f"{bc}\t{count}\n")


def process_region(args):
    """
    Process a specific region of the BAM file for telomere reads.
    Returns the last file position processed along with other results.
    """
    (
        bam_path,
        region_info,
        patterns_regex_forward,
        patterns_regex_reverse,
        consecutive_flag,
        repeat_threshold_calc,
        mapq_threshold,
        remove_duplicates,
        band_info,
        temp_dir,
        singlecell_mode,
    ) = args

    chrom, start, end = region_info  # unpack tuple
    region_str = f"{chrom}__{start}__{end}"
    temp_bam = os.path.join(temp_dir, f"region_{region_str}_filtered.bam")

    gc_content = {}
    read_counts = {"unmapped": {"unmapped": 0}}
    last_position = 0
    filtered_read_count = 0

    open_mode = "rb" if bam_path.endswith(".bam") else "rc"
    with pysam.AlignmentFile(bam_path, open_mode) as bamfile:
        with pysam.AlignmentFile(temp_bam, "wb", template=bamfile) as filtered_file:
            try:
                barcode_counts = defaultdict(int) if singlecell_mode else None
                for read in bamfile.fetch(contig=chrom, start=start - 1, end=end):
                    try:
                        # Track the last file position
                        current_pos = bamfile.tell()
                        last_position = max(last_position, current_pos)

                        is_unmapped = read.is_unmapped
                        mapping_quality = read.mapping_quality
                        if (
                                read.is_secondary
                                or read.is_supplementary
                                or (remove_duplicates and read.is_duplicate)
                        ):
                            continue

                        # Read check
                        sequence = read.query_sequence
                        try:
                            read_length = len(sequence)
                        except (
                                TypeError
                        ):  # skip if there is no sequence for read in BAM file
                            continue

                        n_count = sequence.count("N")
                        read_length_no_N = read_length - n_count
                        if read_length_no_N > 0 and (n_count / read_length) <= 0.2:
                            gc_percent = int(
                                round(
                                    (sequence.count("C") + sequence.count("G"))
                                    / read_length_no_N
                                    * 100
                                )
                            )
                            gc_content[gc_percent] = gc_content.get(gc_percent, 0) + 1

                        # Process band information
                        if is_unmapped or mapping_quality < mapq_threshold:
                            read_counts["unmapped"]["unmapped"] += 1
                        else:
                            ref_name = read.reference_name
                            # Remove 'chr' prefix if present for matching
                            if ref_name.startswith("chr"):
                                ref_name = ref_name[3:]

                            if ref_name not in band_info["bands"]:
                                read_counts["unmapped"]["unmapped"] += 1
                            else:
                                pos = read.reference_start
                                bands = band_info["bands"][ref_name]["bands"]

                                # Binary search for the correct band
                                left, right = 0, len(bands) - 1
                                found_band = None

                                while left <= right:
                                    mid = (left + right) // 2
                                    if pos <= bands[mid]["end"]:
                                        if mid == 0 or pos > bands[mid - 1]["end"]:
                                            found_band = bands[mid]
                                            break
                                        right = mid - 1
                                    else:
                                        left = mid + 1

                                if found_band is None:
                                    found_band = bands[-1]

                                if ref_name not in read_counts:
                                    read_counts[ref_name] = {}
                                if found_band["name"] not in read_counts[ref_name]:
                                    read_counts[ref_name][found_band["name"]] = 0
                                read_counts[ref_name][found_band["name"]] += 1

                        # Barcode counting
                        if singlecell_mode:
                            bc = read.get_tag("CB") if read.has_tag("CB") else None
                            if bc:
                                barcode_counts[bc] += 1

                        # Check if it's a telomere read
                        if is_telomere_read(
                                consecutive_flag,
                                patterns_regex_forward,
                                patterns_regex_reverse,
                                sequence,
                                repeat_threshold_calc,
                        ):
                            filtered_file.write(read)
                            filtered_read_count += 1
                    except Exception as e:
                        print(f"Error processing read in region {region_str}: {e}")
                        traceback.print_exc()
            except (ValueError, KeyError, MemoryError, OSError) as e:
                print(f"Error processing region {region_str}: {e}")
                traceback.print_exc()
            except Exception as e:
                print(f"Unexpected error in region {region_str}: {e}")
                traceback.print_exc()

    return {
        "region": region_info,
        "gc_content": gc_content,
        "read_counts": read_counts,
        "temp_bam": temp_bam,
        "last_position": last_position,
        "filtered_read_count": filtered_read_count,
        "barcode_counts": dict(barcode_counts) if barcode_counts is not None else {},
    }


def process_unmapped_reads(args):
    """Process unmapped reads from a BAM file starting from a specific position."""
    print("Processing unmapped reads")
    (
        bam_path,
        patterns_regex_forward,
        patterns_regex_reverse,
        consecutive_flag,
        repeat_threshold_calc,
        mapq_threshold,
        remove_duplicates,
        temp_dir,
        start_position,
        singlecell_mode,  # Add singlecell_mode to args
    ) = args

    print(f"Filtering parameters: consecutive_flag={consecutive_flag}, repeat_threshold_calc={repeat_threshold_calc}")
    print(f"Patterns (forward): {patterns_regex_forward.pattern}")
    print(f"Patterns (reverse): {patterns_regex_reverse.pattern}")

    region_name = "unmapped"
    temp_bam = os.path.join(temp_dir, f"region_{region_name}_filtered.bam")

    gc_content = {}
    read_counts = {"unmapped": {"unmapped": 0}}
    filtered_read_count = 0
    barcode_counts = defaultdict(int) if singlecell_mode else None
    total_reads_processed = 0

    open_mode = "rb" if bam_path.endswith(".bam") else "rc"
    with pysam.AlignmentFile(bam_path, open_mode) as bamfile:
        with pysam.AlignmentFile(temp_bam, "wb", template=bamfile) as filtered_file:
            try:
                # Bam uses seek, Cram needs last reference
                if open_mode == "rb":
                    if start_position > 0:
                        bamfile.seek(start_position)
                    else:
                        bamfile.reset()
                    fetch_reads = bamfile.fetch(until_eof=True)
                else:
                    fetch_reads = bamfile.fetch(contig="*")

                for read in fetch_reads:
                    try:
                        if not read.is_unmapped:
                            continue
                        if (
                                read.is_secondary
                                or read.is_supplementary
                                or (remove_duplicates and read.is_duplicate)
                        ):
                            continue

                        sequence = read.query_sequence
                        if sequence is None:
                            continue

                        read_length = len(sequence)

                        n_count = sequence.count("N")
                        read_length_no_N = read_length - n_count
                        if read_length_no_N > 0 and (n_count / read_length) <= 0.2:
                            gc_percent = int(
                                round(
                                    (sequence.count("C") + sequence.count("G"))
                                    / read_length_no_N
                                    * 100
                                )
                            )
                            gc_content[gc_percent] = gc_content.get(gc_percent, 0) + 1

                        read_counts["unmapped"]["unmapped"] += 1

                        # Barcode counting
                        if singlecell_mode:
                            bc = read.get_tag("CB") if read.has_tag("CB") else None
                            if bc:
                                barcode_counts[bc] += 1

                        # Check if it's a telomere read
                        if is_telomere_read(
                                consecutive_flag,
                                patterns_regex_forward,
                                patterns_regex_reverse,
                                sequence,
                                repeat_threshold_calc,
                        ):
                            filtered_file.write(read)
                            filtered_read_count += 1
                    except Exception as e:
                        print(f"Error processing unmapped read: {e}")
                        traceback.print_exc()
            except (MemoryError, OSError) as e:
                print(f"Error processing unmapped reads: {e}")
                traceback.print_exc()
            except Exception as e:
                print(f"Unexpected error in unmapped reads: {e}")
                traceback.print_exc()

    print(f"Total reads processed: {total_reads_processed}")
    print(f"Total telomeric reads found: {filtered_read_count}")

    return {
        "region": "unmapped",
        "gc_content": gc_content,
        "read_counts": read_counts,
        "temp_bam": temp_bam,
        "filtered_read_count": filtered_read_count,
        "barcode_counts": dict(barcode_counts) if barcode_counts is not None else {},
    }


@measure_time
def parallel_filter_telomere_reads(
        bam_path,
        out_dir,
        pid,
        sample,
        repeat_threshold_calc,
        mapq_threshold,
        repeats,
        consecutive_flag,
        remove_duplicates,
        band_file=None,
        num_processes=None,
        singlecell_mode=None,
        fast_mode=False,
):
    """
    Region-based parallel implementation of telomere read filtering with improved unmapped reads handling.
    Temporary files are stored in the output directory and cleaned up after processing.
    """

    # Create a temporary directory within the output directory
    temp_dir = os.path.join(out_dir, f"temp_{pid}")
    os.makedirs(temp_dir, exist_ok=True)

    try:
        if num_processes is None:
            num_processes = mp.cpu_count()
        print(f"Using {num_processes} cores")

        # Initialize chromosome and band data
        open_mode = "rb" if bam_path.endswith(".bam") else "rc"
        with pysam.AlignmentFile(bam_path, mode=open_mode) as bamfile:
            band_info = (
                initialize_chromosome_and_band_data(bamfile, band_file)
                if band_file
                else None
            )
            references = bamfile.references
            lengths = bamfile.lengths

        # Compile regex patterns
        patterns_regex_forward, patterns_regex_reverse = compile_patterns(repeats)

        results = []
        max_position = 0
        total_filtered_reads = 0  # Initialize total filtered reads counter
        barcode_counts_merged = defaultdict(int)
        if singlecell_mode:
            print("Single-cell mode activated: Barcode counting enabled.")

        # FAST MODE: Only process unmapped reads, skip region-based processing
        if fast_mode:
            print("Fast mode: Skipping region-based processing, only processing unmapped reads.")
            unmapped_args = (
                bam_path,
                patterns_regex_forward,
                patterns_regex_reverse,
                consecutive_flag,
                repeat_threshold_calc,
                mapq_threshold,
                remove_duplicates,
                temp_dir,
                0,  # start_position
                singlecell_mode,
            )
            try:
                try:
                    unmapped_result = process_unmapped_reads(unmapped_args)
                    if unmapped_result is not None:
                        results.append(unmapped_result)
                        total_filtered_reads += unmapped_result["filtered_read_count"]
                        for bc, count in unmapped_result.get("barcode_counts", {}).items():
                            barcode_counts_merged[bc] += count
                        print(
                            f"Unmapped reads processing completed - {unmapped_result['filtered_read_count']} reads filtered"
                        )
                except BrokenProcessPool as bpe:
                    print("ERROR: BrokenProcessPool encountered during unmapped region processing.")
                    print("This may be due to memory issues, corrupted BAM/CRAM, or a bug in process_unmapped_reads.")
                    print(f"Details: {bpe}")
                    raise
            except Exception as e:
                print(f"Error processing unmapped reads: {e}")
        else:
            # First process all mapped regions
            regions = [(ref, 1, length) for ref, length in zip(references, lengths)]
            with ProcessPoolExecutor(max_workers=num_processes) as executor:
                futures = []
                print(f"Processing {len(regions)} regions")
                for region_info in regions:
                    args = (
                        bam_path,
                        region_info,
                        patterns_regex_forward,
                        patterns_regex_reverse,
                        consecutive_flag,
                        repeat_threshold_calc,
                        mapq_threshold,
                        remove_duplicates,
                        band_info,
                        temp_dir,
                        singlecell_mode,
                    )
                    futures.append(executor.submit(process_region, args))
                try:
                    for future in as_completed(futures):
                        try:
                            result = future.result()
                            # Collect results and track the maximum file position
                            if result is not None:
                                results.append(result)
                                max_position = max(max_position, result.get("last_position", 0))
                                total_filtered_reads += result["filtered_read_count"]
                                for bc, count in result.get("barcode_counts", {}).items():
                                    barcode_counts_merged[bc] += count
                                print(
                                    f"Region {result['region']} completed - {result['filtered_read_count']} reads filtered"
                                )
                        except Exception as e:
                            print(f"Error in region processing: {e}")
                except BrokenProcessPool as bpe:
                    print("ERROR: BrokenProcessPool encountered. One of the subprocesses crashed.")
                    print("This may be due to memory issues, corrupted BAM/CRAM, or a bug in process_region.")
                    print(f"Details: {bpe}")
                    # Optionally, print which regions were completed
                    completed_regions = [r.get("region") for r in results if "region" in r]
                    print(f"Regions completed before crash: {completed_regions}")
                    # Optionally, re-raise or exit
                    raise

            # Process unmapped reads
            print(f"Processing unmapped reads from position: {max_position}")
            unmapped_args = (
                bam_path,
                patterns_regex_forward,
                patterns_regex_reverse,
                consecutive_flag,
                repeat_threshold_calc,
                mapq_threshold,
                remove_duplicates,
                temp_dir,
                max_position,
                singlecell_mode,
            )

            try:
                try:
                    unmapped_result = process_unmapped_reads(unmapped_args)
                    if unmapped_result is not None:
                        results.append(unmapped_result)
                        total_filtered_reads += unmapped_result["filtered_read_count"]
                        for bc, count in unmapped_result.get("barcode_counts", {}).items():
                            barcode_counts_merged[bc] += count
                        print(
                            f"Unmapped reads processing completed - {unmapped_result['filtered_read_count']} reads filtered"
                        )
                except BrokenProcessPool as bpe:
                    print("ERROR: BrokenProcessPool encountered during unmapped region processing.")
                    print("This may be due to memory issues, corrupted BAM/CRAM, or a bug in process_unmapped_reads.")
                    print(f"Details: {bpe}")
                    raise
            except Exception as e:
                print(f"Error processing unmapped reads: {e}")

        # Print total filtered reads
        print(f"\nTotal number of reads that passed filtering: {total_filtered_reads}")

        # Merge results
        gc_content_merged = defaultdict(int)
        read_counts_merged = defaultdict(lambda: defaultdict(int))

        for result in results:
            # Merge GC content
            for gc, count in result["gc_content"].items():
                gc_content_merged[gc] += count

            # Merge read counts
            for ref, bands in result["read_counts"].items():
                for band, count in bands.items():
                    read_counts_merged[ref][band] += count

        # Merge BAM files
        output_bam = os.path.join(out_dir, f"{pid}_filtered.bam")
        temp_bams = [
            result["temp_bam"]
            for result in results
            if os.path.exists(result["temp_bam"])
               and os.path.getsize(result["temp_bam"]) > 0
        ]

        if temp_bams:
            if len(temp_bams) == 1:
                shutil.copy(temp_bams[0], output_bam)
            else:
                # Merge BAMs
                pysam.merge("-f", output_bam, *temp_bams)

            # Sort by name and index the filtered file
            pysam.index(output_bam)
            pysam.sort(
                "-n",
                "-o",
                os.path.join(out_dir, f"{pid}_filtered_name_sorted.bam"),
                output_bam,
            )

            # Clean up temporary BAM files
            for temp_bam in temp_bams:
                try:
                    os.remove(temp_bam)
                except OSError as e:
                    print(f"Error removing temporary BAM file {temp_bam}: {e}")
        else:
            print("Warning: No reads passed filtering criteria")

        if singlecell_mode and isinstance(barcode_counts_merged, dict) and not barcode_counts_merged:
            print(
                "Warning: single-cell mode is active but no barcodes were found. This may indicate an error in barcode extraction or input data.")

        # Write results to files
        write_output(
            out_dir,
            pid,
            sample,
            dict(gc_content_merged),
            dict(read_counts_merged),
            band_info,
            barcode_counts=dict(barcode_counts_merged) if singlecell_mode else None,
        )

    finally:
        # Clean up temporary directory
        try:
            shutil.rmtree(temp_dir)
        except OSError as e:
            print(f"Error removing temporary directory {temp_dir}: {e}")

    return None

def main():
    pass

if __name__ == "__main__":
    main()