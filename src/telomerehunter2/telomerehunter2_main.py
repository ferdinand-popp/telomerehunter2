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


import argparse
import multiprocessing
import os
import re
import shutil
import stat
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pysam

from telomerehunter2 import (TVR_context, TVR_context_summary_tables,
                             TVR_screen, estimate_telomere_content,
                             filter_telomere_reads, get_repeat_threshold,
                             get_summed_intratelomeric_read_length, merge_pdfs,
                             normalize_TVR_counts, plot_functions,
                             repeat_frequency_intratelomeric,
                             sort_telomere_reads)

source_directory = os.path.abspath(os.path.dirname(__file__))


def run_sample(
    sample_name,
    sample_bam,
    outdir_sample,
    sample_id,
    read_length,
    repeat_threshold_calc,
    repeat_threshold_str,
    filter_telomere_reads_flag,
    args,
):
    if getattr(args, "singlecell_mode", False):
        print(f"[TH2] Running sample '{sample_name}' in single-cell mode.")

    if filter_telomere_reads_flag:
        print(f"------ {sample_name}: started filtering telomere reads ------")
        filter_telomere_reads.parallel_filter_telomere_reads(
            bam_path=sample_bam,
            out_dir=outdir_sample,
            pid=args.pid,
            sample=sample_id,
            repeat_threshold_calc=repeat_threshold_calc,
            mapq_threshold=args.mapq_threshold,
            repeats=args.repeats,
            consecutive_flag=args.consecutive,
            remove_duplicates=args.remove_duplicates,
            subsample=args.subsample,
            band_file=args.banding_file,
            num_processes=args.cores,
            singlecell_mode=getattr(args, "singlecell_mode", False),
        )

    if args.sort_telomere_reads_flag:
        print("------ " + sample_name + ": started sorting telomere reads ------")
        sort_telomere_reads.sort_telomere_reads(
            input_dir=outdir_sample,
            out_dir=outdir_sample,
            band_file=args.banding_file,
            pid=args.pid,
            mapq_threshold=args.mapq_threshold,
            repeats=args.repeats,
        )

        # get a table with repeat frequencies per intratelomeric read
        repeat_frequency_intratelomeric.repeat_frequency_intratelomeric(
            input_path=outdir_sample,
            out_dir=outdir_sample,
            pid=args.pid,
            repeats=args.repeats,
        )

    if args.estimate_telomere_content_flag:
        print("------ " + sample_name + ": started estimating telomere content ------")

        estimate_telomere_content.get_gc_content_distribution(
            bam_file=os.path.join(
                outdir_sample, f"{args.pid}_filtered_intratelomeric.bam"
            ),
            out_dir=outdir_sample,
            pid=f"{args.pid}_intratelomeric",
            sample=sample_id,
            remove_duplicates=args.remove_duplicates,
        )

        estimate_telomere_content.estimate_telomere_content(
            input_dir=outdir_sample,
            out_dir=outdir_sample,
            sample=sample_id,
            read_length=read_length,
            repeat_threshold_str=repeat_threshold_str,
            pid=args.pid,
            repeat_threshold_set=args.repeat_threshold_set,
            per_read_length=args.per_read_length,
            gc_lower=args.gc_lower,
            gc_upper=args.gc_upper,
        )

    if args.TVR_screen_flag:
        print("------ " + sample_name + ": started TVR screen ------")

        # screen for TVRs
        TVR_screen.tvr_screen(
            main_path=outdir_sample,
            sample=sample_id,
            telomere_pattern="(.{3})(GGG)",
            min_base_quality=20,
            pid=args.pid,
            repeat_threshold_set=args.repeat_threshold_set,
        )

        # get summed intratelomeric read length
        get_summed_intratelomeric_read_length.summed_intratelomeric_read_length(
            main_path=outdir_sample, sample=sample_id, pid=args.pid
        )

    if args.TVR_context_flag:
        print("------ " + sample_name + ": started TVR context ------")

        # get context of selected TVRs
        # Todo not possible for flexible hexamers yet due to GGG and number of positions interesting
        for TVR in args.TVRs_for_context:
            TVR_context.tvr_context(
                main_path=outdir_sample,
                sample=sample_id,
                pattern=TVR,
                min_base_quality=20,
                telomere_pattern="GGG",
                tel_file="filtered_intratelomeric",
                pid=args.pid,
                context_before=args.bp_context,
                context_after=args.bp_context,
                repeat_threshold_set=args.repeat_threshold_set,
            )


def print_copyright_message():
    print("\n")
    print(
        "\tTelomereHunter2 Copyright 2024 Ferdinand Popp, Lina Sieverling, Philip Ginsbach, Chen Hong, Lars Feuerbach"
    )
    print("\tThis program comes with ABSOLUTELY NO WARRANTY.")
    print("\tThis is free software, and you are welcome to redistribute it")
    print(
        "\tunder certain conditions. For details see the GNU General Public License v3.0"
    )
    print(
        "\tin the license copy received with TelomereHunter2 or <http://www.gnu.org/licenses/>."
    )
    print("\n")
    print("TelomereHunter2 1.0.0")
    print("\n")


def file_exists(parser, path):
    if path is None:
        return None
    if not os.path.exists(path):
        parser.error("The file %s does not exist!" % path)
    return path


def check_banding_file(banding_file, outdir=None):
    """
    Tests if the supplied tsv file has format like: chr1	0	2300000	p36.33	...
    """
    try:
        # Read the file (assuming tab-separated or CSV)
        try:
            df = pd.read_csv(banding_file, sep="\t", header=None)
        except:
            df = pd.read_csv(banding_file, header=None)
            raise ValueError(f"Banding file is not a tsv file.")

        # Check if file has at least 4 columns
        if df.shape[1] < 4:
            raise ValueError(
                f"File has fewer than 4 columns. Found {df.shape[1]} columns."
            )

        # Check column types for the first 4 columns
        type_errors = []

        # Check first column (should be string)
        if not all(isinstance(x, str) for x in df[0].dropna()):
            type_errors.append("Column 1 should contain strings (chromosome names)")

        # Check columns 2 and 3 (should be integers)
        for col_idx in [1, 2]:
            # Convert to numeric if possible
            df[col_idx] = pd.to_numeric(df[col_idx], errors="coerce")
            if df[col_idx].isna().any():
                type_errors.append(f"Column {col_idx + 1} should contain integers")

        # Check column 4 (should be string)
        if not all(isinstance(x, str) for x in df[3].dropna()):
            type_errors.append("Column 4 should contain strings (band names)")

        # If any type errors found, raise exception
        if type_errors:
            error_msg = "Invalid column types:\n" + "\n".join(
                f"- {err}" for err in type_errors
            )
            raise ValueError(error_msg)

        # Check if first column contains only unique values (chromosome-level data)
        chromosome_level_only = df[0].nunique() == len(df)

        if chromosome_level_only:
            print("!!! Only chromosome-level data found, no band information. !!!")
            print("The banding file should be discarded.")

            # Generate band information if requested
            subs_cytobands_df = generate_banding_file(df)
            print(
                "\nGenerated substitution cytoband information (5% from chromosome ends):"
            )
            print(subs_cytobands_df.head())
            print(f"Total generated bands: {len(subs_cytobands_df)}")
            generated_banding_file_path = (
                f"{outdir}/generated_bandingfile_5percent_ends.txt"
            )
            subs_cytobands_df.to_csv(generated_banding_file_path, sep="\t", index=False)
            print(
                f"Saved cytoband information (5% from chromosome ends): {generated_banding_file_path}"
            )

            return generated_banding_file_path
        else:
            print("Band-level information found. File is valid for detailed analysis.")
            return banding_file

    except Exception as e:
        print(f"ERROR: {str(e)}")
        sys.exit(1)


def generate_banding_file(df):
    """
    Generate band information for chromosomes by adding bands at 5% from the ends

    Args:
        df: DataFrame with chromosome-level data

    Returns:
        DataFrame with generated band information
    """
    bands_df = []

    for _, row in df.iterrows():
        chrom = row[0]
        start = row[1]
        end = row[2]

        # Calculate positions at 5% from each end
        band_width = end - start
        p_terminal = start + int(0.05 * band_width)
        q_terminal = end - int(0.05 * band_width)

        # Create bands: p-terminal, middle, q-terminal
        bands_df.append([chrom, start, p_terminal, f"{chrom}pter"])
        bands_df.append([chrom, p_terminal, q_terminal, f"{chrom}pq"])
        bands_df.append([chrom, q_terminal, end, f"{chrom}qter"])

    return pd.DataFrame(bands_df)


def get_bamfile_reference_genome(bam_file_path):
    """
    Extracts a specific string from the second line of the BAM file header.
    """
    try:
        mode = "rb" if bam_file_path.endswith(".bam") else "rc"
        with pysam.AlignmentFile(bam_file_path, mode=mode) as bam_file:
            header = bam_file.header.to_dict()
            reference_genome_string = header["SQ"][0]
            print(f"Input file is based on reference genome: {reference_genome_string}")
    except KeyError:
        return "Error when extracting reference genome from bam file."


def validate_args(args):
    # check which bam files were specified
    if not args.plot_mode and not any([args.tumor_bam, args.control_bam]):
        raise ValueError(
            "argument -ibt/--inputBamTumor or -ibc/--inputBamControl or -icc or -ict is required"
        )

    if not args.tumor_bam:
        args.tumor_flag = False
        print(
            "Tumor BAM/CRAM file was not specified. Only running control BAM/CRAM file."
        )
    elif not args.control_bam:
        args.control_flag = False
        print(
            "Control BAM/CRAM file was not specified. Only running tumor BAM/CRAM file."
        )

    # Check band file format
    if args.banding_file:
        if args.tumor_bam or args.control_bam:
            get_bamfile_reference_genome(args.tumor_bam or args.control_bam)

        # check banding file and if only chromosomes, create 5% bands
        args.banding_file = check_banding_file(args.banding_file, args.outdir)
    else:
        print(
            "No banding file supplied, so all banding steps will be skipped and specific plots omitted."
        )

    # check if repeats only contains ACGT and are between 4 and 9 bases
    for repeat in args.repeats:
        if len(repeat) < 4 or len(repeat) > 9:
            raise ValueError(
                "argument -r/--repeats should be between 4 and 9 bases long."
            )
        x = re.search(r"[^ACGT]", repeat)
        if x is not None:
            raise ValueError(
                "argument -r/--repeats should only contain the letters ACGT."
            )

    for repeat in args.TVRs_for_context:
        x = re.search(r"[^ACGT]", repeat)
        if x is not None:
            raise ValueError(
                "argument -rc/--repeatsContext should only contain the letters ACGT."
            )

    # check if bp for sequencing context is divisible by 6 or base telomere sequence given
    if args.bp_context % len(args.repeats[0]) != 0:
        raise ValueError(
            f"argument -bp/--bpContext must be a multiple of length of first telomere repeat sequence /"
            f'{len(args.repeats[0])}. By default "TTAGGG" so multiple of 6 (e.g. 6, 12, 18, ...).'
        )

    if args.mapq_threshold < 0 or args.mapq_threshold > 40:
        raise ValueError(
            "argument -mqt/--mappingQualityThreshold must be an integer between 0 and 40."
        )

    if args.gc_lower < 0 or args.gc_lower > 100:
        raise ValueError(
            "argument -gc1/--lowerGC must be an integer between 0 and 100."
        )

    if args.gc_upper < 0 or args.gc_upper > 100:
        raise ValueError(
            "argument -gc2/--upperGC must be an integer between 0 and 100."
        )

    if args.gc_lower >= args.gc_upper:
        raise ValueError(
            "argument -gc1/--lowerGC must be less than argument -gc2/--upperGC."
        )

    if args.subsample:
        if args.subsample < 0 or args.subsample > 1:
            raise argparse.ArgumentTypeError(
                f"Subsample value must be between 0 and 1. Got {args.subsample} instead."
            )


def validate_plotting_options(args):
    if args.plotNone and (
        args.plotChr
        or args.plotFractions
        or args.plotTelContent
        or args.plotGC
        or args.plotRepeatFreq
        or args.plotTVR
        or args.plotSingleton
    ):
        raise ValueError(
            "argument -pno/--plotNone should not be specified when other plotting options are selected."
        )

    # can be str or list so change to list
    if isinstance(args.plotFileFormat, str):
        args.plotFileFormat = [args.plotFileFormat]
    if args.plotFileFormat == "all":
        args.plotFileFormat = ["pdf", "png", "svg"]

    # if no plotting options are selected: plot all diagrams.
    if (
        not args.plotChr
        and not args.plotFractions
        and not args.plotTelContent
        and not args.plotGC
        and not args.plotRepeatFreq
        and not args.plotTVR
        and not args.plotSingleton
        and not args.plotNone
    ):
        args.plotChr = True
        args.plotFractions = True
        args.plotTelContent = True
        args.plotGC = True
        args.plotRepeatFreq = True
        args.plotTVR = True
        args.plotSingleton = True
    # no banding info remove plots
    if not args.banding_file:
        args.plotChr = False
        args.plotFractions = False


def set_execution_flags(args):
    args.filter_telomere_reads_flag = not args.noFiltering
    args.sort_telomere_reads_flag = args.estimate_telomere_content_flag = True
    args.TVR_screen_flag = args.TVR_context_flag = True

    if args.repeats[0] != "TTAGGG":
        print(
            "The first repeat to screen for is not canonical TTAGGG! Please double check."
        )
        args.TVR_screen_flag = args.TVR_context_flag = args.plotTVR = (
            args.plotSingleton
        ) = False


def parse_command_line_arguments():
    # Cmd line input.
    parser = argparse.ArgumentParser(
        description="Estimation of telomere content from WGS data of a tumor and/or a control sample.",
        epilog="Contact Ferdinand Popp (f.popp@dkfz-heidelberg.de) for questions and support.",
    )

    input_files_group = parser.add_argument_group("Input Files")
    input_files_group.add_argument(
        "-p",
        "--pid",
        type=str,
        dest="pid",
        required=True,
        help="Sample name used in output files and diagrams (required).",
    )
    input_files_group.add_argument(
        "-o",
        "--outPath",
        type=str,
        dest="parent_outdir",
        required=True,
        help="Path to the output directory into which all results are written.",
    )
    input_files_group.add_argument(
        "-ibt",
        "--inputBamTumor",
        type=lambda x: file_exists(parser, x),
        dest="tumor_bam",
        help="Path to the indexed input BAM/CRAM file of the tumor sample.",
    )
    input_files_group.add_argument(
        "-ibc",
        "--inputBamControl",
        type=lambda x: file_exists(parser, x),
        dest="control_bam",
        help="Path to the indexed input BAM/CRAM file of the control sample.",
    )
    input_files_group.add_argument(
        "-b",
        "--bandingFile",
        type=lambda x: file_exists(parser, x),
        dest="banding_file",
        default=None,
        help="Optional: Path to a tab-separated file with information on chromosome banding. \
                                            The first four columns of the table have to contain the chromosome name, \
                                            the start and end position and the band name. The table should not have a header. \
                                            If no banding file is specified, scripts and plots will omit banding \
                                            information. Reference files are included for hg19 and hg38 in \
                                            telomerehunter2/cytoband_file. For example hg19_cytoBand.txt",
    )

    threshold_filtering_group = parser.add_argument_group(
        "Thresholds and Filtering Options"
    )
    threshold_filtering_group.add_argument(
        "-nf",
        "--noFiltering",
        dest="noFiltering",
        action="store_true",
        help="If the filtering step of TelomereHunter has already been run previously, skip this step.",
    )
    threshold_filtering_group.add_argument(
        "-rt",
        "--repeatThreshold",
        type=int,
        dest="repeat_threshold_set",
        help="The number of repeats needed for a read to be classified as telomeric. Minimum is 4. Values below 4 will be set to 4 automatically.",
    )
    threshold_filtering_group.add_argument(
        "-mqt",
        "--mappingQualityThreshold",
        type=int,
        dest="mapq_threshold",
        default=8,
        help="The mapping quality needed for a read to be considered as mapped (default = 8).",
    )
    threshold_filtering_group.add_argument(
        "-gc1",
        "--lowerGC",
        type=int,
        dest="gc_lower",
        default=48,
        help="Lower limit used for GC correction of telomere content.",
    )
    threshold_filtering_group.add_argument(
        "-gc2",
        "--upperGC",
        type=int,
        dest="gc_upper",
        default=52,
        help="Upper limit used for GC correction of telomere content.",
    )
    threshold_filtering_group.add_argument(
        "-con",
        "--consecutive",
        dest="consecutive",
        action="store_true",
        help="Search for consecutive repeats.",
    )
    threshold_filtering_group.add_argument(
        "-d",
        "--removeDuplicates",
        dest="remove_duplicates",
        action="store_true",
        help="Reads marked as duplicates in the input bam file(s) are removed in the filtering step.",
    )
    threshold_filtering_group.add_argument(
        "-rl",
        "--perReadLength",
        dest="per_read_length",
        action="store_true",
        help="Repeat threshold is set per 100 bp read length.",
    )
    threshold_filtering_group.add_argument(
        "--min-reads-per-barcode",
        type=int,
        dest="min_reads_per_barcode",
        default=10000,
        help="Minimum reads per barcode for single-cell summary only needed by telomerehunter2_sc.py not the bulk application (default: 10000).",
    )

    repeats_context_group = parser.add_argument_group("Repeats and Context Options")
    repeats_context_group.add_argument(
        "-r",
        "--repeats",
        nargs="+",
        dest="repeats",
        type=str,
        help="Base sequences to inspect like TTAGGG and its TVRs. First sequence is /"
        "used as base telomeric sequence. Default: ['TTAGGG', 'TGAGGG', 'TCAGGG', /"
        "'TTCGGG', 'TTGGGG'].",
        default=["TTAGGG", "TGAGGG", "TCAGGG", "TTGGGG", "TTCGGG", "TTTGGG"],
    )
    repeats_context_group.add_argument(
        "-rc",
        "--repeatsContext",
        nargs="+",
        dest="TVRs_for_context",
        type=str,
        default=[
            "TCAGGG",
            "TGAGGG",
            "TTGGGG",
            "TTCGGG",
            "TTTGGG",
            "ATAGGG",
            "CATGGG",
            "CTAGGG",
            "GTAGGG",
            "TAAGGG",
        ],
        help="List of telomere variant repeats for which to analyze the sequence context.",
    )
    repeats_context_group.add_argument(
        "-bp",
        "--bpContext",
        type=int,
        dest="bp_context",
        default=18,
        help="Number of base pairs on either side of the telomere variant repeat to investigate. Please use a number that is divisible by 6.",
    )

    running_group = parser.add_argument_group("Running Options")
    running_group.add_argument(
        "--subsample",
        dest="subsample",
        default=False,
        type=float,
        help="Specify a subsample fraction (float between 0 and 1). If not provided, defaults to False.",
    )
    running_group.add_argument(
        "-c",
        "--cores",
        dest="cores",
        default=None,
        type=int,
        help="Specify the number of cores. If not provided, defaults to None and uses max cores.",
    )
    running_group.add_argument(
        "-pl",
        "--parallel",
        dest="run_parallel",
        action="store_true",
        help="The filtering, sorting, and estimating steps of the tumor and control sample are run in parallel. This will speed up the computation time of telomerehunter.",
    )
    running_group.add_argument(
        "--plot_mode",
        action="store_true",
        dest="plot_mode",
        help="Run only plotting steps on an existing analysis folder",
    )
    running_group.add_argument(
        "-sc",
        "--singlecell_mode",
        dest="singlecell_mode",
        action="store_true",
        help="Enable barcode counting/output for single-cell BAMs (CB tag). Default: auto-detect if CB tags present.",
    )

    plotting_group = parser.add_argument_group("Plotting Options")
    plotting_group.add_argument(
        "-pno",
        "--plotNone",
        dest="plotNone",
        action="store_true",
        default=False,
        help="Do not make any diagrams.",
    )
    plotting_group.add_argument(
        "-p1",
        "--plotChr",
        dest="plotChr",
        action="store_true",
        help="Make diagrams with telomeric reads mapping to each chromosome.",
    )
    plotting_group.add_argument(
        "-p2",
        "--plotFractions",
        dest="plotFractions",
        action="store_true",
        help="Make a diagram with telomere reads in each fraction.",
    )
    plotting_group.add_argument(
        "-p3",
        "--plotTelContent",
        dest="plotTelContent",
        action="store_true",
        help="Make a diagram with the gc corrected telomere content in the analyzed samples.",
    )
    plotting_group.add_argument(
        "-p4",
        "--plotGC",
        dest="plotGC",
        action="store_true",
        help="Make a diagram with GC content distributions in all reads and in intratelomeric reads.",
    )
    plotting_group.add_argument(
        "-p5",
        "--plotRepeatFreq",
        dest="plotRepeatFreq",
        action="store_true",
        help="Make histograms of the repeat frequencies per intratelomeric read.",
    )
    plotting_group.add_argument(
        "-p6",
        "--plotTVR",
        dest="plotTVR",
        action="store_true",
        help="Make plots for telomere variant repeats.",
    )
    plotting_group.add_argument(
        "-p7",
        "--plotSingleton",
        dest="plotSingleton",
        action="store_true",
        help="Make plots for singleton telomere variant repeats.",
    )
    plotting_group.add_argument(
        "-prc",
        "--plotRevCompl",
        dest="plotRevCompl",
        action="store_true",
        help="Distinguish between forward and reverse complement telomere repeats in diagrams.",
    )
    plotting_group.add_argument(
        "-pff",
        "--plotFileFormat",
        dest="plotFileFormat",
        default="pdf",
        choices=["pdf", "png", "svg", "all"],
        help="File format of output diagrams. Choose from pdf (default), png, svg or all (pdf, png and svg).",
    )

    # Create dict from parser args and save keys as variables
    args = parser.parse_args()

    # Set flags based on input arguments
    args.tumor_flag = bool(args.tumor_bam)
    args.control_flag = bool(args.control_bam)

    # Create output directory path
    args.outdir = os.path.join(args.parent_outdir, args.pid)

    # Validate arguments
    validate_args(args)
    validate_plotting_options(args)

    # Set execution flags
    set_execution_flags(args)

    # Input validation for repeat_threshold_set
    if args.repeat_threshold_set is not None and args.repeat_threshold_set < 4:
        print(
            "Warning: repeatThreshold was set below 4. Setting to minimum value of 4."
        )
        args.repeat_threshold_set = 4

    return args


def has_filtering_output(outdir, pid, sample_id):
    sample_dir = os.path.join(outdir, f"{sample_id}_TelomerCnt_{pid}")
    files = [
        os.path.join(sample_dir, f"{pid}_filtered.bam"),
        os.path.join(sample_dir, f"{pid}_filtered_name_sorted.bam"),
        os.path.join(sample_dir, f"{pid}_readcount.tsv"),
        os.path.join(sample_dir, f"{pid}_{sample_id}_gc_content.tsv"),
    ]
    return all(os.path.exists(file) for file in files)


def check_and_prompt_filtering(
    filter_telomere_reads_flag, tumor_flag, control_flag, outdir, pid
):
    tumor_output_missing = (
        not has_filtering_output(outdir, pid, "tumor") if tumor_flag else False
    )
    control_output_missing = (
        not has_filtering_output(outdir, pid, "control") if control_flag else False
    )

    # Both missing
    if tumor_output_missing and control_output_missing:
        print("Running filtering step as no output files are present")
        # Alert wrong flag use!
        if filter_telomere_reads_flag:
            print(
                "!Denied skipping the filtering with --noFiltering, as no output files were present!"
            )
        filter_T, filter_C = True, True
    # One missing due to only one run or only one present
    elif control_output_missing or tumor_output_missing:
        print("Running filtering step as not all output files are already present")
        if filter_telomere_reads_flag:
            print(
                "!Denied skipping the filtering with --noFiltering, as not all output files were present!"
            )
        filter_T, filter_C = tumor_output_missing, control_output_missing
    # Both present
    else:
        if filter_telomere_reads_flag:
            print(
                "Filtering samples, but output files are already present. Consider using the --noFiltering flag"
            )
            filter_T, filter_C = tumor_flag, control_flag
        else:
            print("Skipped filtering samples with --noFiltering flag")
            filter_T, filter_C = False, False

    return filter_T, filter_C


def get_repeat_threshold_from_summary(args):
    summary_file = os.path.join(args.outdir, f"{args.pid}_summary.tsv")
    try:
        # Read the summary file
        summary = pd.read_csv(summary_file, sep="\t")

        # Ensure the column exists
        if "repeat_threshold_used" in summary.columns:
            # Get unique values in the column
            unique_thresholds = summary["repeat_threshold_used"].dropna().unique()

            # Check the number of unique values
            if len(unique_thresholds) == 1:
                return str(int(unique_thresholds[0]))
            elif len(unique_thresholds) > 1:
                return "n"
        else:
            print(f"Error: Column 'repeat_threshold_used' not found in {summary_file}.")
            return "n"

    except Exception as e:
        print(f"Error reading file {summary_file}: {e}")
        return "n"


def combine_summary_files(outdir, pid, tumor_flag, control_flag):
    # Define the path for the combined summary file
    summary_path = os.path.join(outdir, f"{pid}_summary.tsv")

    tumor_summary_path = os.path.join(
        outdir, f"tumor_TelomerCnt_{pid}", f"{pid}_tumor_summary.tsv"
    )
    control_summary_path = os.path.join(
        outdir, f"control_TelomerCnt_{pid}", f"{pid}_control_summary.tsv"
    )

    # Check if both tumor and control flags are True
    if tumor_flag and control_flag:
        # Copy the tumor summary file to the combined summary file
        shutil.copyfile(tumor_summary_path, summary_path)

        # Append the last line of the control summary file to the combined summary file
        with open(control_summary_path, "r") as control_file, open(
            summary_path, "a"
        ) as combined_file:
            last_line = control_file.readlines()[-1]
            combined_file.write(last_line)

    # Check if only the tumor flag is True
    elif tumor_flag:
        # Copy the tumor summary file to the combined summary file
        shutil.copyfile(tumor_summary_path, summary_path)
    # Check if only the control flag is True
    elif control_flag:
        # Copy the control summary file to the combined summary file
        shutil.copyfile(control_summary_path, summary_path)


def check_and_create_index(bam_path):
    # make path absolute
    bam_path = os.path.abspath(bam_path)

    # Check if index file exists
    index_path = bam_path + ".bai"
    if os.path.exists(index_path):
        return_bam = bam_path  # Index file already exists, no need to create
    else:
        # If the directory is writeable
        if os.access(os.path.dirname(bam_path), os.W_OK):
            pysam.index(bam_path)
            return_bam = bam_path
        else:
            # Not writeable needs temp dir
            temp_dir = os.path.abspath(tempfile.mkdtemp("tmp_TH2"))

            # Set read and write permissions for the temporary directory
            os.chmod(temp_dir, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

            temp_bam_path = os.path.join(temp_dir, os.path.basename(bam_path))

            # symbolic link bam file in temp folder
            try:
                os.symlink(bam_path, temp_bam_path)
            except OSError as e:
                print(f"Symbolic link creation failed: {e}")

            pysam.index(temp_bam_path)
            return_bam = temp_bam_path

    return return_bam


def get_read_lengths_and_repeat_thresholds(args, control_bam, tumor_bam):
    """Calculate read lengths and repeat thresholds for tumor and control samples."""
    # Initialize all return variables
    read_lengths_str_control = None
    read_lengths_str_tumor = None
    repeat_thresholds_control = None
    repeat_thresholds_plot = None
    repeat_thresholds_str_control = None
    repeat_thresholds_str_tumor = None
    repeat_thresholds_tumor = None

    # Set default repeat threshold if not provided
    if not args.repeat_threshold_set:
        args.repeat_threshold_set = 6
        args.per_read_length = True
        print(
            "Repeat threshold per 100 bp was not set by the user. Setting it to 6 reads per 100 bp read length."
        )

    if args.per_read_length:
        # Calculate tumor thresholds if needed
        if args.tumor_flag:
            # Get read lengths and calculate thresholds for tumor
            read_lengths_str_tumor, tumor_read_length_counts = (
                get_repeat_threshold.get_read_lengths(tumor_bam)
            )
            print("Calculating repeat threshold for the tumor sample: ")
            repeat_thresholds_tumor, repeat_thresholds_str_tumor = (
                get_repeat_threshold.get_repeat_threshold(
                    read_lengths_str_tumor,
                    tumor_read_length_counts,
                    args.repeat_threshold_set,
                )
            )

        # Calculate control thresholds if needed
        if args.control_flag:
            # Get read lengths and calculate thresholds for control
            read_lengths_str_control, control_read_length_counts = (
                get_repeat_threshold.get_read_lengths(control_bam)
            )
            print("Calculating repeat threshold for the control sample: ")
            repeat_thresholds_control, repeat_thresholds_str_control = (
                get_repeat_threshold.get_repeat_threshold(
                    read_lengths_str_control,
                    control_read_length_counts,
                    args.repeat_threshold_set,
                )
            )

        # Determine which threshold to use for plotting
        if args.tumor_flag and args.control_flag:
            repeat_thresholds_plot = (
                repeat_thresholds_tumor
                if repeat_thresholds_tumor == repeat_thresholds_control
                else "n"
            )
        elif args.tumor_flag:
            repeat_thresholds_plot = repeat_thresholds_tumor
        elif args.control_flag:
            repeat_thresholds_plot = repeat_thresholds_control

    else:
        # Use fixed threshold for all cases
        repeat_thresholds_tumor = args.repeat_threshold_set
        repeat_thresholds_control = args.repeat_threshold_set
        repeat_thresholds_plot = args.repeat_threshold_set
        repeat_thresholds_str_tumor = str(args.repeat_threshold_set)
        repeat_thresholds_str_control = str(args.repeat_threshold_set)

    print(
        f"Repeat Thresholds: Tumor={repeat_thresholds_tumor}, Control={repeat_thresholds_control}"
    )
    print("\n")

    return (
        read_lengths_str_control,
        read_lengths_str_tumor,
        repeat_thresholds_control,
        repeat_thresholds_plot,
        repeat_thresholds_str_control,
        repeat_thresholds_str_tumor,
        repeat_thresholds_tumor,
    )


def run_telomerehunter(
    args,
    control_bam,
    filter_telomere_reads_control,
    filter_telomere_reads_tumor,
    read_lengths_control,
    read_lengths_tumor,
    repeat_thresholds_control,
    repeat_thresholds_str_control,
    repeat_thresholds_str_tumor,
    repeat_thresholds_tumor,
    tumor_bam,
):
    # Determine maximum number of workers based on parallel execution setting
    max_workers = 1 if not args.run_parallel else min(2, multiprocessing.cpu_count())

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        submitted_futures = []

        # run tumor sample for PID
        if args.tumor_flag and (
            filter_telomere_reads_tumor
            or args.sort_telomere_reads_flag
            or args.estimate_telomere_content_flag
        ):
            tumor_sample_id = "tumor"
            tumor_output_dir = os.path.join(
                args.outdir, f"{tumor_sample_id}_TelomerCnt_{args.pid}"
            )
            os.makedirs(tumor_output_dir, exist_ok=True)

            # Submit to run_sample main function
            tumor_future = executor.submit(
                run_sample,
                sample_name="Tumor Sample",
                sample_bam=tumor_bam,
                outdir_sample=tumor_output_dir,
                sample_id=tumor_sample_id,
                read_length=read_lengths_tumor,
                repeat_threshold_calc=repeat_thresholds_tumor,
                repeat_threshold_str=repeat_thresholds_str_tumor,
                filter_telomere_reads_flag=filter_telomere_reads_tumor,
                args=args,
            )

            submitted_futures.append(tumor_future)

            if not args.run_parallel:
                tumor_future.result()

        # run control samples for PID
        if args.control_flag and (
            filter_telomere_reads_control
            or args.sort_telomere_reads_flag
            or args.estimate_telomere_content_flag
        ):
            control_sample_id = "control"
            control_output_dir = os.path.join(
                args.outdir, f"{control_sample_id}_TelomerCnt_{args.pid}"
            )
            os.makedirs(control_output_dir, exist_ok=True)

            # Submit to run_sample main function
            control_future = executor.submit(
                run_sample,
                sample_name="Control Sample",
                sample_bam=control_bam,
                outdir_sample=control_output_dir,
                sample_id=control_sample_id,
                read_length=read_lengths_control,
                repeat_threshold_calc=repeat_thresholds_control,
                repeat_threshold_str=repeat_thresholds_str_control,
                filter_telomere_reads_flag=filter_telomere_reads_control,
                args=args,
            )

            submitted_futures.append(control_future)

            if not args.run_parallel:
                control_future.result()

        # Wait for all futures to complete (only needed in parallel mode)
        if args.run_parallel:
            for future in submitted_futures:
                future.result()


def delete_temp_dir(tumor_bam, control_bam, temp_dir):
    temp_exist = [
        os.path.dirname(path)
        for path in [tumor_bam, control_bam]
        if path is not None and temp_dir in path
    ]

    if temp_exist:
        try:
            # Attempt to remove the folder and its contents (same for both)
            shutil.rmtree(temp_exist[0])
            print(f"Folder '{temp_exist}' successfully deleted.")
        except Exception as e:
            print(f"Error deleting folder '{temp_exist}': {e}")


def summary_log2(main_path, pid):
    # Construct the summary file path
    summary_file = os.path.join(main_path, f"{pid}_summary.tsv")

    # Read the summary table
    summary = pd.read_csv(summary_file, sep="\t", header=0)

    # Check which samples were run
    samples = summary["sample"]

    # If both "tumor" and "control" were run, calculate summary
    if "tumor" in samples.values and "control" in samples.values:
        start_log_index = 13  # tel_content and TRPM
        tumor_row = (
            summary.loc[samples == "tumor"].iloc[:, start_log_index:].fillna(0).values
        )
        control_row = (
            summary.loc[samples == "control"].iloc[:, start_log_index:].fillna(0).values
        )

        log2_values = np.empty_like(
            tumor_row, dtype=float
        )  # Creating an empty array to hold the result
        log2_values.fill(np.nan)
        mask = (control_row != 0.0) & (tumor_row != 0.0)
        log2_values[mask] = np.log2(np.divide(tumor_row[mask], control_row[mask]))

        # Replace inf and -inf with np.nan
        log2_values[np.isinf(log2_values)] = np.nan

        # Create a new row with [pid, "log2(tumor/control)"] and non-log empty columns
        new_row = pd.DataFrame(
            [[pid, "log2(tumor/control)"] + [np.nan] * (len(summary.columns) - 2)],
            columns=summary.columns,
        )

        # Add the calculated log2 values to the new row
        new_row.iloc[0, start_log_index:] = log2_values.flatten()

        # Append the new row to the summary table
        summary = pd.concat([summary, new_row], ignore_index=True)

    # convert cols to right type
    summary = summary.convert_dtypes()

    # Save extended summary table
    summary.to_csv(summary_file, sep="\t", index=False)


def run_plots(args, outdir, pid, repeat_thresholds_plot):
    if not args.plotNone:
        print("------ making plots ------\n")
        # Create plot directory
        os.makedirs(
            os.path.join(outdir, "plots"), exist_ok=True
        )  # all plots there except summary in outdir
        os.makedirs(
            os.path.join(outdir, "html_reports"), exist_ok=True
        )  # plotly folder

        if args.plotChr:
            print("Plot spectrum (plotChr)")
            plot_functions.plot_spectrum(
                outdir,
                pid,
                repeat_thresholds_plot,
                args.consecutive,
                args.mapq_threshold,
                args.banding_file,
                args.plotRevCompl,
                args.plotFileFormat,
            )

        if args.plotFractions:
            print("Plot spectrum summary (plotFractions)")
            plot_functions.plot_spectrum_summary(
                outdir,
                pid,
                repeat_thresholds_plot,
                args.consecutive,
                args.mapq_threshold,
                args.banding_file,
                args.plotRevCompl,
                args.plotFileFormat,
            )

        if args.plotTelContent:
            print("Plot telomere content (plotTelContent)")
            plot_functions.plot_tel_content(
                outdir,
                pid,
                repeat_thresholds_plot,
                args.consecutive,
                args.mapq_threshold,
                args.plotRevCompl,
                args.plotFileFormat,
                args.gc_lower,
                args.gc_upper,
            )

        if args.plotGC:
            print("Plot GC content (plotGC)")
            plot_functions.plot_gc_content(
                outdir, pid, args.plotFileFormat, args.gc_lower, args.gc_upper
            )

        if args.plotRepeatFreq:
            print("Plot intratelomeric repeat frequency (plotRepeatFreq)")
            plot_functions.plot_repeat_frequency_intratelomeric(
                outdir,
                pid,
                repeat_thresholds_plot,
                args.consecutive,
                args.mapq_threshold,
                args.repeats,
                args.plotFileFormat,
            )

        if args.plotTVR:
            print("Plot TVRs (plotTVR)")
            plot_functions.plot_TVR_plot(outdir, pid, args.plotFileFormat)

        if args.plotSingleton:
            print("Plot singletons (plotSingleton)")
            plot_functions.plot_singletons(outdir, pid, args.plotFileFormat)

        # Open to change due to plotly solutions
        # merge PDF plots if PDF files are created
        print("Combining results")
        merge_pdfs.merge_telomere_hunter_results(pid, outdir, args.banding_file)


def main():
    # print welcome
    print_copyright_message()

    # parse and check input parameters
    args = parse_command_line_arguments()
    if getattr(args, "singlecell_mode", False):
        print("Single-cell mode detected. Running full file analysis and barcode-specific insights.")

    # If plot_mode is True, directly run plots and exit
    if args.plot_mode:
        repeat_thresholds_plot = get_repeat_threshold_from_summary(args)

        run_plots(args, args.outdir, args.pid, repeat_thresholds_plot)
        print("Completed plotting.")
        return

    # check if index is present, create it and if not writeable dir then do it in temp dir
    tumor_bam = (
        check_and_create_index(args.tumor_bam) if args.tumor_bam is not None else None
    )
    control_bam = (
        check_and_create_index(args.control_bam)
        if args.control_bam is not None
        else None
    )

    ################################################
    ### check if filtering has already been done ###
    ################################################

    filter_telomere_reads_T, filter_telomere_reads_C = check_and_prompt_filtering(
        args.filter_telomere_reads_flag,
        args.tumor_flag,
        args.control_flag,
        args.outdir,
        args.pid,
    )

    ###############################################
    #### Get Repeat Thresholds and Read Lengths ###
    ###############################################

    (
        read_lengths_control,
        read_lengths_tumor,
        repeat_thresholds_control,
        repeat_thresholds_plot,
        repeat_thresholds_str_control,
        repeat_thresholds_str_tumor,
        repeat_thresholds_tumor,
    ) = get_read_lengths_and_repeat_thresholds(args, control_bam, tumor_bam)

    ###############################
    ## run tumor sample for PID ###
    ###############################

    run_telomerehunter(
        args,
        control_bam,
        filter_telomere_reads_C,
        filter_telomere_reads_T,
        read_lengths_control,
        read_lengths_tumor,
        repeat_thresholds_control,
        repeat_thresholds_str_control,
        repeat_thresholds_str_tumor,
        repeat_thresholds_tumor,
        tumor_bam,
    )

    ##############
    ## Summary ###
    ##############

    # make a combined summary file of tumor and control results
    combine_summary_files(args.outdir, args.pid, args.tumor_flag, args.control_flag)

    # close possible temp dir
    delete_temp_dir(tumor_bam, control_bam, temp_dir="tmp_TH2")

    print("TVR_screen_flag")
    # make a table with normalized TVR counts
    if args.TVR_screen_flag:
        normalize_TVR_counts.normalize_TVR_counts(
            args.outdir, args.pid, ",".join(args.TVRs_for_context)
        )

    print("TVR_context_flag")
    # make TVR context summary tables
    if args.TVR_context_flag:
        TVR_context_summary_tables.tvr_context_singletons_tables(
            args.outdir,
            args.pid,
            args.repeats[0],
            args.bp_context,
            ",".join(args.TVRs_for_context),
        )

    # add log2 ratio to summary file
    if args.tumor_flag and args.control_flag:
        summary_log2(args.outdir, args.pid)

    ###############
    ## Plotting ###
    ###############

    run_plots(args, args.outdir, args.pid, repeat_thresholds_plot)

    print("Completed run.")


if __name__ == "__main__":
    main()
