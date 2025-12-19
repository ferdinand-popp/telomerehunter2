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
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from telomerehunter2 import (
    TVR_context,
    TVR_context_summary_tables,
    TVR_screen,
    estimate_telomere_content,
    filter_telomere_reads,
    get_repeat_threshold,
    get_summed_intratelomeric_read_length,
    merge_pdfs,
    normalize_TVR_counts,
    plot_functions,
    repeat_frequency_intratelomeric,
    sort_telomere_reads,
)
from telomerehunter2.fast_mode import run_fast_mode
from telomerehunter2.utils import (
    check_and_create_index,
    check_and_prompt_filtering,
    combine_summary_files,
    delete_temp_dir,
    file_exists,
    get_repeat_threshold_from_summary,
    print_copyright_message,
    set_execution_flags,
    validate_args,
    validate_plotting_options,
)


def run_sample(
        sample_name,
        sample_bam,
        outdir_sample,
        sample_id,
        read_length,
        repeat_threshold_calc,
        repeat_threshold_str,
        filter_telomere_reads_flag,
        region_cores,
        args,
):
    try:
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
                band_file=args.banding_file,
                num_processes=region_cores,
                singlecell_mode=getattr(args, "singlecell_mode", False),
                barcode_tag=getattr(args, "barcode_tag", "CB"),
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
            print(
                "------ " + sample_name + ": started estimating telomere content ------"
            )

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
                )
        return {"status": "success"}
    except Exception as e:
        print(f"[ERROR] Exception in run_sample for {sample_name}: {e}")
        import traceback

        traceback.print_exc()
        return {"status": "error", "error": str(e)}


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
        "--repeatThresholdRelative",
        type=int,
        default=6,
        dest="repeat_threshold_set",
        help="The number of repeats needed for a read to be classified as telomeric and will be used as x per 100bp \
            input. Default is 6 repeats per 100 bp. Absolute repeats per read needed is minimum 4 after calculation. \
            Can be overwritten by --fixedRepeatThresholds.",
    )
    threshold_filtering_group.add_argument(
        "-frt",
        "--fixedRepeatThresholds",
        dest="fixed_repeat_thresholds",
        type=str,
        help=(
            "If read length distribution is known or you want to overwrite the x per 100bp calculation input via --repeatThreshold, "
            "use fixed thresholds for tumor and control samples. Input can be e.g. '6' (set for both samples) "
            "or '6,5' (for tumor and control respectively)."
        ),
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
        "--min-reads-per-barcode",
        type=int,
        dest="min_reads_per_barcode",
        default=10000,
        help="Minimum reads per barcode for single-cell summary only needed by telomerehunter2_sc.py not the bulk application (default: 10000).",
    )
    threshold_filtering_group.add_argument(
        "--barcode_tag",
        type=str,
        dest="barcode_tag",
        default="CB",
        help="SAM tag that contains the single-cell barcode (default: CB).",
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
    running_group.add_argument(
        "-fast",
        "--fast_mode",
        dest="fast_mode",
        action="store_true",
        help="Run in fast mode: extract and filter unmapped reads only.",
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
        default=["pdf", "html"],
        choices=["pdf", "png", "svg", "html", "all"],
        help=(
            "File format of output diagrams. Choose from: "
            "'pdf' (default), 'png', 'svg', 'html', or 'all' (generates all formats). "
            "The 'html' option creates interactive Plotly HTML files."
        ),
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

    return args


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
    # Determine available CPU cores
    available_cores = multiprocessing.cpu_count()
    # Determine maximum number of workers based on parallel execution setting and user input
    if args.run_parallel:
        # Divide available cores between tumor and control
        if args.cores is not None:
            max_workers = min(args.cores, 2)
            region_cores = max(1, args.cores // max_workers)
        else:
            max_workers = min(2, available_cores)
            region_cores = max(1, available_cores // max_workers)
    else:
        max_workers = 1
        region_cores = args.cores if args.cores is not None else available_cores
    # Warn if region_cores exceeds available_cores
    if region_cores > available_cores:
        print(
            f"Warning: Requested region cores ({region_cores}) exceed available cores ({available_cores}). Limiting to {available_cores}."
        )
        region_cores = available_cores
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
                region_cores=region_cores,
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
                region_cores=region_cores,
                args=args,
            )

            submitted_futures.append(control_future)

            if not args.run_parallel:
                control_future.result()

        # Wait for all futures to complete (only needed in parallel mode)
        if args.run_parallel:
            for future in submitted_futures:
                future.result()


def summary_log2(main_path, pid):
    summary_file = os.path.join(main_path, f"{pid}_summary.tsv")
    summary = pd.read_csv(summary_file, sep="\t", header=0)

    samples = summary["sample"]

    # Select columns for log2 calculation
    log2_cols = [
        col
        for col in summary.columns
        if col == "tel_content"
           or "_arbitrary_context_norm_by_intratel_reads" in col
           or "_singletons_norm_by_all_reads" in col
    ]

    # Indices for log2 columns
    log2_indices = [summary.columns.get_loc(col) for col in log2_cols]

    if "tumor" in samples.values and "control" in samples.values:
        tumor_row = (
            summary.loc[samples == "tumor"].iloc[:, log2_indices].fillna(0).values
        )
        control_row = (
            summary.loc[samples == "control"].iloc[:, log2_indices].fillna(0).values
        )

        log2_values = np.empty_like(tumor_row, dtype=float)
        log2_values.fill(np.nan)
        mask = (control_row != 0.0) & (tumor_row != 0.0)
        log2_values[mask] = np.log2(np.divide(tumor_row[mask], control_row[mask]))
        log2_values[np.isinf(log2_values)] = np.nan

        # Create a new row with [pid, "log2(tumor/control)"] and non-log empty columns
        new_row = pd.DataFrame(
            [[pid, "log2(tumor/control)"] + [np.nan] * (len(summary.columns) - 2)],
            columns=summary.columns,
        )
        # Add the calculated log2 values to the new row
        for idx, col in zip(log2_indices, log2_cols):
            new_row.iloc[0, idx] = log2_values[0, log2_cols.index(col)]

        summary = pd.concat([summary, new_row], ignore_index=True)

    summary = summary.convert_dtypes()
    summary.to_csv(summary_file, sep="\t", index=False)


def prepare_summary_file(main_path, pid):
    summary_file = os.path.join(main_path, f"{pid}_summary.tsv")
    summary = pd.read_csv(summary_file, sep="\t", header=0)

    # Rename 'read_length' to 'read_lengths'
    if "read_length" in summary.columns:
        summary = summary.rename(columns={"read_length": "read_lengths"})

    # Remove 'TRPM' column if present
    if "TRPM" in summary.columns:
        summary = summary.drop(columns=["TRPM"])

    # Move 'tel_content' to third column if present
    cols = list(summary.columns)
    if "tel_content" in cols:
        cols.remove("tel_content")
        cols.insert(2, "tel_content")
        summary = summary[cols]

    summary.to_csv(summary_file, sep="\t", index=False)


def run_plots(args, outdir, pid, repeat_thresholds_plot):
    print("------ making plots ------\n")
    if args.plotNone:
        print("All plots skipped as 'plotNone' is active.")
    else:
        # Create plots directory
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
        print(
            "Single-cell mode detected. Running full file analysis and barcode-specific insights."
        )

    # Fast mode branch
    if getattr(args, "fast_mode", False):
        run_fast_mode(args)
        return

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
        args.noFiltering,
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
    ) = get_repeat_threshold.get_read_lengths_and_repeat_thresholds(
        args, control_bam, tumor_bam
    )

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

    # Prepare summary file: rename, reorder, remove columns
    prepare_summary_file(args.outdir, args.pid)

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
