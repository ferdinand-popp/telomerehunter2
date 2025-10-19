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

from collections import Counter

import os
import numpy as np
import pysam
import plotly.graph_objs as go


def get_read_lengths(bam_file_path, reads_to_parse=1000):
    """
    Get read lengths and their counts from the first N non-supplementary or secondary alignments in a BAM file.

    :param bam_file_path: Path to the input BAM file.
    :param reads_to_parse: Number of reads to parse, default is 1000.
    :return: Tuple containing (Comma-separated string of unique read lengths, read length counts dict)
    """
    # Open input bam_file for reading
    mode = "rb" if bam_file_path.endswith(".bam") else "rc"
    bamfile = pysam.AlignmentFile(bam_file_path, mode)

    # print unique read lengths of the first N non-supplementary or secondary alignments
    cntr = 0
    read_lengths = []

    for read in bamfile.fetch(until_eof=True):
        # Skip secondary alignments
        if read.is_secondary:
            continue

        # Skip supplementary alignments
        if read.flag >= 2048:
            continue

        read_lengths.append(len(read.query_sequence))

        cntr += 1
        if cntr == reads_to_parse:
            break

    # Count occurrences of each read length
    read_length_counts = dict(Counter(read_lengths))
    read_lengths = sorted(list(set(read_lengths)))
    read_lengths_str = ",".join(str(i) for i in read_lengths)

    return read_lengths_str, read_length_counts


def get_repeat_threshold(
        sorted_read_length_str, read_length_counts, repeat_threshold_per_100_bp
):
    """
    Calculate the repeat threshold based on read lengths.
    The threshold can never be less than 4.

    :param sorted_read_length_str: Comma-separated string of read lengths
    :param read_length_counts: Dictionary of read length counts
    :param repeat_threshold_per_100_bp: Threshold per 100bp (int, minimum 4)
    :return: Tuple (repeat_threshold, repeat_threshold_str)
    """
    try:
        # Input validation
        if repeat_threshold_per_100_bp is None:
            raise ValueError("repeat_threshold_per_100_bp must not be None.")
        if not isinstance(repeat_threshold_per_100_bp, int):
            raise ValueError("repeat_threshold_per_100_bp must be an integer.")

        read_lengths = list(map(int, sorted_read_length_str.split(",")))

        repeat_thresholds = [
            int(round(float(i) * repeat_threshold_per_100_bp / 100))
            for i in read_lengths
        ]

        unique_repeat_thresholds = sorted(set(repeat_thresholds))

        if len(unique_repeat_thresholds) == 1:
            repeat_threshold = unique_repeat_thresholds[0]
            print(f"Single repeat threshold: {repeat_threshold}")

        elif len(unique_repeat_thresholds) > 1:
            weights = [read_length_counts.get(length, 1) for length in read_lengths]

            repeat_threshold = int(
                round(np.average(repeat_thresholds, weights=weights))
            )

            print(f"Calculated the Weighted Repeat Threshold: {repeat_threshold} for multiple repeat lengths: {read_lengths}")


        else:
            print("Error: Unable to calculate repeat threshold.")
            repeat_threshold = None

        if repeat_threshold < 4:
            print(
                f"!! Repeat threshold {repeat_threshold} is less than minimum of 4. Setting to 4."
            )
            repeat_threshold = 4

        return repeat_threshold, (
            str(repeat_threshold) if repeat_threshold is not None else None
        )

    except Exception as e:
        print(f"Unexpected error in repeat threshold calculation: {e}")
        return None, None


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

    def _write_readlength_diagnostics(
            sample_name,
            read_length_counts,
            outdir,
            threshold_selected=None,
            control_read_length_counts=None,
            threshold_selected_control=None,
    ):
        # Prepare data for the histogram plot
        data = []
        if read_length_counts and len(read_length_counts) > 1:
            sample_lengths = np.repeat(
                list(read_length_counts.keys()),
                list(read_length_counts.values())
            )
            data.append(go.Histogram(
                x=sample_lengths,
                name=sample_name,
                opacity=0.65,
                marker_color='blue',
            ))

        if control_read_length_counts and len(control_read_length_counts) > 0:
            control_lengths = np.repeat(
                list(control_read_length_counts.keys()),
                list(control_read_length_counts.values())
            )
            data.append(go.Histogram(
                x=control_lengths,
                name="control",
                opacity=0.65,
                marker_color='red',
            ))

        if len(data) == 0:
            print(f"No multiple read length counts for {sample_name} so no diagnostic plot.")
            return

        # Create the figure
        fig = go.Figure(data=data)

        # Add average lines
        if read_length_counts and len(read_length_counts) > 0:
            avg_sample = np.average(
                list(read_length_counts.keys()),
                weights=list(read_length_counts.values())
            )
            fig.add_shape(
                type="line",
                x0=avg_sample,
                x1=avg_sample,
                y0=0,
                y1=1,
                xref="x",
                yref="paper",
                line=dict(color="blue", dash="dash"),
                name=f'{sample_name} Average'
            )

        if control_read_length_counts and len(control_read_length_counts) > 0:
            avg_control = np.average(
                list(control_read_length_counts.keys()),
                weights=list(control_read_length_counts.values())
            )
            fig.add_shape(
                type="line",
                x0=avg_control,
                x1=avg_control,
                y0=0,
                y1=1,
                xref="x",
                yref="paper",
                line=dict(color="red", dash="dash"),
                name='control Average'
            )

        # Update layout
        fig.update_layout(
            barmode='overlay',
            title=f"Read Length Distribution for {sample_name}{' & control' if control_read_length_counts else ''}",
            xaxis_title="Read Length",
            yaxis_title="Count",
            legend_title="Legend"
        )

        # Compose subtitle based on threshold(s) selected
        subtitle = None
        if threshold_selected is not None and threshold_selected_control is not None:
            subtitle = (
                f"Thresholds selected: {sample_name}: {threshold_selected} repeats/read, "
                f"Control: {threshold_selected_control} repeats/read"
            )
        elif threshold_selected is not None:
            subtitle = f"Threshold selected: {sample_name}: {threshold_selected} repeats/read"
        elif threshold_selected_control is not None:
            subtitle = f"Threshold selected: control: {threshold_selected_control} repeats/read"

        if subtitle:
            fig.update_layout(
                title={
                    'text': fig.layout.title.text + f"<br><sup>{subtitle}</sup>",
                    'x': 0.5
                }
            )

        # Save the plot
        os.makedirs(os.path.join(outdir, "html_reports"), exist_ok=True)
        html_path = os.path.join(outdir, "html_reports", f"{sample_name}_read_length_histogram.html")
        fig.write_html(html_path)
        del fig
        print(f"Histogram plot written to {html_path}")

    # Get read lengths for tumor and control files
    if args.tumor_flag:
        # Get read lengths and calculate thresholds for tumor
        read_lengths_str_tumor, tumor_read_length_counts = (
            get_read_lengths(tumor_bam)
        )
    if args.control_flag:
        read_lengths_str_control, control_read_length_counts = (
            get_read_lengths(control_bam)
        )

    # use either fixed repeats or calculate thresholds per 100 bp
    if args.fixed_repeat_thresholds:
        # Override with experienced user input
        thresholds = list(map(int, args.fixed_repeat_thresholds.split(',')))
        if len(thresholds) == 1:
            repeat_thresholds_tumor = repeat_thresholds_control = thresholds[0]
        elif len(thresholds) == 2:  # tumor and control different
            repeat_thresholds_tumor = thresholds[0]
            repeat_thresholds_control = thresholds[1]
        repeat_thresholds_str_tumor = str(repeat_thresholds_tumor)
        repeat_thresholds_str_control = str(repeat_thresholds_control)
        repeat_thresholds_plot = (
            repeat_thresholds_tumor
            if repeat_thresholds_tumor == repeat_thresholds_control
            else "n"
        )
        print(
            f"Using user-specified fixed repeat thresholds: Tumor={repeat_thresholds_tumor}, Control={repeat_thresholds_control}")
    else:
        # Calculate tumor thresholds if needed per 100 bp
        if args.tumor_flag:
            print("Calculating repeat threshold for the tumor sample: ")
            repeat_thresholds_tumor, repeat_thresholds_str_tumor = (
                get_repeat_threshold(
                    read_lengths_str_tumor,
                    tumor_read_length_counts,
                    args.repeat_threshold_set,
                )
            )

        # Calculate control thresholds if needed
        if args.control_flag:
            print("Calculating repeat threshold for the control sample: ")
            repeat_thresholds_control, repeat_thresholds_str_control = (
                get_repeat_threshold(
                    read_lengths_str_control,
                    control_read_length_counts,
                    args.repeat_threshold_set,
                )
            )

    # Determine which threshold to use for plotting and plot distribution and choice threshold
    if args.tumor_flag and args.control_flag:
        repeat_thresholds_plot = (
            repeat_thresholds_tumor
            if repeat_thresholds_tumor == repeat_thresholds_control
            else "n"
        )
        _write_readlength_diagnostics(
            'tumor',
            tumor_read_length_counts,
            args.outdir,
            threshold_selected=repeat_thresholds_tumor,
            control_read_length_counts=control_read_length_counts,
            threshold_selected_control=repeat_thresholds_control,
        )
    elif args.tumor_flag:
        repeat_thresholds_plot = repeat_thresholds_tumor
        _write_readlength_diagnostics(
            'tumor',
            tumor_read_length_counts,
            args.outdir,
            threshold_selected=repeat_thresholds_tumor,
        )
    elif args.control_flag:
        repeat_thresholds_plot = repeat_thresholds_control
        _write_readlength_diagnostics(
            'control',
            control_read_length_counts,
            args.outdir,
            threshold_selected=repeat_thresholds_control,
        )

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
