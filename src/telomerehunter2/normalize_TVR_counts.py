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

import numpy as np
import pandas as pd
import pysam


def count_reads_in_bam(bam_path):
    """
    Count the number of reads in a BAM file using pysam.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        return bam_file.count(until_eof=True)


def normalize_pattern_counts(sample_dir, pid, sample, pattern_table):
    """
    Normalize pattern counts for a given sample.

    Args:
        sample_dir (str): Path to the sample directory.
        pid (str): Process ID.
        pattern_table (pd.DataFrame): Original pattern table.

    Returns:
        pd.DataFrame: Normalized pattern table.
    """
    total_reads = (
        pd.read_table(os.path.join(sample_dir, f"{pid}_readcount.tsv"))["reads"]
        .astype(float)
        .sum()
    )
    pattern_table["Count_norm_by_all_reads"] = pattern_table["Count"] / total_reads

    bam_file_path = os.path.join(sample_dir, f"{pid}_filtered_intratelomeric.bam")
    total_intratelomeric = count_reads_in_bam(bam_file_path)
    pattern_table["Count_norm_by_intratel_reads"] = (
        pattern_table["Count"] / total_intratelomeric
    )

    summed_read_length_table = pd.read_table(
        os.path.join(sample_dir, "TVRs", f"{pid}_{sample}_summed_read_length.tsv"),
        sep="\t",
    )
    pattern_table["Count_per_100_bp_intratel_read"] = (
        pattern_table["Count"]
        * 100
        / summed_read_length_table.at[0, "summed_intratel_read_length"]
    )

    return pattern_table[
        [
            "Pattern",
            "Count_norm_by_all_reads",
            "Count_norm_by_intratel_reads",
            "Count_per_100_bp_intratel_read",
        ]
    ]


def merge_and_sort_tables(pattern_table_tumor, pattern_table_control):
    """
    Merge and sort pattern tables for tumor and control samples.

    Args:
        pattern_table_tumor (pd.DataFrame): Pattern table for the tumor sample.
        pattern_table_control (pd.DataFrame): Pattern table for the control sample.

    Returns:
        pd.DataFrame: Merged and sorted pattern table.
    """
    table_merged = pd.merge(
        pattern_table_tumor,
        pattern_table_control,
        on="Pattern",
        how="outer",
        suffixes=("_T", "_C"),
    )
    return table_merged.sort_values(
        by=["Count_norm_by_intratel_reads_C", "Count_norm_by_intratel_reads_T"],
        ascending=[False, False],
    )


def calculate_log2_ratios(table_merged):
    # Ensure both columns for normalized counts are valid (not NaN and not 0) for the division
    condition_norm = (
        (~table_merged["Count_norm_by_intratel_reads_T"].isna())
        & (~table_merged["Count_norm_by_intratel_reads_C"].isna())
        & (table_merged["Count_norm_by_intratel_reads_T"] != 0)
        & (table_merged["Count_norm_by_intratel_reads_C"] != 0)
    )

    # Calculate log2 ratio for normalized counts
    table_merged["log2_ratio_count_norm_by_intratel_reads"] = (
        np.nan
    )  # Initialize column with NaN

    valid_norm_values = table_merged.loc[
        condition_norm,
        ["Count_norm_by_intratel_reads_T", "Count_norm_by_intratel_reads_C"],
    ]

    table_merged.loc[condition_norm, "log2_ratio_count_norm_by_intratel_reads"] = (
        np.log2(
            valid_norm_values["Count_norm_by_intratel_reads_T"]
            / valid_norm_values["Count_norm_by_intratel_reads_C"]
        )
    )

    # Ensure both columns for per-100-bp counts are valid (not NaN and not 0) for the division
    condition_100_bp = (
        (~table_merged["Count_per_100_bp_intratel_read_T"].isna())
        & (~table_merged["Count_per_100_bp_intratel_read_C"].isna())
        & (table_merged["Count_per_100_bp_intratel_read_T"] != 0)
        & (table_merged["Count_per_100_bp_intratel_read_C"] != 0)
    )

    # Calculate log2 ratio for per-100-bp counts
    table_merged["log2_ratio_count_per_100_bp_intratel_read"] = (
        np.nan
    )  # Initialize column with NaN

    valid_100_bp_values = table_merged.loc[
        condition_100_bp,
        ["Count_per_100_bp_intratel_read_T", "Count_per_100_bp_intratel_read_C"],
    ]

    table_merged.loc[condition_100_bp, "log2_ratio_count_per_100_bp_intratel_read"] = (
        np.log2(
            valid_100_bp_values["Count_per_100_bp_intratel_read_T"]
            / valid_100_bp_values["Count_per_100_bp_intratel_read_C"]
        )
    )

    return table_merged


def update_summary_table(main_path, pid, TVRs_for_summary, table_merged):
    """
    Update the summary table with normalized counts.

    Args:
        main_path (str): Path to the main directory.
        pid (str): Process ID.
        TVRs_for_summary (str): Comma-separated list of TVRs for summary.
        table_merged (pd.DataFrame): Merged and sorted pattern table.

    Returns:
        None
    """
    summary_file = os.path.join(main_path, f"{pid}_summary.tsv")
    summary = pd.read_table(summary_file, sep="\t")

    TVR_colnames = [
        f"{t}_arbitrary_context_norm_by_intratel_reads" for t in TVRs_for_summary
    ]
    summary[TVR_colnames] = None

    for sample in ["tumor", "control"]:
        pattern_col = (
            f"Count_norm_by_intratel_reads_{'T' if sample == 'tumor' else 'C'}"
        )
        pattern_indexed = table_merged.set_index("Pattern")[pattern_col]

        # add TVRs not found but searched for
        pattern_indexed = pattern_indexed.reindex(TVRs_for_summary, fill_value=0)
        summary.loc[summary["sample"] == sample, TVR_colnames] = pattern_indexed.values

    summary.to_csv(summary_file, sep="\t", index=False)


def normalize_TVR_counts(main_path, pid, TVRs_for_summary):
    """
    Normalize TVR counts and update summary table.

    Args:
        main_path (str): Path to the main directory.
        pid (str): Process ID.
        TVRs_for_summary (str): Comma-separated list of TVRs for summary.

    Returns:
        None
    """
    TVRs_for_summary = TVRs_for_summary.split(",")

    pattern_table_list = {"tumor": None, "control": None}

    for sample in ["tumor", "control"]:
        sample_dir = os.path.join(main_path, f"{sample}_TelomerCnt_{pid}")
        pattern_table_path = os.path.join(
            sample_dir, "TVRs", f"{pid}_{sample}_TVRs.txt"
        )

        if not os.path.exists(pattern_table_path):
            # Add placeholder if no table is present
            pattern_table_list[sample] = pd.DataFrame(
                {
                    "Pattern": ["Placeholder"],
                    "Count_norm_by_all_reads": [np.nan],
                    "Count_norm_by_intratel_reads": [np.nan],
                    "Count_per_100_bp_intratel_read": [np.nan],
                }
            )
            continue

        original_pattern_table = pd.read_table(pattern_table_path)
        pattern_table_list[sample] = normalize_pattern_counts(
            sample_dir, pid, sample, original_pattern_table
        )

    table_merged = merge_and_sort_tables(
        pattern_table_list["tumor"], pattern_table_list["control"]
    )
    table_merged = calculate_log2_ratios(table_merged)

    # remove placeholder
    table_merged.drop(
        table_merged[table_merged["Pattern"] == "Placeholder"].index, inplace=True
    )

    # save files
    table_merged.to_csv(
        os.path.join(main_path, f"{pid}_normalized_TVR_counts.tsv"),
        sep="\t",
        index=False,
    )
    update_summary_table(main_path, pid, TVRs_for_summary, table_merged)
