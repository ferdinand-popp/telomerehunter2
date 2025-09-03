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


def process_pattern(
    main_path, pid, standard_repeat, bp_context, TVRs_for_context, summary
):
    df_list_all = {}
    multimer_repeats_count = bp_context // len(standard_repeat)

    for pattern in TVRs_for_context:
        # Initialize DataFrame with proper dtypes
        top_contexts = pd.DataFrame()

        for sample in summary["sample"].unique():
            TVR_context_dir = os.path.join(
                main_path, f"{sample}_TelomerCnt_{pid}", "TVR_context"
            )
            neighborhood_table_file = os.path.join(
                TVR_context_dir,
                f"{pid}_{sample}_{pattern}_{bp_context}bp_{bp_context}bp_neighborhood.tsv",
            )

            try:
                neighborhood_table = pd.read_table(neighborhood_table_file, sep="\t")
            except FileNotFoundError:
                print(f"Didnt find {neighborhood_table_file}")
                continue

            row_name = f"{pid}_{sample}_{bp_context}"

            # Basic fields
            top_contexts.loc[row_name, "PID"] = pid
            top_contexts.loc[row_name, "Sample"] = sample
            top_contexts.loc[row_name, "Context_bp"] = bp_context
            top_contexts.loc[row_name, "pattern"] = pattern

            # Handle empty neighborhood table
            if neighborhood_table.empty:
                print(f"{neighborhood_table_file} was empty")
                top_contexts.loc[row_name, "Bases"] = ""
                top_contexts.loc[row_name, "Count"] = 0
                top_contexts.loc[row_name, "Percent"] = 0.0
                top_contexts.loc[row_name, "Bases_t_type"] = ""
                top_contexts.loc[row_name, "Count_t_type"] = 0
                top_contexts.loc[row_name, "Percent_t_type"] = 0.0
                continue

            # Safe value assignments with type conversion
            top_contexts.loc[row_name, "Bases"] = str(
                neighborhood_table.loc[0, "Bases"]
            )
            top_contexts.loc[row_name, "Count"] = int(
                neighborhood_table.loc[0, "Count"]
            )
            top_contexts.loc[row_name, "Percent"] = float(
                neighborhood_table.loc[0, "Percent"]
            )

            # find singletons
            pattern_t_type_context = f"{standard_repeat * multimer_repeats_count}-{pattern}-{standard_repeat * multimer_repeats_count}"

            if pattern_t_type_context in neighborhood_table["Bases"].values:
                idx = neighborhood_table.index[
                    neighborhood_table["Bases"] == pattern_t_type_context
                ][0]
                top_contexts.loc[row_name, "Bases_t_type"] = str(
                    neighborhood_table.loc[idx, "Bases"]
                )
                top_contexts.loc[row_name, "Count_t_type"] = int(
                    neighborhood_table.loc[idx, "Count"]
                )
                top_contexts.loc[row_name, "Percent_t_type"] = float(
                    neighborhood_table.loc[idx, "Percent"]
                )
            else:
                top_contexts.loc[row_name, "Bases_t_type"] = ""
                top_contexts.loc[row_name, "Count_t_type"] = 0
                top_contexts.loc[row_name, "Percent_t_type"] = 0.0

            # Handle telomere content ratio
        if (
            "tumor" in summary["sample"].values
            and "control" in summary["sample"].values
        ):
            tel_content_tumor = summary.loc[
                summary["sample"] == "tumor", "tel_content"
            ].values[0]
            tel_content_control = summary.loc[
                summary["sample"] == "control", "tel_content"
            ].values[0]

            if tel_content_tumor == 0 or tel_content_control == 0:
                top_contexts["tel_content_log2_ratio"] = 0
            else:
                top_contexts["tel_content_log2_ratio"] = np.log2(
                    tel_content_tumor / tel_content_control
                )
        else:
            top_contexts["tel_content_log2_ratio"] = 0

        df_list_all[pattern] = top_contexts

    # Check if df_list_all is not empty before attempting to concatenate
    if df_list_all:
        pattern_contexts = pd.concat(df_list_all.values(), keys=df_list_all.keys())
    else:
        pattern_contexts = pd.DataFrame()  # Return an empty DataFrame
        print("Pattern context DF empty.")

    return pattern_contexts


def process_singleton_table(pattern_contexts, summary):
    singleton_table = pattern_contexts[
        ["PID", "pattern", "tel_content_log2_ratio", "Sample", "Count_t_type"]
    ]

    # Pivot the table using pivot_table to get Count_t_type per sample
    singleton_table = pd.pivot_table(
        singleton_table,
        values="Count_t_type",
        index=["PID", "pattern", "tel_content_log2_ratio"],
        columns="Sample",
        aggfunc="sum",
    )
    # Reset the index to make it flat
    singleton_table.reset_index(inplace=True)
    # Handle NaN values
    singleton_table.fillna(0, inplace=True)

    # Rename columns if 'control' and 'tumor' are present
    if "control" in singleton_table.columns:
        singleton_table.rename(
            columns={"control": "singleton_count_control"}, inplace=True
        )
        singleton_table["singleton_count_control"] = singleton_table[
            "singleton_count_control"
        ].astype(int)
    else:
        singleton_table["singleton_count_control"] = np.nan

    if "tumor" in singleton_table.columns:
        singleton_table.rename(columns={"tumor": "singleton_count_tumor"}, inplace=True)
        singleton_table["singleton_count_tumor"] = singleton_table[
            "singleton_count_tumor"
        ].astype(int)
    else:
        singleton_table["singleton_count_tumor"] = np.nan

    # Get log2 ratio of counts
    singleton_table["singleton_count_log2_ratio"] = np.log2(
        singleton_table["singleton_count_tumor"]
        / singleton_table["singleton_count_control"]
    )

    # Get total reads each for normalization
    total_reads_control = (
        summary.loc[summary["sample"] == "control", "total_reads"].values[0]
        if "control" in summary["sample"].values
        else 0
    )
    singleton_table["singleton_count_control_norm"] = (
        singleton_table["singleton_count_control"] / total_reads_control
    )

    total_reads_tumor = (
        summary.loc[summary["sample"] == "tumor", "total_reads"].values[0]
        if "tumor" in summary["sample"].values
        else 0
    )
    singleton_table["singleton_count_tumor_norm"] = (
        singleton_table["singleton_count_tumor"] / total_reads_tumor
    )

    singleton_table["singleton_count_log2_ratio_norm"] = np.log2(
        singleton_table["singleton_count_tumor_norm"]
        / singleton_table["singleton_count_control_norm"]
    )

    singleton_table["distance_to_expected_singleton_log2_ratio"] = (
        singleton_table["singleton_count_log2_ratio_norm"]
        - singleton_table["tel_content_log2_ratio"]
    )

    # move tel_ value to distance to make it better understandable
    tel_content_log2_ratio = singleton_table.pop("tel_content_log2_ratio")
    # Insert it at the second-to-last position
    singleton_table.insert(
        len(singleton_table.columns) - 1,
        "tel_content_log2_ratio",
        tel_content_log2_ratio,
    )

    return singleton_table


def add_tvr_normalized_counts(summary, TVRs_for_context, singleton_table):
    # Add column names
    singleton_colnames = [
        f"{pattern}_singletons_norm_by_all_reads" for pattern in TVRs_for_context
    ]
    for col in singleton_colnames:
        summary[col] = np.nan

    # Check if singleton_table is empty
    if singleton_table.empty:
        print(
            "Warning: singleton_table is empty. No TVR normalized counts will be added to the summary."
        )
        return summary, singleton_table

    # Set singleton table index to patterns
    singleton_table.index = singleton_table["pattern"]

    for sample_type in ["tumor", "control"]:
        if (singleton_table[f"singleton_count_{sample_type}"] != 0).any():
            mask = summary["sample"] == sample_type

            # Iterating through the present TVRs since only the ones found are in singleton table
            for pattern in TVRs_for_context:
                if pattern in singleton_table.index:
                    summary.loc[mask, f"{pattern}_singletons_norm_by_all_reads"] = (
                        singleton_table.loc[
                            pattern, f"singleton_count_{sample_type}_norm"
                        ]
                    )
        else:
            print(
                f"Warning: singleton_table {sample_type} has no non-zero counts. No TVR normalized counts will be added."
            )

    return summary, singleton_table


def save_results(summary, pattern_contexts, singleton_table, main_path, pid):
    try:
        summary.to_csv(
            os.path.join(main_path, f"{pid}_summary.tsv"), sep="\t", index=False
        )
    except Exception as e:
        print(f"Error saving extended summary table: {e}")

    if pattern_contexts.empty:
        print("Warning: pattern_contexts is empty. Saving an empty file.")

    patterns_top_contexts = pattern_contexts[
        ["PID", "Sample", "pattern", "Context_bp", "Bases", "Count", "Percent"]
    ]

    patterns_top_contexts.to_csv(
        os.path.join(main_path, f"{pid}_TVR_top_contexts.tsv"), sep="\t", index=False
    )

    if singleton_table.empty:
        print("Warning: singleton_table is empty. Saving an empty file.")
    singleton_table.to_csv(
        os.path.join(main_path, f"{pid}_singletons.tsv"), sep="\t", index=False
    )


def process_and_save_results(
    main_path, pid, summary, TVRs_for_context, pattern_contexts, singleton_table
):
    # add TVR counts to summary
    summary, singleton_table = add_tvr_normalized_counts(
        summary, TVRs_for_context, singleton_table
    )

    save_results(summary, pattern_contexts, singleton_table, main_path, pid)


def tvr_context_singletons_tables(
    main_path, pid, standard_repeat, bp_context, TVRs_for_context
):
    TVRs_for_context = TVRs_for_context.split(",")

    summary_file = os.path.join(main_path, f"{pid}_summary.tsv")
    summary = pd.read_table(summary_file, sep="\t")

    pattern_contexts = process_pattern(
        main_path, pid, standard_repeat, bp_context, TVRs_for_context, summary
    )

    if not pattern_contexts.empty and pattern_contexts["Count"].sum() > 0:
        singleton_table = process_singleton_table(pattern_contexts, summary)
    else:
        singleton_table = pd.DataFrame()

    process_and_save_results(
        main_path, pid, summary, TVRs_for_context, pattern_contexts, singleton_table
    )
