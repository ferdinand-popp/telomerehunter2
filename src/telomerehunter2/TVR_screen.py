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
import re

import pysam

from telomerehunter2.utils import (assure_dir_exists, get_reverse_complement)

##################################################################################################################################################################
### script loops through filtered BAM file containing intratelomeric reads and searches for patterns of the type XXXGGG (and reverse complement)               ###
### only counts patterns if the base qualities are all greater than 20                                                                                         ###
### pattern counts are written to output tables (with the frequency of the pattern and the average base quality at each position)                              ###
##################################################################################################################################################################


def find_prioritized_telomeric_indices(sequence, telomere_pattern="(.{3})(GGG)"):
    indices = [m.start() for m in re.finditer(rf"(?={telomere_pattern})", sequence)]
    valid_indices = []
    i = 0

    while i < len(indices):
        current_index = indices[i]

        if i + 1 < len(indices):
            next_index = indices[i + 1]

            if next_index - current_index == 1:
                # Handle overlapping sequences
                current_hexamer = sequence[current_index : current_index + 6]
                next_hexamer = sequence[next_index : next_index + 6]

                # Determine neighbors
                current_neighbors = [
                    idx
                    for idx in indices
                    if idx != current_index and abs(idx - current_index) == 6
                ]
                next_neighbors = [
                    idx
                    for idx in indices
                    if idx != next_index and abs(idx - next_index) == 6
                ]

                if "TTAGGG" in current_hexamer and "TTAGGG" not in next_hexamer:
                    valid_indices.append(current_index)
                    i += 2  # Skip the next index since we are keeping the current one
                elif "TTAGGG" not in current_hexamer and "TTAGGG" in next_hexamer:
                    valid_indices.append(next_index)
                    i += 2  # Skip the next index since we are keeping it
                else:
                    # Both are TTAGGG or neither are TTAGGG; choose based on neighbors
                    if len(current_neighbors) > len(next_neighbors):
                        valid_indices.append(current_index)
                    else:
                        valid_indices.append(next_index)
                    i += 2  # Skip the next index after deciding
            else:
                valid_indices.append(current_index)
                i += 1
        else:
            valid_indices.append(current_index)
            i += 1

    return valid_indices


def screen_bam_file(
    bam_file, repeat_threshold_set, telomere_pattern="(.{3})(GGG)", qual_threshold=20
):
    """Screen BAM file for patterns and update counts and qualities."""

    patterns = {}  # dict for counting of patterns
    qualities = (
        {}
    )  # dict: patterns are keys, for each pattern 3 lists containing qualities at position 1, 2 and 3

    bamfile = pysam.AlignmentFile(bam_file, "rb")

    for read in bamfile.fetch(until_eof=True):
        seq = read.query_sequence
        # Use the prioritized telomeric indices function
        indices = find_prioritized_telomeric_indices(seq, telomere_pattern)

        # Get reverse complement if the pattern was not found often enough
        if len(indices) < repeat_threshold_set:
            seq = get_reverse_complement(seq)
            qual = read.query_qualities[::-1]
            indices = find_prioritized_telomeric_indices(seq, telomere_pattern)
        else:
            qual = read.query_qualities

        for i in indices:
            pattern_segment = seq[i : i + 6]

            if "N" in pattern_segment:  # skip if pattern contains an "N"
                continue

            # Extract the quality scores for the pattern segment
            q_scores = qual[i : i + 6]

            # Skip if one of the positions has a base quality lower than the threshold
            if any(q < qual_threshold for q in q_scores):
                continue

            # Record or initialize pattern counts and qualities
            if pattern_segment in patterns:
                patterns[pattern_segment] += 1
                for j in range(6):
                    qualities[pattern_segment][j].append(q_scores[j])
            else:
                patterns[pattern_segment] = 1
                qualities[pattern_segment] = [[q_scores[j]] for j in range(6)]

    return patterns, qualities


def show_patterns(pattern_dict, qualities_dict):
    """Sorts and displays patterns and qualites"""
    keys = [
        key
        for _, key in sorted(((v, k) for k, v in pattern_dict.items()), reverse=True)
    ]

    output = (
        "\t".join(
            [
                "Pattern",
                "Count",
                "Frequency_in_Percent",
                "Avg_Qual_pos1",
                "Avg_Qual_pos2",
                "Avg_Qual_pos3",
                "Avg_Qual_pos4",
                "Avg_Qual_pos5",
                "Avg_Qual_pos6",
            ]
        )
        + "\n"
    )

    counts = sum(pattern_dict.values())

    for key in keys:
        (qs1, qs2, qs3, qs4, qs5, qs6) = qualities_dict[
            key
        ]  # lists of qualites for position 1, 2 and 3
        q1_mean = sum(qs1) / len(qs1)
        q2_mean = sum(qs2) / len(qs2)
        q3_mean = sum(qs3) / len(qs3)
        q4_mean = sum(qs4) / len(qs4)
        q5_mean = sum(qs5) / len(qs5)
        q6_mean = sum(qs6) / len(qs6)
        output += (
            "\t".join(
                [
                    key,
                    str(pattern_dict[key]),
                    str(pattern_dict[key] * 100.0 / counts),
                    str(q1_mean),
                    str(q2_mean),
                    str(q3_mean),
                    str(q4_mean),
                    str(q5_mean),
                    str(q6_mean),
                ]
            )
            + "\n"
        )

    return output


def tvr_screen(
    main_path,
    pid,
    sample,
    repeat_threshold_set,
    min_base_quality=20,
    telomere_pattern="(.{3})(GGG)",
):
    """Main function for TVR screening."""

    path_bam_intratelomeric = os.path.join(
        main_path, f"{pid}_filtered_intratelomeric.bam"
    )

    pattern_dict, qualities_dict = screen_bam_file(
        path_bam_intratelomeric,
        repeat_threshold_set,
        telomere_pattern=telomere_pattern,
        qual_threshold=min_base_quality,
    )

    outfile_path = os.path.join(main_path, "TVRs", f"{pid}_{sample}_TVRs.txt")
    assure_dir_exists(outfile_path)

    with open(outfile_path, "w") as outfile:
        outfile.write(show_patterns(pattern_dict, qualities_dict))
