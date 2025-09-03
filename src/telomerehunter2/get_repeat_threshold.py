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

import numpy as np
import pysam


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

        # Enforce minimum threshold of 4
        if repeat_threshold_per_100_bp < 4:
            print(
                "Warning: repeat_threshold_per_100_bp was set below 4. Setting to minimum value of 4."
            )
            repeat_threshold_per_100_bp = 4

        read_lengths = list(map(int, sorted_read_length_str.split(",")))

        repeat_thresholds = [
            int(round(float(i) * repeat_threshold_per_100_bp / 100))
            for i in read_lengths
        ]

        # Threshold can never be less than 4
        repeat_thresholds = [max(4, t) for t in repeat_thresholds]

        unique_repeat_thresholds = sorted(set(repeat_thresholds))

        if len(unique_repeat_thresholds) == 1:
            repeat_threshold = unique_repeat_thresholds[0]
            print(f"Single unique repeat threshold: {repeat_threshold}")

        elif len(unique_repeat_thresholds) > 1:
            print(f"Multiple Unique Read Lengths: {read_lengths}")

            weights = [read_length_counts.get(length, 1) for length in read_lengths]

            repeat_threshold = int(
                round(np.average(repeat_thresholds, weights=weights))
            )

            # Also ensure at least 4 here
            repeat_threshold = max(4, repeat_threshold)

            print(f"Calculating the Weighted Repeat Threshold: {repeat_threshold}")

        else:
            print("Error: Unable to calculate repeat threshold.")
            repeat_threshold = None

        return repeat_threshold, (
            str(repeat_threshold) if repeat_threshold is not None else None
        )

    except Exception as e:
        print(f"Unexpected error in repeat threshold calculation: {e}")
        return None, None
