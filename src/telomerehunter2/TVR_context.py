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

import pandas as pd
import pysam

from telomerehunter2.utils import get_reverse_complement

#####################################################################################
###
### get the context of telomere variant repeats (TVRs)
###
#####################################################################################


def dictionary_to_table(dictionary, outfile_path, cutoff):
    try:
        # write sorted dictionary as a table to file
        df = pd.DataFrame.from_dict(dictionary, orient="index").reset_index()
        df.columns = ["Bases", "Count"]
        df = df.sort_values(by="Count", ascending=False)

        # add per all instances
        df["Percent"] = (df["Count"] / df["Count"].sum()) * 100

        # remove ones below cutoff
        df = df[df["Count"] >= cutoff]

        df = df.convert_dtypes()
        df.to_csv(outfile_path, sep="\t", index=False)

    except (ValueError, KeyError):
        # Handle the case where the dictionary is empty
        df = pd.DataFrame(columns=["Bases", "Count", "Percent"])
        df.to_csv(outfile_path, sep="\t", index=False)


def tvr_context(
    main_path,
    pid,
    sample,
    pattern,
    min_base_quality=20,
    context_before=18,
    context_after=18,
    telomere_pattern="GGG",
    repeat_threshold_set=6,
    cutoff=0,
    tel_file="filtered_intratelomeric",
):
    """
    # main_path: Path to TelomereHunter results
    # pattern: Pattern for which to get context
    # min_base_quality: Minimum base quality required for pattern
    # context_before: Number of bases before start of pattern
    # context_after: Number of bases after end of pattern
    # telomere_pattern: Pattern with which to identify G-rich reads
    # repeat_threshold_set: int number of telomeric reads on read threshold
    # cutoff: Count cutoff for displaying neighborhood in output table
    # tel_file: the telomere file in which to search for TVRs (filtered, filtered_intratelomeric, ...)
    """
    outdir = os.path.join(main_path, "TVR_context")

    os.makedirs(outdir, exist_ok=True)

    # Get BAM reads
    bam_file_path = os.path.join(main_path, f"{pid}_{tel_file}.bam")
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")

    # directories for counting of patterns in neighborhood
    neighborhood_before = {}
    neighborhood_after = {}
    neighborhood = {}

    for read in bamfile.fetch(until_eof=True):
        seq = read.query_sequence
        qual = read.query_qualities

        # find GGG counted once for reverse complement decision
        indices_telomeric = [m.start() for m in re.finditer(telomere_pattern, seq)]

        # get reverse complement of sequence if the telomeric pattern was not found often enough
        if len(indices_telomeric) < repeat_threshold_set:
            seq = get_reverse_complement(seq)
            qual = read.query_qualities[::-1]

        indices_pattern = [m.start() for m in re.finditer(rf"(?={pattern})", seq)]

        for i in indices_pattern:
            bases_before = ""
            bases_after = ""

            # get base qualities of pattern
            base_qualities = [qual[i + j] for j in range(6)]
            if any(q < min_base_quality for q in base_qualities):
                continue

            # todo check if keyerror needs to be caught
            if i - context_before >= 0:
                bases_before = seq[i - context_before : i]
                neighborhood_before[bases_before] = (
                    neighborhood_before.get(bases_before, 0) + 1
                )

            if i + len(pattern) + context_after <= len(seq):
                bases_after = seq[i + len(pattern) : i + len(pattern) + context_after]
                neighborhood_after[bases_after] = (
                    neighborhood_after.get(bases_after, 0) + 1
                )

            if bases_before and bases_after:
                bases_all = f"{bases_before}-{pattern}-{bases_after}"
                neighborhood[bases_all] = neighborhood.get(bases_all, 0) + 1

    dictionary_to_table(
        neighborhood_before,
        os.path.join(
            outdir,
            f"{pid}_{sample}_{pattern}_{context_before}bp_neighborhood_before.tsv",
        ),
        cutoff=cutoff,
    )
    dictionary_to_table(
        neighborhood_after,
        os.path.join(
            outdir, f"{pid}_{sample}_{pattern}_{context_after}bp_neighborhood_after.tsv"
        ),
        cutoff=cutoff,
    )
    dictionary_to_table(
        neighborhood,
        os.path.join(
            outdir,
            f"{pid}_{sample}_{pattern}_{context_before}bp_{context_after}bp_neighborhood.tsv",
        ),
        cutoff=cutoff,
    )
