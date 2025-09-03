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

from telomerehunter2.utils import get_reverse_complement, measure_time

###############################################################################################
### get the distribution of telomeric repeats per intratelomeric read in the input BAM file ###
###############################################################################################


@measure_time
def repeat_frequency_intratelomeric(input_path, out_dir, pid, repeats):
    ################################################
    ### get patterns and make regular expression ###
    ################################################

    # Create a regular expression pattern for forward and reverse repeats
    patterns_regex_forward = "|".join(repeats)
    patterns_regex_reverse = "|".join(
        get_reverse_complement(repeat) for repeat in repeats
    )

    #########################
    ### open file handles ###
    #########################

    # open input bam_file for reading
    bamfile = pysam.AlignmentFile(
        os.path.join(input_path, f"{pid}_filtered_intratelomeric.bam"), "rb"
    )

    ##################################
    ### initialize frequency table ###
    ##################################

    frequency_table = {repeats: 0 for repeats in range(0, 2 + 1)}

    ######################################
    ### loop through filtered BAM file ###
    ######################################

    for read in bamfile.fetch(until_eof=True):
        sequence = read.query_sequence
        number_repeats_forward = len(re.findall(patterns_regex_forward, sequence))
        number_repeats_reverse = len(re.findall(patterns_regex_reverse, sequence))
        number_repeats = max(number_repeats_forward, number_repeats_reverse)

        # Increment count for the number of repeats or set to 1 if the key is not present
        try:
            frequency_table[number_repeats] += 1
        except KeyError:
            frequency_table[number_repeats] = 1

    ##################################
    ### write frequency table file ###
    ##################################

    frequency_table_filepath = os.path.join(
        out_dir, f"{pid}_repeat_frequency_per_intratelomeric_read.tsv"
    )
    with open(frequency_table_filepath, "w") as frequency_table_file:
        # Header
        frequency_table_file.write("number_repeats\tcount\n")

        # Write line for each frequency
        for frequency, count in sorted(frequency_table.items()):
            frequency_table_file.write(f"{frequency}\t{count}\n")

    ##########################
    ### close file handles ###
    ##########################

    bamfile.close()
