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

import pysam

#########################################################################
### gives the sum of all intratelomeric read lengths (base pairs)     ###
### secondary and supplementary alignments are skipped                ###
#########################################################################


def summed_intratelomeric_read_length(main_path, pid, sample):
    """
    Calculate the sum of intratelomeric read lengths and save the result in a TSV file.

    :param main_path: The main directory path.
    :param pid: Patient ID.
    :param sample: Sample ID.
    """

    # Output file path
    outfile_path = os.path.join(
        main_path, "TVRs", f"{pid}_{sample}_summed_read_length.tsv"
    )

    # Open output file for writing
    with open(outfile_path, "w") as outfile:
        # Write header
        outfile.write("PID\tsample\tsummed_intratel_read_length\n")

        # Input BAM file path
        bam_file_path = os.path.join(main_path, f"{pid}_filtered_intratelomeric.bam")

        # Open input BAM file for reading
        with pysam.AlignmentFile(bam_file_path, "rb") as bamfile:
            summed_read_length = 0

            # Iterate through reads
            for read in bamfile.fetch(until_eof=True):
                # Accumulate read length
                summed_read_length += len(read.query_sequence)

            # Write result to the output file
            outfile.write(f"{pid}\t{sample}\t{summed_read_length}\n")
