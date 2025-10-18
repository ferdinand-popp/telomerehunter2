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

import csv
import os

import numpy
import pandas as pd
import pysam


##################################################################
### get the gc content distribution of the reads in a bam file ###
##################################################################


def get_gc_content_distribution(bam_file, out_dir, pid, sample, remove_duplicates):
    """
    Get the GC content distribution of the reads in a BAM file.
    """

    # open input bam_file for reading
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    # make GC content list
    gc_content_list = {
        gc_content: 0 for gc_content in range(0, 100 + 1)
    }  # TODO use samtools stat gcf method ?

    for read in bamfile.fetch(until_eof=True):

        if (
                read.is_secondary
        ):  # skips all secondary alignments (only needed for RNA analysis!)
            continue

        if (
                remove_duplicates and read.is_duplicate
        ):  # if duplicate flag is set: skip all reads that are flagged as optical or PCR duplicate
            continue

        if read.flag >= 2048:  # skip supplementary alignments
            continue

        sequence = read.query_sequence
        read_length = len(read.query_sequence)

        n_count = sequence.count("N")

        if float(n_count) / float(read_length) <= 0.2:
            gc_content = int(
                round(
                    float(sequence.count("C") + sequence.count("G"))
                    / float(read_length - n_count)
                    * 100
                )
            )
            gc_content_list[gc_content] += 1

    bamfile.close()

    #############################
    ### write gc content file ###
    #############################

    # Write GC content file
    gc_content_file_path = os.path.join(out_dir, f"{pid}_{sample}_gc_content.tsv")
    with open(gc_content_file_path, "w") as gc_content_file:
        gc_content_file.write("gc_content_percent\tread_count\n")
        gc_content_file.writelines(
            f"{gc}\t{count}\n" for gc, count in gc_content_list.items()
        )


####################################################################################################
### estimate the telomere content in telomeric reads per million reads with a similar gc content ###
####################################################################################################


def estimate_telomere_content(
        input_dir,
        out_dir,
        pid,
        sample,
        read_length,
        repeat_threshold_set,
        repeat_threshold_str,
        gc_lower,
        gc_upper,
):
    """
    Estimate telomere content in telomeric reads per million reads with similar GC content.
    """

    # gc bins used for correction
    gc_bins_correction = list(range(gc_lower, gc_upper + 1))

    # read in total gc content counts
    gc_content_file = os.path.join(input_dir, f"{pid}_{sample}_gc_content.tsv")
    gc_content_list = {}

    with open(gc_content_file, "rt") as tsvin:
        tsvin = csv.reader(tsvin, delimiter="\t")
        next(tsvin, None)  # skip the headers
        for row in tsvin:
            gc_content_list[int(row[0])] = int(row[1])

    # get total number of reads in these bins
    sum_over_threshold = sum(gc_content_list[gc] for gc in gc_bins_correction)

    total_read_count = get_total_read_count(input_dir, pid)

    counts_df = count_telomere_read_files(input_dir, pid)

    tel_read_count, intratel_read_count = get_telomere_counts(input_dir, pid)

    # estimate telomere content (intratelomeric reads per million reads with similar gc content)
    tel_content = (
        float(intratel_read_count) / float(sum_over_threshold) * 1000000
        if sum_over_threshold > 0
        else 0.0
    )

    # TRPM
    trpm = (
        float(tel_read_count) / float(total_read_count) * 1000000
        if total_read_count > 0
        else 0.0
    )

    # take special repeat threshold situations into account
    if repeat_threshold_str == "n":
        repeat_threshold_str = "heterogeneous"

    repeat_threshold_set = str(repeat_threshold_set) + " per 100 bp"

    ##########################
    ### write summary file ###
    ##########################

    summary_file_path = os.path.join(out_dir, f"{pid}_{sample}_summary.tsv")

    # Write summary file
    with open(summary_file_path, "w") as summary_file:
        # Header
        header = "\t".join(
            [
                "PID",
                "sample",
                "total_reads",
                "read_length",
                "repeat_threshold_set",
                "repeat_threshold_used",
                "intratelomeric_reads",
                "junctionspanning_reads",
                "subtelomeric_reads",
                "intrachromosomal_reads",
                "tel_read_count",
                "gc_bins_for_correction",
                "total_reads_with_tel_gc",
                "tel_content",
                "TRPM",
            ]
        )
        summary_file.write(header + "\n")

        # Results
        gc_bins = f"{gc_lower}-{gc_upper}"
        results_line = (
            f"{pid}\t{sample}\t{total_read_count}\t{read_length}\t{repeat_threshold_set}\t"
            f"{repeat_threshold_str}\t{int(counts_df['intratelomeric'].values)}\t{int(counts_df['junctionspanning'].values)}\t{int(counts_df['subtelomeric'].values)}\t{int(counts_df['intrachromosomal'].values)}\t{tel_read_count}\t{gc_bins}\t"
            f"{sum_over_threshold}\t{tel_content:.6f}\t{trpm:.6f}\n"
        )
        summary_file.write(results_line)


def get_telomere_counts(input_dir, pid):
    """
    Grab reads_with_patterns from _spectrum.tsv
    """
    spectrum_file = os.path.join(input_dir, f"{pid}_spectrum.tsv")
    try:
        data_tel_reads = numpy.loadtxt(spectrum_file, usecols=2, skiprows=1)

        if data_tel_reads.size == 0:
            return 0, 0

        tel_read_count = int(data_tel_reads.sum())

        # get unmapped_reads count
        if data_tel_reads.size != 1:
            intratel_read_count = int(data_tel_reads[-1])
        else:
            intratel_read_count = int(data_tel_reads.item())

        return tel_read_count, intratel_read_count

    except Exception as e:
        print(f"Error processing file {spectrum_file}: {str(e)}")
        return 0, 0


def get_total_read_count(input_dir, pid):
    # get total read count
    read_count_file = os.path.join(input_dir, f"{pid}_readcount.tsv")
    read_count_array = [
        x.split("\t")[2].strip() for x in open(read_count_file).readlines()
    ]
    read_count_array = [int(x) for x in read_count_array[1: len(read_count_array)]]
    total_read_count = sum(read_count_array)
    return total_read_count


def count_telomere_read_files(input_dir, pid):
    """
    Opens and indexes each of the 4 telomere BAM files, counts their reads,
    and returns a DataFrame with the counts.

    Parameters:
        input_dir (str): Path to the directory containing the BAM files
        pid (str): Process ID or identifier for the BAM files

    Returns:
        pandas.DataFrame: DataFrame containing read counts for each BAM file
    """
    # Define the file categories and their corresponding file names
    categories = [
        "intratelomeric",
        "junctionspanning",
        "subtelomeric",
        "intrachromosomal",
    ]
    file_names = [f"{pid}_filtered_{category}.bam" for category in categories]

    # Initialize a dictionary to store the counts
    counts = {}

    # Process each BAM file
    for category, file_name in zip(categories, file_names):
        file_path = os.path.join(input_dir, file_name)

        try:
            # Check if the file exists
            if not os.path.exists(file_path):
                print(f"Warning: File {file_path} does not exist")
                counts[category] = 0
                continue

            # # Check if the file is indexed; if not, index it
            # index_path = file_path + ".bai"
            # if not os.path.exists(index_path):
            #     print(f"Indexing {file_path}...")
            #     pysam.index(file_path)

            # Open the BAM file and count the reads
            with pysam.AlignmentFile(file_path, "rb") as bamfile:
                # count all reads in the file
                read_count = bamfile.count(until_eof=True)
                counts[category] = read_count
                print(f"{category}: {read_count} reads")

        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            counts[category] = 0

    # Create a DataFrame from the counts dictionary
    df = pd.DataFrame([counts])

    return df
