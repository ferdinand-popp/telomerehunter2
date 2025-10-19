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

from telomerehunter2.utils import get_band_info, get_reverse_complement


def make_band_and_spectrum_lists(band_file, patterns):
    """
    Create band and spectrum lists based on the information from a banding file.

    Parameters:
    - band_file (str): Path to the banding file.
    - patterns (list): List of patterns for spectrum initialization.

    Returns:
    - bands_list (dict): Dictionary containing band information for each chromosome.
    - spectrum_list (dict): Dictionary containing spectrum information for each chromosome and band.
    - chromosome_list (list): List of chromosome names.

    """
    # Get chromosome information
    _, chromosome_list = get_band_info(band_file)

    # Initialize bands_list with unmapped entry
    bands_list = {"unmapped": {"band_name": ["unmapped"], "end": [0]}}
    bands_list.update({chr: {"band_name": [], "end": []} for chr in chromosome_list})

    # Initialize spectrum_list dictionary
    spectrum_list = {
        chr: {
            band_name: {
                "other": 0.0,
                "reads_with_pattern": 0,
                **{pattern: 0 for pattern in patterns},
            }
            for band_name in ["junction1", "junction2"]
        }
        for chr in chromosome_list
    }
    spectrum_list["unmapped"] = {
        "unmapped": {
            "other": 0.0,
            "reads_with_pattern": 0,
            **{pattern: 0 for pattern in patterns},
        }
    }

    # Process banding file
    with open(band_file, "r") as band_file:
        for line in band_file:
            try:
                # Parse banding file line
                chrom, _, end, band_name, *additional_items = line.strip().split()
                chrom = chrom[3:] if chrom.startswith("chr") else chrom

                # Update bands_list with band information
                bands_list[chrom]["band_name"].append(band_name)
                bands_list[chrom]["end"].append(int(end))

                # Initialize spectrum_list for the current band
                spectrum_list[chrom][band_name] = {
                    "other": 0.0,
                    "reads_with_pattern": 0,
                    **{pattern: 0 for pattern in patterns},
                }

            except Exception as e:
                print(f"Invalid line in banding file: '{line}', {e}")

    return bands_list, spectrum_list, chromosome_list


def write_spectrum_file(
    out_dir, pid, chromosome_list, bands_list, spectrum_list, patterns
):
    spectrum_file_path = os.path.join(out_dir, f"{pid}_spectrum.tsv")

    if not spectrum_list:
        spectrum_list = get_default_spectrum_list(patterns, spectrum_list)

    with open(spectrum_file_path, "w") as spectrum_file:
        spectrum_file.write("chr\tband\treads_with_pattern")

        for pattern in patterns:
            spectrum_file.write("\t" + pattern)

        spectrum_file.write("\tother\n")

        for chromosome in chromosome_list:
            for band in (
                ["junction1"] + bands_list[chromosome]["band_name"] + ["junction2"]
            ):
                spectrum_file.write(
                    f"{chromosome}\t{band}\t{spectrum_list[chromosome][band]['reads_with_pattern']}"
                )

                for pattern in patterns:
                    spectrum_file.write(f"\t{spectrum_list[chromosome][band][pattern]}")

                spectrum_file.write(
                    f"\t{int(round(spectrum_list[chromosome][band]['other']))}\n"
                )

        spectrum_file.write(
            f"unmapped\tunmapped\t{spectrum_list['unmapped']['unmapped']['reads_with_pattern']}"
        )

        for pattern in patterns:
            spectrum_file.write(f"\t{spectrum_list['unmapped']['unmapped'][pattern]}")

        spectrum_file.write(
            f"\t{int(round(spectrum_list['unmapped']['unmapped']['other']))}\n"
        )


def get_default_spectrum_list(patterns, spectrum_list):
    spectrum_list = {
        "unmapped": {
            "unmapped": {
                "reads_with_pattern": 0,
                "other": 0,
                **{pattern: 0 for pattern in patterns},
            }
        }
    }
    return spectrum_list


def close_file_handles(
    bamfile,
    intratelomeric_file,
    junctionspanning_file,
    subtelomeric_file,
    intrachromosomal_file,
):
    bamfile.close()
    intratelomeric_file.close()
    junctionspanning_file.close()
    subtelomeric_file.close()
    intrachromosomal_file.close()


# function for checking reads without mate and sorting
def read_check(read, references, bands_list, chromosome_list, mapq_threshold):
    """
    Check properties of a read and return relevant information.

    Parameters:
    - read: object, BAM read object.
    - references: list, list of reference sequences.
    - bands_list: dict, dictionary containing band information for each chromosome.
    - chromosome_list: list, list of valid chromosome names.
    - mapq_threshold: int, minimum mapping quality threshold.

    Returns:
    Tuple containing information about the read:
    (chromosome, band, read_is_unmapped, read_junctionspanning, read_junction)
    """

    # get reference
    read_tid = read.tid
    read_ref = references[read_tid].replace("chr", "") if read_tid != -1 else ""

    read_junctionspanning = False
    read_junction = ""

    # check if read is considered unmapped
    if (
        read.is_unmapped
        or read_ref not in chromosome_list
        or read.mapping_quality < mapq_threshold
    ):
        read_is_unmapped = True
        chromosome, band = "unmapped", "unmapped"
    # map read
    else:
        read_is_unmapped = False
        chromosome = read_ref

        # find band
        i = 0
        while read.pos > bands_list[chromosome]["end"][i] and i < (
            len(bands_list[chromosome]["end"]) - 1
        ):
            i += 1
        else:
            band = bands_list[chromosome]["band_name"][i]

            # check if read is mapped to first or last chromosome band
            if i == 0:
                read_junctionspanning = True
                read_junction = "junction1"
            elif i == (len(bands_list[chromosome]["end"]) - 1):
                read_junctionspanning = True
                read_junction = "junction2"

    return chromosome, band, read_is_unmapped, read_junctionspanning, read_junction


# function for sorting reads without a mate into correct fraction and add counts to spectrum list
def sort_reads_without_mate(
    read,
    references,
    bands_list,
    chromosome_list,
    mapq_threshold,
    spectrum_list,
    patterns,
    intratelomeric_file,
    subtelomeric_file,
    intrachromosomal_file,
):
    """
    Sort reads without a mate into the correct fraction and update counts.

    Parameters:
    - read: object, BAM read object.
    - references: list, list of reference sequences.
    - bands_list: dict, dictionary containing band information for each chromosome.
    - chromosome_list: list, list of valid chromosome names.
    - mapq_threshold: int, minimum mapping quality threshold.
    - spectrum_list: dict, dictionary to store counts of different patterns for each chromosome and band.
    - patterns: list, with base patterns
    - intratelomeric_file: file, output file for intratelomeric reads.
    - subtelomeric_file: file, output file for subtelomeric reads.
    - intrachromosomal_file: file, output file for intrachromosomal reads.
    """

    # do read check
    chromosome, band, read_is_unmapped, read_junctionspanning, read_junction = (
        read_check(read, references, bands_list, chromosome_list, mapq_threshold)
    )

    # write read to correct fraction
    if read_is_unmapped:
        intratelomeric_file.write(read)
    elif read_junctionspanning:
        subtelomeric_file.write(read)
    else:
        intrachromosomal_file.write(read)

    # add read counts to spectrum list
    spectrum_list[chromosome][band]["reads_with_pattern"] += 1

    read_total_pattern_count = 0

    for pattern in patterns:
        pattern_count = read.query_sequence.count(pattern)
        spectrum_list[chromosome][band][pattern] += pattern_count
        read_total_pattern_count += pattern_count

    # get count of other e.g. hexamers found
    k_mer_length = len(patterns[0])
    spectrum_list[chromosome][band]["other"] += (
        float(len(read.query_sequence)) / k_mer_length
    ) - read_total_pattern_count


# function for writing reads with mates to fraction and adding counts to spectrum list
def sort_read_with_mate(read, fraction_file, chromosome, band, spectrum_list, patterns):
    """
    Sort reads with a mate into the correct fraction and update counts.

    Parameters:
    - read: object, BAM read object.
    - fraction_file: file, output file for the specific fraction.
    - chromosome: str, chromosome name.
    - band: str, band name.
    - spectrum_list: dict, dictionary to store counts of different patterns for each chromosome and band.
    - patterns: list of target base patterns
    """

    fraction_file.write(read)

    spectrum_list[chromosome][band]["reads_with_pattern"] += 1

    read_total_pattern_count = 0

    # count occurrence of pattern types
    for pattern in patterns:
        pattern_count = read.query_sequence.count(pattern)
        spectrum_list[chromosome][band][pattern] += pattern_count
        read_total_pattern_count += pattern_count

    # get count of other e.g. hexamers found
    k_mer_length = len(patterns[0])
    spectrum_list[chromosome][band]["other"] += (
        float(len(read.query_sequence)) / k_mer_length
    ) - read_total_pattern_count


def process_reads(
    bamfile,
    chromosome_list,
    bands_list,
    spectrum_list,
    patterns,
    mapq_threshold,
    intratelomeric_file,
    junctionspanning_file,
    subtelomeric_file,
    intrachromosomal_file,
):
    """
    Process paired-end reads and categorize them into different fractions based on mapping characteristics.

    Parameters:
    - bamfile: iterator, iterator over BAM file containing paired-end reads.
    - chromosome_list: list, list of valid chromosome names.
    - bands_list: dict, dictionary containing band information for each chromosome.
    - spectrum_list: dict, dictionary to store counts of different patterns for each chromosome and band.
    - patterns: list, target base sequences
    - mapq_threshold: int, minimum mapping quality threshold.
    - intratelomeric_file: file, output file for intratelomeric reads.
    - junctionspanning_file: file, output file for junction spanning reads.
    - subtelomeric_file: file, output file for subtelomeric reads.
    - intrachromosomal_file: file, output file for intrachromosomal reads.

    Returns:
    Tuple: chromosome_list, bands_list, spectrum_list
    """

    break_flag = False
    references = bamfile.references

    while True:
        try:
            read1 = next(bamfile)
        except StopIteration:
            break

        try:
            read2 = next(bamfile)
        except StopIteration:
            sort_reads_without_mate(
                read1,
                references,
                bands_list,
                chromosome_list,
                mapq_threshold,
                spectrum_list,
                patterns,
                intratelomeric_file,
                subtelomeric_file,
                intrachromosomal_file,
            )
            break

        read1_first_flag = True

        # while not mates
        while read1.query_name != read2.query_name:
            if read1_first_flag:
                sort_reads_without_mate(
                    read1,
                    references,
                    bands_list,
                    chromosome_list,
                    mapq_threshold,
                    spectrum_list,
                    patterns,
                    intratelomeric_file,
                    subtelomeric_file,
                    intrachromosomal_file,
                )
                try:
                    read1 = next(bamfile)
                except StopIteration:
                    sort_reads_without_mate(
                        read2,
                        references,
                        bands_list,
                        chromosome_list,
                        mapq_threshold,
                        spectrum_list,
                        patterns,
                        intratelomeric_file,
                        subtelomeric_file,
                        intrachromosomal_file,
                    )
                    break_flag = True
                    break
                read1_first_flag = False
            else:
                sort_reads_without_mate(
                    read2,
                    references,
                    bands_list,
                    chromosome_list,
                    mapq_threshold,
                    spectrum_list,
                    patterns,
                    intratelomeric_file,
                    subtelomeric_file,
                    intrachromosomal_file,
                )
                try:
                    read2 = next(bamfile)
                except StopIteration:
                    sort_reads_without_mate(
                        read1,
                        references,
                        bands_list,
                        chromosome_list,
                        mapq_threshold,
                        spectrum_list,
                        patterns,
                        intratelomeric_file,
                        subtelomeric_file,
                        intrachromosomal_file,
                    )
                    break_flag = True
                    break
                read1_first_flag = True

        if break_flag:
            break

        # Complete Mates
        (
            read1_chromosome,
            read1_band,
            read1_is_unmapped,
            read1_junctionspanning,
            read1_junction,
        ) = read_check(
            read1, bamfile.references, bands_list, chromosome_list, mapq_threshold
        )
        (
            read2_chromosome,
            read2_band,
            read2_is_unmapped,
            read2_junctionspanning,
            read2_junction,
        ) = read_check(
            read2, bamfile.references, bands_list, chromosome_list, mapq_threshold
        )

        # INTRATELOMERIC: both mates are unmapped
        if read1_is_unmapped and read2_is_unmapped:
            sort_read_with_mate(
                read=read1,
                fraction_file=intratelomeric_file,
                chromosome="unmapped",
                band="unmapped",
                spectrum_list=spectrum_list,
                patterns=patterns,
            )
            sort_read_with_mate(
                read=read2,
                fraction_file=intratelomeric_file,
                chromosome="unmapped",
                band="unmapped",
                spectrum_list=spectrum_list,
                patterns=patterns,
            )

        # JUNCTION SPANNING: one mate is unmapped and the other mate is mapped to first or last band of chromosome
        elif (read1_is_unmapped and read2_junctionspanning) or (
            read2_is_unmapped and read1_junctionspanning
        ):
            if read1_is_unmapped:
                chromosome = read2_chromosome
                band = read2_junction
            else:
                chromosome = read1_chromosome
                band = read1_junction

            sort_read_with_mate(
                read=read1,
                fraction_file=junctionspanning_file,
                chromosome=chromosome,
                band=band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )
            sort_read_with_mate(
                read=read2,
                fraction_file=junctionspanning_file,
                chromosome=chromosome,
                band=band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )

        # SUBTELOMERIC: both mates are mapped to first or last band or chromosome (and have similar mapping positions)
        elif read1_junctionspanning and read2_junctionspanning:
            sort_read_with_mate(
                read=read1,
                fraction_file=subtelomeric_file,
                chromosome=read1_chromosome,
                band=read1_band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )
            sort_read_with_mate(
                read=read2,
                fraction_file=subtelomeric_file,
                chromosome=read2_chromosome,
                band=read2_band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )

        # SUBTELOMERIC/INTRACHROMOSOMAL: one read is subtelomeric, the other is intra-chromosomal => count and sort individually
        elif read1_junctionspanning or read2_junctionspanning:
            if read1_junctionspanning:
                read1_file = subtelomeric_file
                read2_file = intrachromosomal_file
            else:
                read2_file = subtelomeric_file
                read1_file = intrachromosomal_file

            sort_read_with_mate(
                read=read1,
                fraction_file=read1_file,
                chromosome=read1_chromosome,
                band=read1_band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )
            sort_read_with_mate(
                read=read2,
                fraction_file=read2_file,
                chromosome=read2_chromosome,
                band=read2_band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )

        # INTRACHROMOSOMAL: one read is intra-chromosomal, the other is unmapped => count both to position of mapped read
        elif (read1_is_unmapped and not read2_junctionspanning) or (
            read2_is_unmapped and not read1_junctionspanning
        ):
            if read1_is_unmapped:
                chromosome = read2_chromosome
                band = read2_band
            else:
                chromosome = read1_chromosome
                band = read1_band

            sort_read_with_mate(
                read=read1,
                fraction_file=intrachromosomal_file,
                chromosome=chromosome,
                band=band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )
            sort_read_with_mate(
                read=read2,
                fraction_file=intrachromosomal_file,
                chromosome=chromosome,
                band=band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )

        # INTRACHROMOSOMAL: both mapped
        else:
            sort_read_with_mate(
                read=read1,
                fraction_file=intrachromosomal_file,
                chromosome=read1_chromosome,
                band=read1_band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )
            sort_read_with_mate(
                read=read2,
                fraction_file=intrachromosomal_file,
                chromosome=read2_chromosome,
                band=read2_band,
                spectrum_list=spectrum_list,
                patterns=patterns,
            )

    return chromosome_list, bands_list, spectrum_list


def write_spectrum_file_no_banding(out_dir, pid, spectrum_list, patterns):
    spectrum_file_path = os.path.join(out_dir, f"{pid}_spectrum.tsv")

    if not spectrum_list:
        spectrum_list = get_default_spectrum_list(patterns, spectrum_list)

    with open(spectrum_file_path, "w") as spectrum_file:
        spectrum_file.write("chr\tband\treads_with_pattern")

        for pattern in patterns:
            spectrum_file.write("\t" + pattern)

        spectrum_file.write("\tother\n")

        for ref in spectrum_list:
            for band in spectrum_list[ref]:
                spectrum_file.write(
                    f"{ref}\t{band}\t{spectrum_list[ref][band]['reads_with_pattern']}"
                )

                for pattern in patterns:
                    spectrum_file.write(f"\t{spectrum_list[ref][band][pattern]}")

                spectrum_file.write(
                    f"\t{int(round(spectrum_list[ref][band]['other']))}\n"
                )


def count_patterns_in_read_no_banding(read, mapq_threshold, patterns, spectrum_list):
    """
    Counts the occurrences of specified patterns in a read and updates the spectrum list.

    Args:
        read (pysam.AlignedSegment): Aligned read object.
        mapq_threshold (int): Minimum mapping quality threshold for considering a read mapped.
        patterns (list): List of patterns to count occurrences of.
        spectrum_list (dict): Dictionary containing counts of patterns for different reference sequences and bands.

    Returns:
        dict: Updated spectrum list after counting occurrences of patterns in the read.
    """

    # Unmapped or low qual reads should be "unmapped"
    is_unmapped = read.is_unmapped or read.mapping_quality < mapq_threshold
    ref = "unmapped" if is_unmapped else read.reference_name
    band = "unmapped" if is_unmapped else "_"

    read_total_pattern_count = 0

    # Intialize the dict counts (patterns + other) to store other hexamers counts
    if ref not in spectrum_list:
        spectrum_list[ref] = {band: {"other": 0, "reads_with_pattern": 0}}

    # count occurrence of pattern types
    for pattern in patterns:
        if pattern not in spectrum_list[ref][band]:
            spectrum_list[ref][band][pattern] = 0

        pattern_count = read.query_sequence.count(pattern)
        spectrum_list[ref][band][pattern] += pattern_count
        read_total_pattern_count += pattern_count

    # get count of other e.g. hexamers found
    k_mer_length = len(patterns[0])
    spectrum_list[ref][band]["other"] += (
        float(len(read.query_sequence)) / k_mer_length
    ) - read_total_pattern_count

    # count in general the read with count
    spectrum_list[ref][band]["reads_with_pattern"] += 1

    return spectrum_list


def sort_read_without_mate_no_banding(
    read,
    mapq_threshold,
    patterns,
    intratelomeric_file,
    intrachromosomal_file,
    spectrum_list,
):
    is_mapped_read = read.is_mapped and read.mapping_quality >= mapq_threshold
    # Categorize reads and write them to respective files
    if not is_mapped_read:
        # not mapped -> telo
        intratelomeric_file.write(read)
    else:
        # mapped -> chr
        intrachromosomal_file.write(read)

    # Count patterns in read sequences
    spectrum_list = count_patterns_in_read_no_banding(
        read, mapq_threshold, patterns, spectrum_list
    )

    return spectrum_list


def process_reads_no_banding(
    bamfile, patterns, mapq_threshold, intratelomeric_file, intrachromosomal_file
):
    break_flag = False
    # dict with levels: reference_name, band, count
    spectrum_list = {}

    while True:
        try:
            read1 = next(bamfile)
        except StopIteration:
            break

        try:
            read2 = next(bamfile)
        except StopIteration:
            sort_read_without_mate_no_banding(
                read1,
                mapq_threshold,
                patterns,
                intratelomeric_file,
                intrachromosomal_file,
                spectrum_list,
            )
            break

        read1_first_flag = True

        # while not mates iterate over read 1 and 2 step by step
        while read1.query_name != read2.query_name:
            # exchange read 1
            if read1_first_flag:
                sort_read_without_mate_no_banding(
                    read1,
                    mapq_threshold,
                    patterns,
                    intratelomeric_file,
                    intrachromosomal_file,
                    spectrum_list,
                )
                try:
                    read1 = next(bamfile)
                except StopIteration:
                    sort_read_without_mate_no_banding(
                        read2,
                        mapq_threshold,
                        patterns,
                        intratelomeric_file,
                        intrachromosomal_file,
                        spectrum_list,
                    )
                    break_flag = True
                    break
                read1_first_flag = False
            # exchange read 2
            else:
                sort_read_without_mate_no_banding(
                    read2,
                    mapq_threshold,
                    patterns,
                    intratelomeric_file,
                    intrachromosomal_file,
                    spectrum_list,
                )
                try:
                    read2 = next(bamfile)
                except StopIteration:
                    sort_read_without_mate_no_banding(
                        read1,
                        mapq_threshold,
                        patterns,
                        intratelomeric_file,
                        intrachromosomal_file,
                        spectrum_list,
                    )
                    break_flag = True
                    break
                read1_first_flag = True

        if break_flag:
            break

        # True mates

        # Check mapping quality and determine if reads are mapped
        is_mapped_read1 = read1.is_mapped and read1.mapping_quality >= mapq_threshold
        is_mapped_read2 = read2.is_mapped and read2.mapping_quality >= mapq_threshold

        # Categorize reads and write them to respective files
        if not is_mapped_read1 and not is_mapped_read2:
            # both not mapped -> telomeric
            intratelomeric_file.write(read1)
            intratelomeric_file.write(read2)
        else:
            # one or both mapped -> chromsomal
            intrachromosomal_file.write(read1)
            intrachromosomal_file.write(read2)

        # Count patterns in read sequences
        spectrum_list = count_patterns_in_read_no_banding(
            read1, mapq_threshold, patterns, spectrum_list
        )
        spectrum_list = count_patterns_in_read_no_banding(
            read2, mapq_threshold, patterns, spectrum_list
        )

    return spectrum_list


def sort_telomere_reads(input_dir, band_file, out_dir, pid, mapq_threshold, repeats):
    """
    Sorts the reads from an input BAM file (needs to be name sorted!) containing telomere reads into 4 different fractions:
    intratelomeric, junction-spanning, subtelomeric, intrachromosomal. Additionally, a table containing the number of
    telomere reads and the pattern occurrences per chromosome band is written.

    Parameters:
        input_dir (str): Path to the directory containing input BAM file.
        band_file (str): Path to the file containing chromosome band information. Pass None if not available.
        out_dir (str): Path to the directory where output files will be written.
        pid (str): Process ID or identifier for the output files.
        mapq_threshold (int): Mapping quality threshold for filtering reads.
        repeats (list): List of telomere repeat sequences.

    Returns:
        None
    """

    #####################
    ### get patterns  ###
    #####################

    # Extract patterns and reverse complements
    patterns = []
    for repeat in repeats:
        patterns.append(repeat)
        patterns.append(get_reverse_complement(repeat))

    # Open input BAM file for reading
    bamfile = pysam.AlignmentFile(
        os.path.join(input_dir, f"{pid}_filtered_name_sorted.bam"), "rb"
    )

    # Open filtered files for writing
    intratelomeric_file = pysam.AlignmentFile(
        os.path.join(out_dir, f"{pid}_filtered_intratelomeric.bam"),
        "wb",
        template=bamfile,
    )
    junctionspanning_file = pysam.AlignmentFile(
        os.path.join(out_dir, f"{pid}_filtered_junctionspanning.bam"),
        "wb",
        template=bamfile,
    )
    subtelomeric_file = pysam.AlignmentFile(
        os.path.join(out_dir, f"{pid}_filtered_subtelomeric.bam"),
        "wb",
        template=bamfile,
    )
    intrachromosomal_file = pysam.AlignmentFile(
        os.path.join(out_dir, f"{pid}_filtered_intrachromosomal.bam"),
        "wb",
        template=bamfile,
    )

    #################################################################
    ### go through name sorted BAM file containing telomere reads ###
    ### => always look at mate pairs (if possible)                ###
    #################################################################

    # run classic version
    if band_file:
        # Make band and spectrum lists
        bands_list, spectrum_list, chromosome_list = make_band_and_spectrum_lists(
            band_file, patterns
        )

        # Process reads writing to the files
        chromosome_list, bands_list, spectrum_list = process_reads(
            bamfile,
            chromosome_list,
            bands_list,
            spectrum_list,
            patterns,
            mapq_threshold,
            intratelomeric_file,
            junctionspanning_file,
            subtelomeric_file,
            intrachromosomal_file,
        )

        # Write spectrum file
        write_spectrum_file(
            out_dir, pid, chromosome_list, bands_list, spectrum_list, patterns
        )

    # run without banding information (fast)
    else:
        # Process reads writing to the files, taking references from reads
        spectrum_list = process_reads_no_banding(
            bamfile,
            patterns,
            mapq_threshold,
            intratelomeric_file,
            intrachromosomal_file,
        )

        # Write spectrum file slim version
        write_spectrum_file_no_banding(out_dir, pid, spectrum_list, patterns)

    # Close file handles
    close_file_handles(
        bamfile,
        intratelomeric_file,
        junctionspanning_file,
        subtelomeric_file,
        intrachromosomal_file,
    )
