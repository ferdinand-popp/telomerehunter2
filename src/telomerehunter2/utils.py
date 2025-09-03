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

import cProfile
import io
import os
import pstats
import time

import pandas as pd


def measure_time(func):
    """Wrapper function to simply print the execution time for a function"""

    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time for {func.__name__}: {elapsed_time} seconds")
        return result

    return wrapper


def profile_function(func, *args, **kwargs):
    """Wrapper function to profile a function"""
    pr = cProfile.Profile()
    pr.enable()
    result = func(*args, **kwargs)
    pr.disable()
    s = io.StringIO()
    sortby = "cumulative"
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(5)  # This will print only the top 10 entries
    print(s.getvalue())
    return result


def get_reverse_complement(dna_sequence):
    """Get the reverse complement of a DNA sequence."""
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    complemented_sequence = (complement_dict[base] for base in reversed(dna_sequence))
    return "".join(complemented_sequence)


def get_band_info(banding_file):
    # Get band lengths
    band_info = pd.read_table(
        banding_file, sep="\t", names=["chr", "start", "end", "band_name", "stain"]
    )  # maybe not specify col names?
    band_info["chr"] = band_info["chr"].str.replace("chr", "")
    band_info["length"] = band_info["end"] - band_info["start"]
    chrs = band_info["chr"].unique()
    sorted_chromosomes = sorted(
        chrs, key=lambda x: (int(x) if x.isdigit() else float("inf"), x)
    )
    return band_info, sorted_chromosomes


def assure_dir_exists(path):
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)
