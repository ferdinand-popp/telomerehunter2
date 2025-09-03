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

from PyPDF2 import PdfMerger, PdfReader

from telomerehunter2.utils import get_band_info

###############################################################################
### merge all PDF plots produced as TelomereHunter output into one PDF file ###
###############################################################################


def merge_telomere_hunter_results(pid, outdir, banding_file=None):
    """
    Merge all TelomereHunter PDF plots into one PDF file.
    """
    # Define the possible PDF file names in correct order for merged PDF
    possible_file_names = [
        f"{pid}_telomere_content",
        f"{pid}_sorted_telomere_read_counts",
        f"{pid}_hist_telomere_repeats_per_intratelomeric_read",
        f"{pid}_gc_content",
        f"{pid}_TVR_barplot",
        f"{pid}_TVR_scatterplot",
        f"{pid}_TVR_plots_absolute_counts",
        f"{pid}_TVR_plots_log2_ratio",
        f"{pid}_pattern_screen_scatterplot",
        f"{pid}_singletons",
    ]

    if banding_file:
        _, chrs = get_band_info(banding_file)
        possible_file_names += [f"{pid}_{chromo}" for chromo in chrs]

    possible_pdf_names = [f"{name}.pdf" for name in possible_file_names]

    # Check which of the possible PDF files exist and sort
    files_dir = os.path.join(outdir, "plots")
    pdf_files = [f for f in os.listdir(files_dir) if f.endswith("pdf")]
    pdf_files_ordered = [pdf for pdf in possible_pdf_names if pdf in pdf_files]

    # Merge files
    merger = PdfMerger()

    for filename in pdf_files_ordered:
        pdf_path = os.path.join(files_dir, filename)
        merger.append(PdfReader(pdf_path, "rb"))

    # Write the merged PDF file
    merged_pdf_path = os.path.join(outdir, f"{pid}_all_plots_merged.pdf")
    merger.write(merged_pdf_path)

    print(f"Merged PDFs saved to: {merged_pdf_path}")

    # Merge all possible HTML files # maybe bring in order like pdfs
    possible_html_names = [f"{name}.html" for name in possible_file_names]
    html_dir = os.path.join(outdir, "html_reports")
    html_files = [f for f in os.listdir(html_dir) if f.endswith("html")]
    html_files_ordered = [html for html in possible_html_names if html in html_files]

    all_html_path = os.path.join(outdir, "all_reports.html")
    with open(all_html_path, "w") as output_html:
        for html_file in html_files_ordered:
            html_path = os.path.join(html_dir, html_file)
            with open(html_path, "r") as input_html:
                output_html.write(input_html.read())

    print(f"Merged HTML reports to: {all_html_path}")
