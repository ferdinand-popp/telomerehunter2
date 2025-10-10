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
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

pio.kaleido.scope.mathjax = None

from telomerehunter2.utils import get_band_info, get_reverse_complement

TVR_COLORS = {
    "TTAGGG": "#CB181D",
    "CCCTAA": "#FF656A",
    "TGAGGG": "#74C476",
    "CCCTCA": "#C1FFC3",
    "TCAGGG": "#2171B5",
    "CCCTGA": "#6EBEFF",
    "TTGGGG": "#FFA500",
    "CCCCAA": "#FFF24D",
    "TTCGGG": "#9370DB",
    "CCCGAA": "#E0BDFF",
    "TTTGGG": "#FF00FF",
    "CCCAAA": "#00FFFF",
    "other": "#000000",
}

PLOT_CONFIG = {"width": 1000, "font_family": "Arial", "font_size": 12}


#    'height': 500,
#    'title_font_size': 16,
# 'xaxis_showgrid': True,
# 'xaxis_gridcolor': 'LightGray',
# 'yaxis_showgrid': True,
# 'yaxis_gridcolor': 'LightGray',
# 'paper_bgcolor': 'white',
# 'plot_bgcolor': 'white',


def save_plot(
    fig, outdir, name, formats="pdf", width=PLOT_CONFIG["width"], height=None
):
    """Save plot in multiple formats."""
    for fmt in formats:
        if fmt == "html":
            # Save interactive HTML version
            html_path = os.path.join(outdir, "html_reports", f"{name}.html")
            try:
                fig.write_html(html_path)
            except Exception as e:
                print(f"Error saving HTML plot: {str(e)}")
        else:
            # Save static image
            output_path = os.path.join(outdir, "plots", f"{name}.{fmt}")
            try:
                (
                    fig.write_image(output_path, width=width)
                    if not height
                    else fig.write_image(output_path, width=width, height=height)
                )
            except Exception as e:
                print(f"Error saving plot in {fmt} format: {str(e)}")



def validate_data(df, required_columns):
    """Validate dataframe has required columns and data."""
    if df is None or df.empty:
        return False
    return all(col in df.columns for col in required_columns)


def combine_height_df(df1, df2):
    """merging two dataframes outer tumor and control"""

    return pd.merge(
        df1,
        df2,
        how="outer",
        left_index=True,
        right_index=True,
        suffixes=("", "_right_control"),
    )


def get_sample_types_from_spectrum(outdir, pid):
    # Get samples
    spectrum_tumor_file = os.path.join(
        outdir, "tumor_TelomerCnt_" + pid, pid + "_spectrum.tsv"
    )
    spectrum_control_file = os.path.join(
        outdir, "control_TelomerCnt_" + pid, pid + "_spectrum.tsv"
    )
    samples = []
    if os.path.exists(spectrum_tumor_file):
        samples.append("tumor")
    if os.path.exists(spectrum_control_file):
        samples.append("control")
    return samples


def combine_forward_patterns(plot_rev_complement, df, group_by_cols):
    # Combine forward and reverse patterns
    if not plot_rev_complement:
        df["repeat_type_forward"] = df["repeat_type"].apply(
            lambda x: x if "GGG" in x or x == "other" else get_reverse_complement(x)
        )
        df = df.groupby(group_by_cols).sum().reset_index()
        df = df.drop(columns=["repeat_type"], errors="ignore")
        df = df.rename(columns={"repeat_type_forward": "repeat_type"})

    df = df.rename(columns={"value": "sum"})
    return df


def read_frequency_table(file_path, sample_type):
    if os.path.exists(file_path):
        frequency_table = pd.read_table(file_path)
        frequency_table = frequency_table.sort_values(by="number_repeats")
        frequency_table["sample"] = sample_type
        frequency_table["percent"] = (
            frequency_table["count"] / frequency_table["count"].sum() * 100
        )
        frequency_table["percent_cumulative"] = frequency_table["percent"].cumsum()
    else:
        frequency_table = pd.DataFrame(
            columns=[
                "count",
                "number_repeats",
                "percent",
                "percent_cumulative",
                "sample",
            ]
        )

    return frequency_table


def update_repeat_types_order_and_col(df, colors):
    # kick unused TVRs
    colors = {
        key: value for key, value in colors.items() if key in df["repeat_type"].unique()
    }

    # get unmatched TVRs
    repeats_no_color = df.loc[
        ~df["repeat_type"].isin(colors.keys()), "repeat_type"
    ].unique()

    if len(repeats_no_color) > 0:
        # make possible TVRs categorical
        possible_TVRs = list(colors.keys())[:-1] + list(repeats_no_color) + ["other"]
        df["repeat_type"] = pd.Categorical(
            df["repeat_type"], categories=possible_TVRs, ordered=True
        )

        # add colors for unmatched TVRs
        palette = px.colors.sequential.Greys
        colors.update(
            {
                repeat: color
                for repeat, color in zip(repeats_no_color, reversed(palette[1:-1]))
            }
        )

    return df, colors


def barplot_repeattype(
    height_matrix,
    outdir,
    plot_name,
    plot_file_format=["pdf", "png", "svg"],
    main="",
    sub_title="",
    ylab="Telomere Reads (per Million Reads)",
    plot_reverse_complement=False,
    repeat_threshold="",
    count_type="",
    mapq_threshold="",
    control_axis=True,
    labels=None,
):
    # sum complements if not separate plotting
    if not plot_reverse_complement:
        height_patterns = height_matrix.iloc[:-1, :]
        height_summed = height_patterns.groupby(
            np.arange(len(height_patterns)) // 2
        ).sum()
        first_row_names = height_patterns.index[::2]
        height_summed.index = first_row_names
        height_summed.loc["other"] = height_matrix.iloc[-1, :]
    else:
        height_summed = height_matrix

    # handle control case: grouped stacked barchart is special in plotly!!
    if control_axis:
        # Get the number of columns and select half of them
        half_num_cols = len(height_summed.columns) // 2

        # Split DataFrame into two based on the number of columns
        height_summed_control = height_summed.iloc[:, half_num_cols:]
        height_summed_tumor = height_summed.iloc[:, :half_num_cols]

        height_summed_control.columns = height_summed_control.columns.str.replace(
            "_right_control", "", regex=False
        )

        # Create x-axis positions for bars
        x_positions = np.arange(len(height_summed_tumor.columns))
        bar_width = 0.35

        fig = go.Figure(
            layout=go.Layout(
                barmode="stack",  # Change to stack mode
                yaxis_showticklabels=True,
                yaxis_range=[0, height_matrix.max().max() * 2.5],
                hovermode="x",
                margin=dict(
                    b=100, t=50, l=0, r=10
                ),  # Increased bottom margin for legend
                title=dict(
                    text=main
                    + f"<br><sup>{repeat_threshold} {count_type} Repeats, Mapq Threshold = {mapq_threshold}</sup>",
                    x=0.5,
                ),
                xaxis_title=sub_title,
                yaxis_title=ylab,
            )
        )

        # Define colors for the TVRs (assuming TVR_COLORS is defined)
        colors = {TVR: TVR_COLORS.get(TVR, "#000000") for TVR in height_matrix.index}

        # Add traces for tumor and control data
        for tvr in height_summed_tumor.index:
            # Tumor bars (shifted left)
            fig.add_trace(
                go.Bar(
                    x=[x - bar_width / 2 for x in x_positions],
                    y=height_summed_tumor.loc[tvr],
                    name=tvr,
                    marker_color=colors[tvr],
                    showlegend=True,
                )
            )

            # Control bars (shifted right)
            fig.add_trace(
                go.Bar(
                    x=[x + bar_width / 2 for x in x_positions],
                    y=height_summed_control.loc[tvr],
                    name=tvr,
                    marker_color=colors[tvr],
                    showlegend=False,
                )
            )

        # Update x-axis with smaller font size and proper labels
        if labels:
            fig.update_xaxes(
                tickangle=45,
                tickmode="array",
                tickvals=x_positions,
                ticktext=labels,
                ticks="outside",
                tickfont=dict(size=10),
            )
        else:
            fig.update_xaxes(
                tickangle=45,
                tickmode="array",
                tickvals=x_positions,
                ticktext=height_summed_tumor.columns,
                ticks="outside",
                tickfont=dict(size=10),
            )

        # Update layout for better readability and legend placement
        fig.update_layout(
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.3,  # Move legend below plot
                xanchor="center",
                x=0.5,  # Center the legend
                traceorder="reversed",
                font=dict(size=10),
                bgcolor="rgba(255,255,255,0.5)",
                bordercolor="rgba(0,0,0,0.5)",
                borderwidth=1,
            ),
            bargap=0,
            bargroupgap=0,
        )

    else:
        # Rest of the code for non-control case
        x = height_summed.columns
        y = height_summed.values
        colors = [TVR_COLORS[TVR] for TVR in height_summed.index]

        traces = []
        for i in range(len(y)):
            trace = go.Bar(
                x=x, y=y[i], name=height_summed.index[i], marker=dict(color=colors[i])
            )
            traces.append(trace)

        layout = go.Layout(
            title=main
            + f"<br><sup>{repeat_threshold} {count_type} Repeats, Mapq Threshold = {mapq_threshold}</sup>",
            barmode="stack",
            xaxis=dict(title=sub_title, tickangle=45, tickfont=dict(size=10)),
            yaxis=dict(title=ylab, rangemode="tozero"),
            margin=dict(b=100, t=50, l=0, r=10),  # Increased bottom margin for legend
        )

        fig = go.Figure(data=traces, layout=layout)

        fig.update_layout(
            legend=dict(
                orientation="h",
                yanchor="top",
                y=-0.3,  # Move legend below plot
                xanchor="center",
                x=0.5,  # Center the legend
                traceorder="reversed",
                font=dict(size=10),
                bgcolor="rgba(255,255,255,0.5)",
                bordercolor="rgba(0,0,0,0.5)",
                borderwidth=1,
            ),
            **PLOT_CONFIG,
        )
        fig.update_xaxes(tickmode="linear", ticks="outside")

    # Save figure
    save_plot(fig, outdir, plot_name, plot_file_format)


def plot_spectrum(
    outdir,
    pid,
    repeat_threshold_plot,
    consecutive,
    mapq_threshold,
    banding_file,
    plot_rev_complement,
    plot_file_format,
):
    # Tumor, control or both
    samples = get_sample_types_from_spectrum(outdir, pid)

    band_info, chrs = get_band_info(banding_file)

    spectrum_list = {}

    for sample in samples:
        # Get spectrum
        spectrum_file = os.path.join(
            outdir, f"{sample}_TelomerCnt_{pid}", f"{pid}_spectrum.tsv"
        )

        try:
            spectrum = pd.read_table(spectrum_file, sep="\t", comment="#")
        except FileNotFoundError:
            print(f"Didnt find {spectrum_file}")
            continue

        # Check if the file has no rows, if so, continue to the next iteration
        if spectrum.empty:
            print(f"Empty: {spectrum_file}")
            continue

        # Get total number of reads
        read_count_file = os.path.join(
            outdir, f"{sample}_TelomerCnt_{pid}", f"{pid}_readcount.tsv"
        )
        try:
            readcount = pd.read_table(read_count_file, sep="\t", comment="#")
        except FileNotFoundError:
            print(f"Didnt find {read_count_file}")
            continue

        # Check if the file has no rows, if so, continue to the next iteration
        if readcount.empty:
            print(f"Empty: {read_count_file}")
            continue

        total_reads = readcount["reads"].sum()

        # Normalize the values in the spectrum data (excluding the first three columns) are multiplied by the ratio of reads_with_pattern to the sum of each row.
        # Then, the values are further multiplied by 1e9 to scale them for normalization by total_reads.

        # Normalize based on reads_with_pattern
        # First ensure reads_with_pattern is float64
        spectrum["reads_with_pattern"] = spectrum["reads_with_pattern"].astype(
            "float64"
        )
        # Convert the numeric columns explicitly
        numeric_cols = spectrum.columns[3:]  # Get names of columns to convert
        spectrum[numeric_cols] = spectrum[numeric_cols].astype("float64")
        # Calculate row sums and normalization factors
        row_sums = spectrum[numeric_cols].sum(axis=1)
        normalization_factors = (spectrum["reads_with_pattern"] / row_sums).astype(
            "float64"
        )
        # Perform the multiplication
        for col in numeric_cols:
            spectrum[col] = spectrum[col].multiply(normalization_factors, axis=0)

        # Further normalize by a constant and total_reads
        spectrum.iloc[:, 3:] = spectrum.iloc[:, 3:] * 1000000000 / total_reads

        # Fill nan
        spectrum.fillna(0, inplace=True)

        spectrum_list[sample] = spectrum

    # Iterate over chromosomes
    for chr in chrs:
        spectrum_chr_list = {}

        for sample in samples:
            spectrum = spectrum_list[sample]
            spectrum_chr = spectrum[spectrum["chr"] == chr].copy()

            # Normalize by band length
            band_info_chr = band_info[band_info["chr"] == chr].copy()

            spectrum_chr.iloc[1:-1, 3:] *= 1e6 / band_info_chr["length"].values[:, None]

            spectrum_chr.set_index("band", inplace=True, drop=False)

            spectrum_chr_list[sample] = spectrum_chr

        bands = list(spectrum_chr_list[sample]["band"])

        if len(samples) == 2:
            height_matrix = combine_height_df(
                spectrum_chr_list["tumor"].iloc[:, 3:].T,
                spectrum_chr_list["control"].iloc[:, 3:].T,
            )
            sub_title = "Left: Tumor, Right: Control"
            axis = True
        else:
            height_matrix = spectrum_chr_list[sample].iloc[:, 3:].T
            sub_title = f"{samples[0]} Sample"
            axis = False

        if not height_matrix.empty and not height_matrix.isnull().values.all():
            main = f"{pid}: Telomere Repeat Types in Chr{chr}"

            barplot_repeattype(
                height_matrix=height_matrix,
                outdir=outdir,
                plot_name=f"{pid}_{chr}",
                plot_file_format=plot_file_format,
                main=main,
                sub_title=sub_title,
                ylab="Normalized Number of Telomere Reads",
                plot_reverse_complement=plot_rev_complement,
                repeat_threshold=repeat_threshold_plot,
                count_type="Consecutive" if consecutive else "Non-consecutive",
                mapq_threshold=mapq_threshold,
                control_axis=axis,
                labels=bands,
            )
        else:
            print(f"Height matrix is empty for chromosome {chr}. Skipping plotting.")


def plot_spectrum_summary(
    outdir,
    pid,
    repeat_threshold_plot,
    consecutive,
    mapq_threshold,
    banding_file,
    plot_rev_complement,
    plot_file_format,
):
    # Tumor, Control or both
    samples = get_sample_types_from_spectrum(outdir, pid)

    # Get band lengths
    band_info, chrs = get_band_info(banding_file)

    spectrum_summary_list = []

    for sample in samples:
        read_count_file = os.path.join(
            outdir, f"{sample}_TelomerCnt_{pid}", f"{pid}_readcount.tsv"
        )
        spectrum_file = os.path.join(
            outdir, f"{sample}_TelomerCnt_{pid}", f"{pid}_spectrum.tsv"
        )

        # Get spectrum
        spectrum = pd.read_table(spectrum_file, comment="#")

        # Make empty spectrum summary table
        spectrum_summary = pd.DataFrame(
            0,
            index=[
                "intra_chromosomal",
                "subtelomeric",
                "junction_spanning",
                "intra_telomeric",
            ],
            columns=spectrum.columns[2:],
        )

        # Get intra-telomeric reads (= unmapped reads)
        spectrum_summary.loc["intra_telomeric", :] = spectrum.loc[
            spectrum["chr"] == "unmapped", spectrum.columns[2:]
        ].values

        for chr_value in chrs:
            spectrum_chr = spectrum.loc[
                spectrum["chr"] == chr_value, spectrum.columns[2:]
            ]

            # Get junction spanning reads
            spectrum_summary.loc["junction_spanning"] += (
                spectrum_chr.iloc[0, :] + spectrum_chr.iloc[-1, :]
            )

            # Get subtelomeric reads (= first and last band)
            spectrum_summary.loc["subtelomeric"] += (
                spectrum_chr.iloc[1, :] + spectrum_chr.iloc[-2, :]
            )

            # Get intra-chromosomal reads (= all other bands)
            spectrum_summary.loc["intra_chromosomal"] += spectrum_chr.iloc[2:-2, :].sum(
                axis=0
            )

        # Get total number of reads
        readcount = pd.read_table(read_count_file, comment="#")
        total_reads = readcount["reads"].sum()

        # Normalize
        ratio = (
            spectrum_summary["reads_with_pattern"]
            / spectrum_summary.iloc[:, 1:].apply(sum, axis=1)
        ).fillna(0)
        spectrum_summary_norm = spectrum_summary.iloc[:, 1:].mul(ratio, axis=0)

        spectrum_summary_norm = spectrum_summary_norm * (1000000 / total_reads)

        # Bring into correct format
        spectrum_summary_norm["region"] = spectrum_summary_norm.index
        spectrum_summary_norm_m = pd.melt(
            spectrum_summary_norm,
            id_vars="region",
            var_name="repeat_type",
            value_name="sum",
        )
        spectrum_summary_norm_m["sample"] = sample

        spectrum_summary_list.append(spectrum_summary_norm_m)

    spectrum_summary = pd.concat(spectrum_summary_list)

    spectrum_summary = combine_forward_patterns(
        plot_rev_complement,
        spectrum_summary,
        ["sample", "region", "repeat_type_forward"],
    )

    # Change ordering of samples and region
    spectrum_summary["sample"] = pd.Categorical(
        spectrum_summary["sample"], categories=["tumor", "control"], ordered=True
    )
    spectrum_summary["region"] = pd.Categorical(
        spectrum_summary["region"],
        categories=[
            "intra_chromosomal",
            "subtelomeric",
            "junction_spanning",
            "intra_telomeric",
        ],
        ordered=True,
    )
    # Using cat.rename_categories to omit future warning
    spectrum_summary["region"] = spectrum_summary["region"].cat.rename_categories(
        {
            "intra_chromosomal": "intrachromosomal",
            "junction_spanning": "junction spanning",
            "intra_telomeric": "intratelomeric",
        }
    )

    # Change ordering of patterns and add colors
    colors = TVR_COLORS.copy()
    spectrum_summary, colors = update_repeat_types_order_and_col(
        spectrum_summary, colors
    )

    # Title
    main = (
        f"{pid}: Repeat types in telomeric reads"
        if len(pid) < 11
        else f"{pid}:\nRepeat types in telomeric reads"
    )

    # Replace NaN with 0
    spectrum_summary["sum"] = spectrum_summary["sum"].fillna(0)

    # Plot
    fig = px.bar(
        spectrum_summary,
        x="sample",
        y="sum",
        color="repeat_type",
        color_discrete_map=colors,
        facet_col="region",
        labels={"sum": "Telomeric reads (per million reads)"},
        title=main,
        category_orders={"sample": ["tumor", "control"]},
    )

    fig.update_layout(barmode="stack", showlegend=True)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))
    fig.update_layout(
        legend=dict(
            title=f'{repeat_threshold_plot} {"Consecutive" if consecutive else "Non-consecutive"} Repeats<br>Mapq Threshold = {mapq_threshold}',
            bgcolor="rgba(255,255,255,0.5)",
            bordercolor="rgba(0,0,0,0.5)",
            borderwidth=1,
            traceorder="reversed",
        ),
        **PLOT_CONFIG,
    )
    # fig.show()

    # Save plots
    save_plot(fig, outdir, f"{pid}_sorted_telomere_read_counts", plot_file_format)


def plot_tel_content(
    outdir,
    pid,
    repeat_threshold_plot,
    consecutive,
    mapq_threshold,
    plot_rev_complement,
    plot_file_format,
    gc_lower,
    gc_upper,
):
    # Tumor, control or both
    samples = get_sample_types_from_spectrum(outdir, pid)

    heights_samples_list = []

    for sample in samples:
        # Get spectrum
        spectrum_file = os.path.join(
            outdir, f"{sample}_TelomerCnt_{pid}", f"{pid}_spectrum.tsv"
        )
        spectrum = pd.read_table(spectrum_file, comment="#")
        spectrum = spectrum[spectrum["chr"] == "unmapped"]

        if spectrum.empty:
            print("No unmapped reads found -- skipping plotting tel_content")
            return

        # Get gc corrected telomere content
        summary_file = os.path.join(
            outdir, f"{sample}_TelomerCnt_{pid}", f"{pid}_{sample}_summary.tsv"
        )
        summary = pd.read_table(summary_file, sep="\t")
        tel_content = summary["tel_content"].values[0]

        # Get relative telomere repeat type occurrences in telomere content
        # Convert the numeric columns explicitly
        spectrum["reads_with_pattern"] = spectrum["reads_with_pattern"].astype(
            "float64"
        )
        numeric_cols = spectrum.columns[3:]  # Get names of columns to convert
        spectrum[numeric_cols] = spectrum[numeric_cols].astype("float64")
        # Calculate row sums
        row_sums = spectrum[numeric_cols].sum(axis=1)
        # Normalize based on reads_with_pattern
        normalization_factors = (spectrum["reads_with_pattern"] / row_sums).astype(
            "float64"
        )
        # Perform the multiplication for each numeric column
        for col in numeric_cols:
            spectrum[col] = (
                spectrum[col].multiply(normalization_factors, axis=0)
                * tel_content
                / spectrum["reads_with_pattern"]
            )

        height = spectrum.iloc[:, 3:].T
        height.columns = ["value"]
        height["repeat_type"] = height.index
        height["sample"] = sample

        heights_samples_list.append(height)

    height = pd.concat(heights_samples_list)

    # Combine forward and reverse patterns
    height = combine_forward_patterns(
        plot_rev_complement, height, ["sample", "repeat_type_forward"]
    )

    # Change ordering of samples
    height["sample"] = pd.Categorical(
        height["sample"], categories=["tumor", "control"], ordered=True
    )

    colors = TVR_COLORS.copy()
    height, colors = update_repeat_types_order_and_col(height, colors)

    # Plot
    fig = px.bar(
        height,
        x="sample",
        y="sum",
        color="repeat_type",
        color_discrete_map=colors,
        labels={
            "sum": f"Telomere content (intratelomeric reads per\nmillion reads <br>with GC Content of {gc_lower}-{gc_upper}%)"
        },
        title=f"{pid}: GC corrected telomere content",
        category_orders={"sample": ["tumor", "control"]},
    )

    fig.update_layout(
        legend=dict(
            title=f'{repeat_threshold_plot} {"Consecutive" if consecutive else "Non-consecutive"} Repeats <br>Mapq Threshold = {mapq_threshold}',
            bgcolor="rgba(255,255,255,0.5)",
            bordercolor="rgba(0,0,0,0.5)",
            borderwidth=1,
            traceorder="reversed",
        ),
        **PLOT_CONFIG,
    )
    # fig.show()

    # Save plots
    save_plot(fig, outdir, f"{pid}_telomere_content", plot_file_format)


def plot_gc_content(outdir, pid, plot_file_format, gc_lower, gc_upper):
    samples = []
    df_all = pd.DataFrame({"gc_content_percent": range(101)})
    df_intratel = pd.DataFrame({"gc_content_percent": range(101)})

    for sample_type in ["tumor", "control"]:
        gc_content_file = os.path.join(
            outdir,
            f"{sample_type}_TelomerCnt_{pid}",
            f"{pid}_{sample_type}_gc_content.tsv",
        )
        gc_content_file_intratel = os.path.join(
            outdir,
            f"{sample_type}_TelomerCnt_{pid}",
            f"{pid}_intratelomeric_{sample_type}_gc_content.tsv",
        )

        if os.path.exists(gc_content_file):
            samples.append(sample_type)
            gc_content = pd.read_csv(gc_content_file, sep="\t")
            gc_content_intratel = pd.read_csv(gc_content_file_intratel, sep="\t")

            gc_content["fraction_of_reads"] = (
                gc_content["read_count"] / gc_content["read_count"].sum()
            )
            gc_content_intratel["fraction_of_reads"] = (
                gc_content_intratel["read_count"]
                / gc_content_intratel["read_count"].sum()
            )

            df_all[sample_type] = gc_content["fraction_of_reads"]
            df_intratel[sample_type] = gc_content_intratel["fraction_of_reads"]

    dfm_all = pd.melt(
        df_all,
        id_vars=["gc_content_percent"],
        var_name="sample",
        value_name="fraction_of_reads",
    )
    dfm_intratel = pd.melt(
        df_intratel,
        id_vars=["gc_content_percent"],
        var_name="sample",
        value_name="fraction_of_reads",
    )

    dfm_all["read_type"] = "All reads"
    dfm_intratel["read_type"] = "Intratelomeric reads"

    dfm = pd.concat([dfm_all, dfm_intratel])

    dfm["group"] = pd.Categorical(dfm["sample"], categories=["tumor", "control"])

    # Plot
    fig = px.line(
        dfm,
        x="gc_content_percent",
        y="fraction_of_reads",
        color="sample",
        line_group="group",
        facet_col="read_type",
        facet_col_wrap=2,
        labels={"fraction_of_reads": "Fraction of reads"},
    )

    fig.add_vrect(
        x0=gc_lower - 0.5,
        x1=gc_upper + 0.5,
        annotation_text="bin used for GC correction",
        annotation_position="top left",
        fillcolor="grey",
        opacity=0.5,
        line_width=0,
        layer="below",
    )

    # Remove individual x-axis titles from facets
    for axis in fig.layout:
        if axis.startswith("xaxis") and "title" in fig.layout[axis]:
            fig.layout[axis].title = "GC content [%]"

    # Set shared X axis title on bottom
    fig.update_layout(
        title_text=f"{pid}: GC content",
        legend_title_text="",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(t=60, b=10),
        showlegend=True,
        yaxis=dict(range=[0, dfm["fraction_of_reads"].max() + 0.1]),
        **PLOT_CONFIG,
    )

    # fig.show()

    # Save plots
    save_plot(fig, outdir, f"{pid}_gc_content", plot_file_format)


def plot_repeat_frequency_intratelomeric(
    outdir,
    pid,
    repeat_threshold_plot,
    consecutive,
    mapq_threshold,
    repeats,
    plot_file_format,
):
    tumor_file_path = os.path.join(
        outdir,
        f"tumor_TelomerCnt_{pid}",
        f"{pid}_repeat_frequency_per_intratelomeric_read.tsv",
    )
    control_file_path = os.path.join(
        outdir,
        f"control_TelomerCnt_{pid}",
        f"{pid}_repeat_frequency_per_intratelomeric_read.tsv",
    )

    # Process data and pass empty dataframe if no values
    frequency_table_T = read_frequency_table(tumor_file_path, "tumor")
    frequency_table_C = read_frequency_table(control_file_path, "control")

    # Concatenate dataframes
    df = pd.concat(
        [df for df in [frequency_table_T, frequency_table_C] if not df.empty],
        ignore_index=True,
    )
    # If both DataFrames are empty, create an empty DataFrame
    if df.empty:
        df = pd.DataFrame()

    df["sample"] = pd.Categorical(df["sample"], categories=["tumor", "control"])

    # Convert columns to numeric
    df["number_repeats"] = pd.to_numeric(df["number_repeats"])

    # Make histogram
    p_hist = px.bar(
        df,
        x="number_repeats",
        y="percent",
        color="sample",
        facet_col="sample",
        labels={
            "number_repeats": "Number of repeats",
            "percent": "Percent of intratelomeric reads",
        },
        title=f"{pid}: Frequency of telomere repeat occurrences in intratelomeric reads",
    )

    # p_hist.show()

    # Make cumulative line plot with two lines
    p_cum = px.line(
        df,
        x="number_repeats",
        y="percent_cumulative",
        color="sample",
        labels={
            "number_repeats": "Number of repeats",
            "percent_cumulative": "Percent (cumulative)",
        },
        title=f"{pid}: Frequency of telomere repeat occurrences in intratelomeric reads",
    )

    # Reverse the x-axis
    p_cum.update_xaxes(autorange="reversed")

    # p_cum.show()

    # Merge plots
    fig = make_subplots(
        rows=1, cols=2, subplot_titles=["Histogram", "Cumulative Line Plot"]
    )
    fig.add_trace(p_hist["data"][0], row=1, col=1)
    fig.add_trace(p_cum["data"][0], row=1, col=2)

    if len(frequency_table_C) > 0:
        fig.add_trace(p_hist["data"][1], row=1, col=1)
        fig.add_trace(p_cum["data"][1], row=1, col=2)

    # Set layout and save plot
    consecutive_txt = "consecutive" if consecutive else "non-consecutive"
    fig.update_layout(
        title=f"{pid}: Frequency of telomere repeats in intratelomeric reads <br>"
        f"<sup>Filtering criteria: {repeat_threshold_plot} {consecutive_txt} repeats, mapq threshold = {mapq_threshold}. Repeat types = {repeats}</sup>",
        # legend=dict(title='', orientation='h', y=1.1, x=0.5),
        xaxis=dict(title="Number of repeats"),
        yaxis=dict(title="Percent (cumulative)"),
        xaxis2=dict(title="Number of repeats"),
        yaxis2=dict(title="Percent of intratelomeric reads"),
        **PLOT_CONFIG,
    )
    fig.update_xaxes(autorange="reversed", row=1, col=2)

    # fig.show()

    # Save plot
    save_plot(
        fig,
        outdir,
        f"{pid}_hist_telomere_repeats_per_intratelomeric_read",
        plot_file_format,
    )


def plot_TVR_plot(outdir, pid, plot_file_format):
    # Read the table file
    table_file = os.path.join(outdir, f"{pid}_normalized_TVR_counts.tsv")
    table_merged = pd.read_table(table_file)

    # Check if the DataFrame is empty
    if table_merged.empty:
        print(f"Skipping plotting TVRs as the file {table_file} is empty.")
        return

    # Only keep patterns if the normalized pattern count is bigger than 0.01 in tumor or control
    table_merged_filter = table_merged[
        (table_merged["Count_norm_by_intratel_reads_T"] >= 0.01)
        | (table_merged["Count_norm_by_intratel_reads_C"] >= 0.01)
    ]

    # Create a list of patterns for color mapping
    patterns = table_merged_filter["Pattern"].unique()

    # Create the figure for absolute counts (main plot)
    fig_absolute_counts = go.Figure()

    # Add bars for Tumor samples
    for pattern in patterns:
        mask = table_merged_filter["Pattern"] == pattern
        fig_absolute_counts.add_trace(
            go.Bar(
                x=table_merged_filter[mask]["Pattern"],
                y=table_merged_filter[mask]["Count_norm_by_intratel_reads_T"],
                name="Tumor",
                marker_color=TVR_COLORS.get(pattern, "#000000"),
                offsetgroup=0,
                showlegend=False,  # Hide individual bars in the legend
            )
        )

    # Add bars for Control samples with 0.5 opacity
    for pattern in patterns:
        mask = table_merged_filter["Pattern"] == pattern
        fig_absolute_counts.add_trace(
            go.Bar(
                x=table_merged_filter[mask]["Pattern"],
                y=table_merged_filter[mask]["Count_norm_by_intratel_reads_C"],
                name="Control",
                marker_color=TVR_COLORS.get(pattern, "#000000"),
                marker_opacity=0.5,
                offsetgroup=1,
                showlegend=False,  # Hide individual bars in the legend
            )
        )

    # Add representative legend entries
    fig_absolute_counts.add_trace(
        go.Bar(
            x=["TTAGGG"],
            y=[0],  # Invisible bar for legend
            name="Tumor (Full Color)",
            marker_color="#212F3D",
            showlegend=True,
        )
    )

    fig_absolute_counts.add_trace(
        go.Bar(
            x=["TTAGGG"],
            y=[0],  # Invisible bar for legend
            name="Control (Lighter, Alpha=0.5)",
            marker_color="#ABB2B9",
            marker_opacity=0.5,
            showlegend=True,
        )
    )

    # Update layout
    fig_absolute_counts.update_layout(
        barmode="group",
        title=f"{pid}: TVRs found in intratelomeric reads",
        xaxis_title="Patterns",
        yaxis_title="Mean TVR counts per intratelomeric read",
        legend_title="Sample",
        **PLOT_CONFIG,
    )

    # Create inset plot without TTAGGG if TTAGGG is present
    if "TTAGGG" in patterns:
        # Filter data without TTAGGG
        table_no_ttaggg = table_merged_filter[
            table_merged_filter["Pattern"] != "TTAGGG"
        ]

        # Get patterns excluding TTAGGG
        patterns_no_ttaggg = [p for p in patterns if p != "TTAGGG"]

        # Add inset traces to the main figure with xaxis2 and yaxis2
        for pattern in patterns_no_ttaggg:
            mask = table_no_ttaggg["Pattern"] == pattern
            fig_absolute_counts.add_trace(
                go.Bar(
                    x=table_no_ttaggg[mask]["Pattern"],
                    y=table_no_ttaggg[mask]["Count_norm_by_intratel_reads_T"],
                    name="Tumor_inset",
                    marker_color=TVR_COLORS.get(pattern, "#000000"),
                    offsetgroup=0,
                    showlegend=False,
                    xaxis="x2",
                    yaxis="y2",
                )
            )

        # Add bars for Control samples (without TTAGGG) to inset
        for pattern in patterns_no_ttaggg:
            mask = table_no_ttaggg["Pattern"] == pattern
            fig_absolute_counts.add_trace(
                go.Bar(
                    x=table_no_ttaggg[mask]["Pattern"],
                    y=table_no_ttaggg[mask]["Count_norm_by_intratel_reads_C"],
                    name="Control_inset",
                    marker_color=TVR_COLORS.get(pattern, "#000000"),
                    marker_opacity=0.5,
                    offsetgroup=1,
                    showlegend=False,
                    xaxis="x2",
                    yaxis="y2",
                )
            )

        # Update layout to include the overlapping inset in upper right corner
        fig_absolute_counts.update_layout(
            barmode="group",
            title=f"{pid}: TVRs found in intratelomeric reads",
            xaxis_title="Patterns",
            yaxis_title="Mean TVR counts per intratelomeric read",
            legend_title="Sample",
            # Main plot axes
            xaxis=dict(domain=[0, 1]),
            yaxis=dict(domain=[0, 1]),
            # Inset plot axes (positioned in upper right corner)
            xaxis2=dict(
                domain=[0.6, 0.95],  # Bigger inset - increased width
                anchor="y2",
                title="",  # Remove axis title
                showticklabels=True,  # Show x-axis tick labels
                tickfont=dict(size=8),
            ),
            yaxis2=dict(
                domain=[0.6, 0.95],  # Bigger inset - increased height
                anchor="x2",
                title="",  # Remove axis title
                showticklabels=True,  # Show y-axis tick labels
                tickfont=dict(size=8),
            ),
            # Add background for inset to make it stand out
            shapes=[
                dict(
                    type="rect",
                    xref="paper",
                    yref="paper",
                    x0=0.57,  # Adjusted for bigger rectangle
                    y0=0.48,  # Adjusted for bigger rectangle
                    x1=1.00,  # Adjusted for bigger rectangle
                    y1=1.00,  # Adjusted for bigger rectangle
                    fillcolor="white",
                    opacity=0.8,
                    line=dict(color="black", width=1),
                    layer="below",  # Put rectangle behind the plot
                )
            ],
            # Add annotation for inset title in upper right corner
            annotations=[
                dict(
                    text="TVRs without TTAGGG",
                    xref="paper",
                    yref="paper",
                    x=0.95,  # Positioned in upper right corner
                    y=0.95,  # Positioned in upper right corner
                    showarrow=False,
                    font=dict(size=10, color="black"),
                    xanchor="right",  # Anchor text to the right
                    yanchor="top",  # Anchor text to the top
                )
            ],
            **PLOT_CONFIG,
        )

        # Save the combined plot with overlapping inset
        save_plot(
            fig_absolute_counts,
            outdir,
            f"{pid}_TVR_plots_absolute_counts",
            plot_file_format,
        )
    else:
        # If TTAGGG is not present, just save the original plot
        save_plot(
            fig_absolute_counts,
            outdir,
            f"{pid}_TVR_plots_absolute_counts",
            plot_file_format,
        )

    # Create bar plot with log2 ratio, only if control
    if table_merged_filter["log2_ratio_count_norm_by_intratel_reads"].any():
        condition_top = (
            table_merged_filter["log2_ratio_count_norm_by_intratel_reads"] > 5
        )
        condition_bottom = (
            table_merged_filter["log2_ratio_count_norm_by_intratel_reads"] < -5
        )

        # Create an explicit copy
        table_merged_filter_copy = table_merged_filter.copy()

        # Perform operations on the copy
        table_merged_filter_copy.loc[condition_top, "label_top"] = (
            table_merged_filter_copy.loc[
                condition_top, "log2_ratio_count_norm_by_intratel_reads"
            ].round(1)
        )

        table_merged_filter_copy.loc[condition_bottom, "label_bottom"] = (
            table_merged_filter_copy.loc[
                condition_bottom, "log2_ratio_count_norm_by_intratel_reads"
            ].round(1)
        )

        # Create the log2 ratio figure with individual bars to ensure consistent colors
        fig_log2_ratio = go.Figure()

        # Add bars with consistent colors
        for pattern in patterns:
            mask = table_merged_filter_copy["Pattern"] == pattern
            fig_log2_ratio.add_trace(
                go.Bar(
                    x=table_merged_filter_copy[mask]["Pattern"],
                    y=table_merged_filter_copy[mask][
                        "log2_ratio_count_norm_by_intratel_reads"
                    ],
                    marker_color=TVR_COLORS.get(pattern, "#000000"),
                    showlegend=False,  # Remove legend as requested
                )
            )

        # Update layout
        fig_log2_ratio.update_layout(
            title="Normalized TVR counts tumor/control (log2)",
            xaxis_title="Patterns",
            yaxis_title="log2 Ratio",
            yaxis=dict(range=[-5, 5]),  # Fix y-axis range as requested
            showlegend=False,  # Remove legend as requested
            **PLOT_CONFIG,
        )

        # Add labels for values greater than 5 or less than -5
        for index, row in table_merged_filter_copy.iterrows():
            if row["log2_ratio_count_norm_by_intratel_reads"] > 5:
                fig_log2_ratio.add_annotation(
                    x=row["Pattern"],
                    y=5,
                    text=str(round(row["log2_ratio_count_norm_by_intratel_reads"], 1)),
                    showarrow=True,
                    arrowhead=1,
                )
            elif row["log2_ratio_count_norm_by_intratel_reads"] < -5:
                fig_log2_ratio.add_annotation(
                    x=row["Pattern"],
                    y=-5,
                    text=str(round(row["log2_ratio_count_norm_by_intratel_reads"], 1)),
                    showarrow=True,
                    arrowhead=1,
                )

        # Save plot
        save_plot(
            fig_log2_ratio, outdir, f"{pid}_TVR_plots_log2_ratio", plot_file_format
        )

    # Additional processing for scatter plot if both T and C are present
    table_for_scatter = table_merged.copy()
    if "Pattern" in table_for_scatter.columns:
        table_for_scatter.set_index("Pattern", inplace=True)

    if (
        "TTAGGG" in table_for_scatter.index
        and not pd.isnull(
            table_for_scatter.loc["TTAGGG", "Count_norm_by_intratel_reads_T"]
        )
        and not pd.isnull(
            table_for_scatter.loc["TTAGGG", "Count_norm_by_intratel_reads_C"]
        )
    ):

        # Reset index for easier use with plotly
        table_for_scatter = table_for_scatter.reset_index()

        table_for_scatter["label"] = "other"
        label_patterns = ["TTAGGG", "TGAGGG", "TCAGGG", "TTGGGG", "TTCGGG"]
        for pattern in label_patterns:
            if pattern in table_for_scatter["Pattern"].values:
                table_for_scatter.loc[
                    table_for_scatter["Pattern"] == pattern, "label"
                ] = pattern

        # Create scatter plot with consistent colors
        fig_scatter = go.Figure()

        # Add points with consistent pattern colors
        for pattern in table_for_scatter["Pattern"].unique():
            mask = table_for_scatter["Pattern"] == pattern
            fig_scatter.add_trace(
                go.Scatter(
                    x=table_for_scatter[mask]["Count_norm_by_intratel_reads_C"],
                    y=table_for_scatter[mask]["Count_norm_by_intratel_reads_T"],
                    mode="markers",
                    name=pattern,
                    marker=dict(color=TVR_COLORS.get(pattern, "#000000"), size=10),
                )
            )

        # Determine range for diagonal line
        x_min = table_for_scatter["Count_norm_by_intratel_reads_C"].min() * 0.9
        x_max = table_for_scatter["Count_norm_by_intratel_reads_C"].max() * 1.1

        # Add diagonal line
        fig_scatter.add_trace(
            go.Scatter(
                x=[x_min, x_max],
                y=[x_min, x_max],
                mode="lines",
                line=dict(dash="dash", color="black"),
                name="y = x",
            )
        )

        # Update layout
        fig_scatter.update_layout(
            title=f"{pid}: TVRs per intratelomeric read",
            xaxis_title="Mean TVR counts per intratel. read (control)",
            yaxis_title="Mean TVR counts per intratel. read (tumor)",
            xaxis_type="log",
            yaxis_type="log",
            **PLOT_CONFIG,
        )

        # Save plot
        save_plot(
            fig_scatter, outdir, f"{pid}_pattern_screen_scatterplot", plot_file_format
        )


def plot_singletons(outdir, pid, plot_file_format):
    """Generate singleton analysis plots with improved visualization.

    Args:
        outdir: Output directory path
        pid: Patient/Sample ID
        plot_file_format: List of output formats (e.g., ['pdf', 'png'])
    """
    try:
        # Read and validate data
        singleton_file = os.path.join(outdir, f"{pid}_singletons.tsv")
        singletons = pd.read_csv(singleton_file, sep="\t")

        required_cols = ["pattern", "singleton_count_tumor", "singleton_count_control"]
        if not validate_data(singletons, required_cols):
            print("Invalid or empty singleton data.")
            return

        # Create main figure with subplots
        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=[
                "Raw Singleton Counts",
                "Singleton Counts Normalized by total reads count",
                "Normalized Log2 Ratio (Tumor/Control)",
                "Distance to Expected Singleton Count",
            ],
            vertical_spacing=0.2,
            horizontal_spacing=0.15,
        )

        # Add individual plots
        add_raw_counts_subplot(fig, singletons, 1, 1)
        add_normalized_counts_subplot(fig, singletons, 1, 2)
        add_log2_ratio_subplot(fig, singletons, 2, 1)
        add_distance_subplot(fig, singletons, 2, 2)

        # Add representative legend entries
        fig.add_trace(
            go.Bar(
                x=[
                    singletons.pattern.iloc[0]
                ],  # grab first singleton entry to show fake legend
                y=[0],
                name="Tumor (Full Color)",
                marker_color="#212F3D",
                showlegend=True,
            )
        )

        fig.add_trace(
            go.Bar(
                x=[
                    singletons.pattern.iloc[0]
                ],  # grab first singleton entry to show fake legend
                y=[0],
                name="Control (Lighter, Alpha=0.5)",
                marker_color="#ABB2B9",
                marker_opacity=0.5,
                showlegend=True,
            )
        )

        # Update layout
        fig.update_layout(
            title=dict(
                text=f"{pid}: Singleton Analysis<br><sup>(TTAGGG)3-NNNGGG-(TTAGGG)3</sup>",
                x=0.5,
                y=0.95,
            ),
            # margin=dict(b=50),  # Moved margin to main layout level
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                xanchor="center",
                x=0.5,
                y=-0.2,  # Moved legend lower
            ),
            **PLOT_CONFIG,
        )
        # fig.show()

        # Save plot
        save_plot(fig, outdir, f"{pid}_singletons", plot_file_format, height=600)

    except Exception as e:
        print(f"Error generating singleton plots: {str(e)}")


def prepare_singleton_data(df):
    """Prepare singleton data for plotting."""
    # Sort by tumor count and set up categorical patterns
    df = df.sort_values("singleton_count_tumor", ascending=False)
    df["pattern"] = pd.Categorical(df["pattern"], categories=df["pattern"].unique())

    # Calculate normalized values (percentages)
    for sample in ["tumor", "control"]:
        col = f"singleton_count_{sample}"
        if col in df.columns:
            total = df[col].sum()
            if total > 0:
                df[f"{col}_norm"] = df[col] / total * 100

    # Calculate log2 ratio and distances
    if (
        "singleton_count_tumor" in df.columns
        and "singleton_count_control" in df.columns
    ):
        df["log2_ratio"] = np.log2(
            df["singleton_count_tumor_norm"]
            / df["singleton_count_control_norm"].replace([0.0, 0], np.nan)
        )

        if "tel_content_log2_ratio" in df.columns:
            df["distance"] = df["log2_ratio"] - df["tel_content_log2_ratio"]

    return df


def add_raw_counts_subplot(fig, df, row, col):
    """Add raw count bars subplot."""
    for sample in ["tumor", "control"]:
        col_name = f"singleton_count_{sample}"
        if col_name in df.columns:
            opacity = 0.5 if sample == "control" else 1.0
            fig.add_trace(
                go.Bar(
                    x=df["pattern"],
                    y=df[col_name],
                    name=f"{sample.capitalize()}",
                    marker_color=[
                        TVR_COLORS.get(p, TVR_COLORS["other"]) for p in df["pattern"]
                    ],
                    marker_opacity=opacity,
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

    fig.update_yaxes(title_text="Count", type="log", row=row, col=col)
    fig.update_xaxes(tickangle=45, row=row, col=col)


def add_normalized_counts_subplot(fig, df, row, col):
    """Add normalized percentage subplot."""
    for sample in ["tumor", "control"]:
        col_name = f"singleton_count_{sample}_norm"
        if col_name in df.columns:
            opacity = 0.5 if sample == "control" else 1.0
            fig.add_trace(
                go.Bar(
                    x=df["pattern"],
                    y=df[col_name],
                    name=f"{sample.capitalize()}",
                    marker_color=[
                        TVR_COLORS.get(p, TVR_COLORS["other"]) for p in df["pattern"]
                    ],
                    marker_opacity=opacity,
                    showlegend=False,
                ),
                row=row,
                col=col,
            )

    fig.update_yaxes(title_text="Normalized singleton counts", row=row, col=col)
    fig.update_xaxes(tickangle=45, row=row, col=col)


def add_log2_ratio_subplot(fig, df, row, col):
    """Add log2 ratio subplot with reference line."""
    if "singleton_count_log2_ratio_norm" in df.columns:
        fig.add_trace(
            go.Bar(
                x=df["pattern"],
                y=df["singleton_count_log2_ratio_norm"],
                marker_color=[
                    TVR_COLORS.get(p, TVR_COLORS["other"]) for p in df["pattern"]
                ],
                showlegend=False,
            ),
            row=row,
            col=col,
        )

        # Add reference line if available
        if "tel_content_log2_ratio" in df.columns:
            fig.add_hline(
                y=df["tel_content_log2_ratio"].iloc[0],
                line_dash="dash",
                line_color="gray",
                row=row,
                col=col,
                annotation_text="Expected ratio",
            )

    fig.update_yaxes(
        title_text="Normalized singleton count <br> tumor/control (log2)",
        zeroline=True,
        range=[-5, 5],
        row=row,
        col=col,
    )
    fig.update_xaxes(tickangle=45, row=row, col=col)


def add_distance_subplot(fig, df, row, col):
    """Add distance from expected subplot."""
    if "distance_to_expected_singleton_log2_ratio" in df.columns:
        fig.add_trace(
            go.Bar(
                x=df["pattern"],
                y=df["distance_to_expected_singleton_log2_ratio"],
                marker_color=[
                    TVR_COLORS.get(p, TVR_COLORS["other"]) for p in df["pattern"]
                ],
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    fig.update_yaxes(
        title_text="Distance to expected <br> singleton count",
        zeroline=True,
        range=[-5, 5],
        row=row,
        col=col,
    )
    fig.update_xaxes(tickangle=45, row=row, col=col)
