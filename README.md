# TelomereHunter2

[![PyPI version](https://img.shields.io/pypi/v/telomerehunter2.svg)](https://pypi.org/project/telomerehunter2/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Build Status](https://img.shields.io/github/actions/workflow/status/ferdinand-popp/telomerehunter2/pypi-release.yml?branch=main)](https://github.com/ferdinand-popp/telomerehunter2/actions)
[![Python Versions](https://img.shields.io/pypi/pyversions/telomerehunter2.svg)](https://pypi.org/project/telomerehunter2/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
[![Last Commit](https://img.shields.io/github/last-commit/ferdinand-popp/telomerehunter2.svg)](https://github.com/ferdinand-popp/telomerehunter2/commits/main)
[![Docker Pulls](https://img.shields.io/docker/pulls/fpopp22/telomerehunter2)](https://hub.docker.com/r/fpopp22/telomerehunter2)

TelomereHunter2 is a Python-based tool for estimating telomere content and analyzing telomeric variant repeats (TVRs)
from genome sequencing data. It supports BAM/CRAM files, flexible telomere repeat and reference genome inputs, and
provides outputs for bulk and single-cell genome sequencing data.

---

## Release Notes

See [RELEASE_NOTES.md](RELEASE_NOTES.md) for the latest changes and version history.

---

## New Features

- Fast, container-friendly Python 3 implementation
- Parallelization and algorithmic steps for drastic speedup
- Supports BAM/CRAM, custom telomeric repeats, and now also non-human genomes
- Static and interactive HTML reports (Plotly)
- Docker and Apptainer/Singularity containers
- Single cell sequencing support (e.g. scATAC-seq; barcode splitting and per-cell analysis)
- Robust input handling and exception management
- Fast mode for quick overview of unmapped reads

## Installation

**Classic setup:**

```bash
pip install telomerehunter2
```

**From source:**

```bash
# With pip:
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2
python -m venv venv
source venv/bin/activate
pip install -e . --no-cache-dir

# With uv:
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2
uv pip install -e . --no-cache-dir
```

**Container usage:**  
See [Container Usage](#container-usage) for Docker/Apptainer instructions.

**Operating systems:**  
Currently tested on Linux and macOS. Windows support via WSL2 and Docker not completely tested (WIP check GitHub Issues)

## Usage Bulk vs single cell Analysis

### Bulk Analysis

```bash   
telomerehunter2 -ibt TUMOR_FILE -ibc CONTROL_FILE -o OUTPUT_DIRECTORY -p ID_OF_SAMPLE -b BANDING_FILE [options]
```

- **Single sample:**  
  `telomerehunter2 -ibt sample.bam -o results/ -p SampleID -b telomerehunter2/cytoband_files/hg19_cytoBand.txt`
- **Tumor vs Control:**  
  `telomerehunter2 -ibt sample.bam -ibc control.bam -o results/ -p PairID -b telomerehunter2/cytoband_files/hg19_cytoBand.txt`
- **Custom repeats/species:**  
  `telomerehunter2 ... --repeats TTTAGGG TTAAGGG --repeatsContext TTAAGGG`
- **Fast mode (quick overview of unmapped reads generating summary with overview):**  
  `telomerehunter2 -ibt sample.bam -o results/ -p SampleID --fast_mode`

### Single cell sequencing Analysis

TelomereHunter2 now supports direct single-cell BAM analysis (with CB barcode tag). Simply run:

`telomerehunter2_sc -ibt sample.bam -o results/ -p SampleID -b telomerehunter2/cytoband_files/cytoband.txt --min-reads-per-barcode 10000`

This will perform barcode-aware telomere analysis and output per-cell results in a summary file. The minimum reads per
barcode threshold can be set with `--min-reads-per-barcode`. To rerun postprocessing with adjusted `--min-reads-per-barcode` 
threshold run command again with `--noFiltering` to skip the expensive filtering step from all reads to telomeric reads.
If the reads have a different barcode tag than `CB`, use `--barcodeTag` to set the correct one.
More information on correcting chromatin state for scATAC follows in (Engel et al., 2024).

See `tests/test_telomerehunter2_sc.py` for example usage and validation.

### Usage full list of option

`telomerehunter2 --help`

## Input & Output

**Input:**

- BAM/CRAM files (aligned reads, <-ibt> for tumor, <-ibc> for control)
- Cytoband file (tab-delimited, e.g. `telomerehunter2/cytoband_files/hg19_cytoBand.txt`, <-b>)
- Identifier for sample/pair (<-p>)
- Optional: custom telomeric repeats

**Output:**

- `summary.tsv`, `TVR_top_contexts.tsv`, `singletons.tsv`
- Plots (`plots/` directory, PNG/HTML)
- Logs (run status/errors)
- For sc-seq: Additionally to the complete bulk run you get per-cell results in sc_summary.tsv and barcode_counts.tsv
  with reads counts per barcode

### Explanation of summary.tsv

| Column                                          | Value example | Description                                                         |
|-------------------------------------------------|---------------|---------------------------------------------------------------------|
| PID                                             | TEST_PATIENT  | Sample name                                                         |
| sample                                          | tumor         | Sample classification (tumor (single), control, log2(t/c))          |
| **tel_content**                                 | 1.8           | Intratelomeric reads / reads in GC correction range * 1e6           |
| total_reads                                     | 120           | Number of reads in the input file                                   |
| read_lengths                                    | 25,36,42,54   | Unique lengths of reads                                             |
| repeat_threshold_set                            | 6 per 100 bp  | Telomeric repeat threshold set                                      |
| repeat_threshold_used                           | 4             | Repeats threshold applied based on avg. read length                 |
| intratelomeric_reads                            | 4             | Filtered Tel reads in unmapped reads                                |
| junctionspanning_reads                          | 0             | Filtered Tel reads spanning junctions into first/last band          |
| subtelomeric_reads                              | 6             | Filtered Tel reads in subtelomeric regions (first/last band)        |
| intrachromosomal_reads                          | 0             | Filtered Tel reads in intrachromosomal regions                      |
| tel_read_count                                  | 10            | Total telomeric reads identified                                    |
| gc_bins_for_correction                          | 48-52         | GC content range used for normalization of reads                    |
| total_reads_with_tel_gc                         | 8             | Total reads within GC bin for normalization                         |
| TCAGGG_arbitrary_context_norm_by_intratel_reads | 1.5           | Telomeric variant repeat count normalized by intratelomeric reads   |
| ...                                             | ...           | ...                                                                 |
| TCAGGG_singletons_norm_by_all_reads             | 0.0           | Singleton (TVR flanked by canonicals) count normalized by all reads |
| ...                                             | ...           | ...                                                                 |

## Dependencies

- Python >=3.6
- pysam, numpy, pandas, plotly, PyPDF2
- For static image export: kaleido (requires chrome/chromium)
- Docker/Apptainer (optional)

Install all dependencies:

```bash
pip install -r requirements.txt
```

## Container Usage

**Docker (recommended):**

*Build locally:*

```bash
docker build -t telomerehunter2 .
docker run --rm -it -v /data:/data telomerehunter2 telomerehunter2 -ibt /data/sample.bam -o /data/results -p SampleID -b /data/hg19_cytoBand.txt
```

*Pull from Docker Hub:*

```bash
docker pull fpopp22/telomerehunter2
```

*Run from Docker Hub:*

```bash
docker run --rm -it -v /data:/data fpopp22/telomerehunter2 telomerehunter2 -ibt /data/sample.bam -o /data/results -p SampleID -b /data/hg19_cytoBand.txt
```

**Apptainer/Singularity:**

*Build locally:*

```bash
apptainer build telomerehunter2.sif Apptainer_TH2.def
# mount data needed
apptainer run telomerehunter2.sif telomerehunter2 -ibt /data/sample.bam -o /data/results -p SampleID -b /data/hg19_cytoBand.txt
```

*Pull from Docker Hub (as Apptainer image):*

```bash
apptainer pull docker://fpopp22/telomerehunter2:latest
apptainer run telomerehunter2_latest.sif telomerehunter2 ...
```

## Troubleshooting

- **Memory errors:** Use more RAM or limit cores used with `-c` flag.
- **Missing dependencies:** Check `requirements.txt`.
- **Banding file missing:** Needs reference genome banding file `-b` otherwise analysis will run without reads mapped to
  subtelomeres.
- **Plotting:** Try disabling with `--plotNone` or use plotting only mode with `--plotNone`.
- **Minor changes to TH1:** Skipping the tvrs normalization per 100 bp, improved detection of GXXGGG TVRs, read lengths
  are estimated from first 1000 reads, added TRPM

For help: [GitHub Issues](https://github.com/fpopp22/telomerehunter2/issues) or our FAQ.

## Documentation & Resources

- [Wiki](https://github.com/fpopp22/telomerehunter2/wiki) (WIP)
- [Example Data](tests/)
- [Tutorial Videos](https://github.com/fpopp22/telomerehunter2/wiki) (WIP)
- [Telomerehunter Website](https://www.dkfz.de/angewandte-bioinformatik/telomerehunter)
- [Original TelomereHunter Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0)

## Citation

If you use TelomereHunter2, please cite:

- Feuerbach, L., et al. "TelomereHunter – in silico estimation of telomere content and composition from cancer genomes."
  BMC Bioinformatics 20, 272 (2019). https://doi.org/10.1186/s12859-019-2851-0
- Application Note for TH2 (in preparation).

## Contributing

Fork, branch, and submit pull requests. Please add tests and follow code style. For major changes, open an issue first.
Before submitting, please install the `tox` package and run the following checks:

1. **Run Unit Tests and Style Checks**:
   ```bash
   tox
   ```

## License

GNU General Public License v3.0. See [LICENSE](LICENSE.txt).

## Contact

- Ferdinand Popp (f.popp@dkfz.de)
- Lars Feuerbach (l.feuerbach@dkfz.de)

## Acknowledgements

Developed by Ferdinand Popp, Lina Sieverling, Philip Ginsbach, Lars Feuerbach. Supported by German Cancer Research
Center (DKFZ) - Division Applied Bioinformatics.

---

Copyright 2025 Ferdinand Popp, Lina Sieverling, Philip Ginsbach, Lars Feuerbach
