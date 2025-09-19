# TelomereHunter2

[![PyPI version](https://img.shields.io/pypi/v/telomerehunter2.svg)](https://pypi.org/project/telomerehunter2/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Build Status](https://img.shields.io/github/actions/workflow/status/ferdinand-popp/telomerehunter2/pypi-release.yml?branch=main)](https://github.com/ferdinand-popp/telomerehunter2/actions)
[![Python Versions](https://img.shields.io/pypi/pyversions/telomerehunter2.svg)](https://pypi.org/project/telomerehunter2/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
[![Last Commit](https://img.shields.io/github/last-commit/ferdinand-popp/telomerehunter2.svg)](https://github.com/ferdinand-popp/telomerehunter2/commits/main)
[![Docker Pulls](https://img.shields.io/docker/pulls/fpopp22/telomerehunter2)](https://hub.docker.com/r/fpopp22/telomerehunter2)

TelomereHunter2 is a Python-based tool for estimating telomere content and analyzing telomeric variant repeats (TVRs) from genome sequencing data. It supports BAM/CRAM files, flexible telomere repeat and reference genome inputs, and provides outputs for bulk and single-cell genome sequencing data.

## Features

- Fast, container-friendly Python 3 implementation
- Parallelization and algorithmic steps for drastic speedup
- Supports BAM/CRAM, custom telomeric repeats, and now also non-human genomes
- Static and interactive HTML reports (Plotly)
- Docker and Apptainer/Singularity containers
- Single cell sequencing support (e.g. scATAC-seq; barcode splitting and per-cell analysis)
- Robust input handling and exception management
- **Fast mode for quick overview of unmapped reads**

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

# With poetry:
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2
poetry env use python3
poetry install

# With uv:
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2
uv pip install -e . --no-cache-dir
```

**Container usage:**  
See [Container Usage](#container-usage) for Docker/Apptainer instructions.

## Quickstart

### Bulk Analysis

```bash
telomerehunter2 -ibt sample.bam -o results/ -p SampleID -b telomerehunter2/cytoband_files/hg19_cytoBand.txt
```
For all options:
```bash
telomerehunter2 --help
```

### Single Cell Analysis

```bash
telomerehunter2_sc -ibt sample.bam -o results/ -p SampleID -b telomerehunter2/cytoband_files/cytoband.txt --min-reads-per-barcode 10000
```

See [Bulk Analysis](#bulk-analysis) and [Single cell sequencing Analysis](#single-cell-sequencing-analysis) below for more details.

## Usage

### Bulk Analysis

- **Single sample:**  
  `telomerehunter2 -ibt tumor.bam -o out/ -p TumorID -b cytoband.txt`
- **Tumor vs Control:**  
  `telomerehunter2 -ibt tumor.bam -ibc control.bam -o out/ -p PairID -b cytoband.txt`
- **Custom repeats/species:**  
  `telomerehunter2 ... --repeats TTTAGGG TTAAGGG --repeatsContext TTAAGGG`
- **Fast mode (quick overview of unmapped reads generating summary with overview):**  
  `telomerehunter2 -ibt sample.bam -o out/ -p SampleID --fast_mode`

### Single cell sequencing Analysis

TelomereHunter2 now supports direct single-cell BAM analysis (with CB barcode tag). Simply run:

```bash
telomerehunter2_sc -ibt sample.bam -o results/ -p SampleID -b telomerehunter2/cytoband_files/cytoband.txt --min-reads-per-barcode 10000
```

This will perform barcode-aware telomere analysis and output per-cell results in a summary file. The minimum reads per barcode threshold can be set with `--min-reads-per-barcode`.

See `tests/test_telomerehunter2_sc.py` for example usage and validation.

## Input & Output

**Input:**  
- BAM/CRAM files (aligned reads)
- Cytoband file (tab-delimited, e.g. `hg19_cytoBand.txt`)
- Optional: custom telomeric repeats

**Output:**  
- `summary.tsv`, `TVR_top_contexts.tsv`, `singletons.tsv`
- Plots (`plots/` directory, PNG/HTML)
- Logs (run status/errors)
- For sc-seq: Additionally to the complete bulk run you get per-cell results in sc_summary.tsv and barcode_counts.tsv with reads counts per barcode

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
- **Banding file missing:** Needs reference genome banding file `-b` otherwise analysis will run without reads mapped to subtelomeres.
- **Plotting:** Try disabling with `--plotNone` or use plotting only mode with `--plotNone`.
- **Minor changes to TH1:** Skipping the tvrs normalization per 100 bp, improved detection of GXXGGG TVRs, read lengths are estimated from first 1000 reads, added TRPM


For help: [GitHub Issues](https://github.com/fpopp22/telomerehunter2/issues) or our FAQ.

## Documentation & Resources

- [Wiki](https://github.com/fpopp22/telomerehunter2/wiki) (WIP)
- [Example Data](tests/) 
- [Tutorial Videos](https://github.com/fpopp22/telomerehunter2/wiki) (WIP)
- [Original TelomereHunter Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0)

## Citation

If you use TelomereHunter2, please cite:
- Feuerbach, L., et al. "TelomereHunter – in silico estimation of telomere content and composition from cancer genomes." BMC Bioinformatics 20, 272 (2019). https://doi.org/10.1186/s12859-019-2851-0
- Application Note for TH2 (in preparation).

## Contributing

Fork, branch, and submit pull requests. Please add tests and follow code style. For major changes, open an issue first.

## License

GNU General Public License v3.0. See [LICENSE](LICENSE.txt).

## Contact

- Ferdinand Popp (f.popp@dkfz.de)
- Lars Feuerbach (l.feuerbach@dkfz.de)

## Acknowledgements

Developed by Ferdinand Popp, Lina Sieverling, Philip Ginsbach, Lars Feuerbach. Supported by German Cancer Research Center (DKFZ) - Division Applied Bioinformatics.

---

Copyright 2025 Ferdinand Popp, Lina Sieverling, Philip Ginsbach, Lars Feuerbach
