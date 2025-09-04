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

## Installation

**Classic setup (bulk analysis):**  
```bash
pip install telomerehunter2
```

**Single-cell setup (sc analysis):**  
Requires [sinto](https://github.com/timoast/sinto) for barcode splitting:  
```bash
pip install 'telomerehunter2[sc]'
```

**From source:**  
```bash
# With pip:
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2
python -m venv venv
source venv/bin/activate
pip install -e . --no-cache-dir

# For single-cell support:
pip install -e .[sc] --no-cache-dir

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

```bash
telomerehunter2 -ibt sample.bam -o results/ -p SampleID -b telomerehunter2/cytoband_files/hg19_cytoBand.txt
```
For all options:  
```bash
telomerehunter2 --help
```

## Usage

### Bulk Analysis

- **Single sample:**  
  `telomerehunter2 -ibt tumor.bam -o out/ -p TumorID -b cytoband.txt`
- **Tumor vs Control:**  
  `telomerehunter2 -ibt tumor.bam -ibc control.bam -o out/ -p PairID -b cytoband.txt`
- **Custom repeats/species:**  
  `telomerehunter2 ... --repeats TTTAGGG TTAAGGG --repeatsContext TTAAGGG`

### Single cell sequencing Analysis

Requires BAM with CB barcode tag and Sinto for splitting:  
Install with `pip install 'telomerehunter2[sc]'`  
```bash
python sc_barcode_splitter_run.py -ibt input.bam -b cytoband.txt -p PatientID -o out/ --keep-bams
```
See `tests/run_sc_atac.py` for examples.

## Input & Output

**Input:**  
- BAM/CRAM files (aligned reads)
- Cytoband file (tab-delimited, e.g. `hg19_cytoBand.txt`)
- Optional: custom telomeric repeats

**Output:**  
- `summary.tsv`, `TVR_top_contexts.tsv`, `singletons.tsv`
- Plots (`plots/` directory, PNG/HTML)
- Logs (run status/errors)
- For sc-seq: per-cell results in subfolders and barcode file with barcodes over read threshold

## Dependencies

- Python >=3.6
- pysam, numpy, pandas, plotly, PyPDF2
- For static image export: kaleido (requires chrome/chromium)
- Sinto (for sc-seq, install with `[sc]` extra)
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

For help: [GitHub Issues](https://github.com/fpopp22/telomerehunter2/issues) or our FAQ.

## Documentation & Resources

- [Wiki](https://github.com/fpopp22/telomerehunter2/wiki)
- [Example Data](tests/)
- [Tutorial Videos](https://github.com/fpopp22/telomerehunter2/wiki)
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
