# Release Notes for v1.0.0

## Highlights
- Initial public release of TelomereHunter2 (Python 3)
- Major features:
  - Telomere content estimation from NGS/WGS data
  - Tumor/control analysis
  - Advanced plotting and PDF merging
  - Support for multiple cytoband files (dog, hg19, hg38, T2T)
  - Barcode splitting for single-cell ATAC-seq
- Dependencies: pysam, pandas, numpy, PyPDF2, plotly, kaleido
- CLI entry points: `telomerehunter2`, `sc-barcode-split-run-telomerehunter2`

## Improvements
- Modular codebase for easier maintenance
- Optional dependency for single-cell analysis

## Notes
- License: GPLv3
- For usage and documentation, see [README.md](README.md)

