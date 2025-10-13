# TelomereHunter2 Release Notes

## v1.0.3 (2025-10-13)

### Highlights
- Added chromium to docker and apptainer files for kaleido to generate pdfs + test case
- enhance parallel processing by dynamically adjusting core usage and adding warnings for exceeding available cores

## v1.0.2 (2025-10-10)

### Highlights
- Added documentation for interpretation of summary.tsv in README.md
- refactored summary.tsv columns order and names
- added warning for runs that show no unmapped reads
- plotting default set to ['pdf', 'html]', of which html can be generated with plotly and pdf if kaleido can access chromium, errors are now caught
- version retrieval in utils.py now with regex as toml dependency varied for python versions
- removed measure_time decorator from development
- general code cleanup


## v1.0.1 (2025-10-02)

### Highlights
- Fixed region info string joining by using tuple and double underscore separator; region fetching now uses contig, start, and stop directly
- Removed subsample code from region processing and streamlined preprocessing
- Implemented fast mode for tumor and control samples with summary generation and readme in separate file
- Improved single-cell mode and barcode counting in telomere read filtering by removing sinto dependency and custom file logic
- Single cell mode is invoked by `telomerehunter2_sc` command and runs pseudo bulk and sc custom code
- General code cleanup and documentation improvements


## v1.0.0 (Initial release)

### Highlights
- Initial public release of TelomereHunter2 (Python 3)
- Major features:
  - Telomere content estimation from NGS/WGS data
  - Tumor/control analysis
  - Advanced plotting and PDF merging
  - Support for multiple cytoband files (dog, hg19, hg38, T2T)
  - Barcode splitting for single-cell ATAC-seq
- Dependencies: pysam, pandas, numpy, PyPDF2, plotly, kaleido
- CLI entry points: `telomerehunter2`, `sc-barcode-split-run-telomerehunter2`
