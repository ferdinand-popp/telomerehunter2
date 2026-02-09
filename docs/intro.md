# Intro

TelomereHunter2 detects telomeric and telomere-adjacent reads in genome sequencing data and summarizes telomere content as a normalized metric. It also reports telomeric variant repeat (TVR) composition.

## Inputs

- **BAM/CRAM**: aligned sequencing reads
- **Cytoband file (`-b`)**: used to determine the first/last band per chromosome (subtelomeric context)
- **Sample / pair ID (`-p`)**: identifier used for file names and report labels

## Bulk vs single-cell

- **Bulk** (`telomerehunter2`): processes one BAM/CRAM (optionally tumor+control) and generates summary tables and plots.
- **Single-cell** (`telomerehunter2-sc`): processes a BAM with barcode tags and reports per-cell metrics.

## Outputs (high level)

- `*_summary.tsv`: telomere content and read category counts
- `*_TVR_top_contexts.tsv`: top TVR contexts
- `*_singletons.tsv`: singleton TVR statistics
- `plots/` and `html_reports/`: figures and HTML reports (if enabled)

