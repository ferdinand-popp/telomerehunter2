# FAQ

## Where do I get a cytoband file?

TelomereHunter2 ships cytoband files under `telomerehunter2/cytoband_files/`. Use the one matching your reference build (e.g. hg19 vs GRCh38) when available.

## Can I run without `-b`?

Yes. The tool can run without a cytoband file, but subtelomeric/junction classification may be reduced.

## The run is slow / uses too much RAM. What can I do?

- Reduce cores via `-c` (if available in your version)
- Use `--fast_mode` for a quick overview
- Ensure input BAM/CRAM is indexed and stored on fast storage

## How do I change the barcode tag for single-cell?

Use `--barcodeTag` to match the tag used in your BAM (default is typically `CB`).

## Plot export fails (kaleido/chrome)

Static image export via kaleido may require chrome/chromium. If you only need TSV outputs, disable plotting (option name depends on your run mode; see `--help`).
# TelomereHunter2

TelomereHunter2 is a Python tool for estimating telomere content and analyzing telomeric variant repeats (TVRs) from genome sequencing data.

## What you can do

- Bulk telomere content estimation from BAM/CRAM inputs
- Tumor vs control comparison (log2 ratio)
- Single-cell BAM analysis (barcode-aware)
- Generate static/interactive reports

## Where to start

- Read the [Intro](intro.md) to understand the concepts.
- Follow the [Quickstart](quickstart.md) to run your first analysis.
- Check [Use cases](usecases.md) for common workflows.
- See [FAQ](faq.md) if you run into issues.

