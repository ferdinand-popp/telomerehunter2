# FAQ

## Technical issues

### Where do I get a cytoband file?

TelomereHunter2 ships cytoband files under `telomerehunter2/cytoband_files/`. Use the one matching your reference build (e.g. hg19 vs GRCh38) when available.

### Can I run without `-b` cytobanding file?

Yes. The tool can run without a cytoband file, but subtelomeric/junction classification may be reduced as we can discern between reads in the different regions.

### The run is slow / uses too much RAM. What can I do?

- Reduce cores via `-c` (if available in your version)
- Use `--fast_mode` for a quick overview
- Ensure input BAM/CRAM is indexed and stored on fast storage

### How do I change the barcode tag for single-cell?

Use `--barcodeTag` to match the tag used in your BAM (default is typically `CB`).

### Plot export fails (kaleido/chrome)

Static image export via kaleido may require chrome/chromium. If you only need TSV outputs, disable plotting (option name depends on your run mode; see `--help`).

## Experimental Issues

### How should I correct log2 T/C values for tumour purity in TelomereHunter2?

#### Does TelomereHunter2 include built-in purity correction?

No. Built-in purity correction is not implemented because the necessary information is typically unavailable in standard bulk WGS data. The blood control sample used as a reference is only a poor proxy for the non-tumour cells actually present in the tumour sample, which include fibroblasts, infiltrating immune cells, and endothelial cells, each with distinct telomere length distributions.

#### Is linear unmixing an appropriate correction approach?

It is an available approach for bulk data, but it oversimplifies the biology. A linear correction assumes the blood control equals the non-tumour telomere content in the tumour sample, which does not hold in practice. We therefore recommend using the blood control primarily to correct for age effects and inherited telomere length differences, rather than for purity correction. Single-cell sequencing methods provide more reliable cell-type-resolved estimates when purity correction is critical (see [preprint](https://www.biorxiv.org/content/10.1101/2024.08.28.609339v2.full)).

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

