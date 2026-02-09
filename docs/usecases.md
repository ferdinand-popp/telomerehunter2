# Use cases

## 1) Quick telomere overview (fast mode)

Use fast mode to quickly summarize unmapped reads.

```bash
telomerehunter2 -ibt sample.bam -o results/ -p SampleID --fast_mode
```

## 2) Non-human / custom repeats

Provide custom telomeric repeats and (optionally) the repeat context.

```bash
telomerehunter2 \
  -ibt sample.bam \
  -o results/ \
  -p SampleID \
  --repeats TTTAGGG TTAAGGG \
  --repeatsContext TTAAGGG
```

## 3) Run without cytoband file

A cytoband file improves subtelomeric classification. If you do not provide `-b`, analysis can still run, but subtelomeric/junction assignments may be limited.

```bash
telomerehunter2 -ibt sample.bam -o results/ -p SampleID
```

## 4) Single-cell: rerun postprocessing with different filtering

If you want to adjust `--min-reads-per-barcode` without repeating the expensive filtering step, rerun with `--noFiltering`.

```bash
telomerehunter2-sc \
  -ibt sc_sample.bam \
  -o results/ \
  -p SampleID \
  -b telomerehunter2/cytoband_files/cytoband.txt \
  --min-reads-per-barcode 5000 \
  --noFiltering
```

## 5) Containers

If you run via Docker/Apptainer, the command-line interface is identical; ensure your input and output paths are mounted into the container.

