# Quickstart

## Install

```bash
pip install telomerehunter2
```

## Bulk (single sample)

```bash
telomerehunter2 \
  -ibt sample.bam \
  -o results/ \
  -p SampleID \
  -b telomerehunter2/cytoband_files/hg19_cytoBand.txt
```

## Bulk (tumor vs control)

```bash
telomerehunter2 \
  -ibt tumor.bam \
  -ibc control.bam \
  -o results/ \
  -p PairID \
  -b telomerehunter2/cytoband_files/hg19_cytoBand.txt
```

## Single-cell

```bash
telomerehunter2-sc \
  -ibt sc_sample.bam \
  -o results/ \
  -p SampleID \
  -b telomerehunter2/cytoband_files/cytoband.txt \
  --min-reads-per-barcode 10000
```

## Next steps

- See [Use cases](usecases.md) for common patterns.
- See [FAQ](faq.md) for troubleshooting and performance tips.

