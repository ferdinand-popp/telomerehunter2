# Advanced Options and Customization

This page covers advanced configuration options, customization techniques, and specialized use cases for TelomereHunter2.

## Command-Line Options Reference

### Input and Output Options

#### Basic Input/Output
```bash
-ibt, --inputBamTumor        # Path to tumor/sample BAM/CRAM file
-ibc, --inputBamControl      # Path to control BAM/CRAM file (optional)
-o, --outPath                # Output directory path  
-p, --pid                    # Sample/project identifier
-b, --bandingFile           # Cytoband file path (optional but recommended)
```

#### Advanced Input Options
```bash
--subsample FRACTION         # Subsample reads (0.0-1.0), e.g., 0.1 for 10%
-nf, --noFiltering          # Skip filtering if already performed
```

### Threshold and Quality Control Options

#### Repeat Thresholds
```bash
-rt, --repeatThreshold N     # Minimum repeats for telomeric classification (default: auto)
-rl, --perReadLength        # Set repeat threshold per 100bp read length
-con, --consecutive         # Search for consecutive repeats only
```

**Examples:**
```bash
# Fixed threshold of 6 repeats
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -rt 6

# Length-dependent threshold (6 per 100bp)
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -rt 6 -rl

# Consecutive repeats only
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -con
```

#### Quality Filtering
```bash
-mqt, --mappingQualityThreshold N   # Mapping quality threshold (default: 8)
-d, --removeDuplicates             # Remove duplicate reads
```

**Example:**
```bash
# High-quality reads only
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt \
  -mqt 20 -d
```

### GC Correction Options

```bash
-gc1, --lowerGC PERCENT      # Lower GC limit for correction
-gc2, --upperGC PERCENT      # Upper GC limit for correction
```

**Example:**
```bash
# Custom GC correction range
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt \
  -gc1 45 -gc2 55
```

### Performance and Parallelization

#### Core and Processing Options
```bash
-c, --cores N               # Number of CPU cores to use
-pl, --parallel            # Enable parallel processing (tumor vs control)
```

**Examples:**
```bash
# Use 8 CPU cores
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -c 8

# Parallel tumor-control processing
telomerehunter2 -ibt tumor.bam -ibc control.bam -o results/ -p ID \
  -b cytoband.txt -pl -c 16
```

## Custom Telomeric Repeats

### Default Repeat Sequences

TelomereHunter2 analyzes these sequences by default:
- `TTAGGG` - Canonical human telomeric repeat
- `TGAGGG` - Telomeric variant repeat (TVR)
- `TCAGGG` - TVR
- `TTCGGG` - TVR  
- `TTGGGG` - TVR

### Custom Repeat Configuration

#### Basic Custom Repeats
```bash
-r, --repeats SEQUENCE [SEQUENCE ...]           # Custom repeat sequences
-rc, --repeatsContext SEQUENCE [SEQUENCE ...]  # TVRs for context analysis
-bp, --bpContext N                             # Context window size (divisible by 6)
```

**Examples:**

```bash
# Plant telomeres (Arabidopsis)
telomerehunter2 -ibt plant.bam -o results/ -p PlantID \
  --repeats TTTAGGG TTAAGGG TTAGGG

# Mouse telomeres with context analysis
telomerehunter2 -ibt mouse.bam -o results/ -p MouseID \
  --repeats TTAGGG TCAGGG TGAGGG \
  --repeatsContext TCAGGG TGAGGG \
  --bpContext 18

# Custom organism with unique repeats
telomerehunter2 -ibt custom.bam -o results/ -p CustomID \
  --repeats GGGTTT AGGTTT CGGTTT \
  --repeatsContext AGGTTT CGGTTT
```

### Non-Human Genome Analysis

#### Using Custom Cytoband Files

For non-human genomes, create custom cytoband files or use existing ones:

```bash
# Dog genome analysis
telomerehunter2 -ibt dog.bam -o results/ -p DogID \
  -b telomerehunter2/cytoband_files/cytoBand_dog.txt \
  --repeats TTAGGG TCAGGG

# Analysis without cytoband (no subtelomeric analysis)
telomerehunter2 -ibt organism.bam -o results/ -p OrganismID \
  --repeats TTAGGG TCAGGG TGAGGG
```

#### Creating Custom Cytoband Files

Format: Tab-separated, no header
```
chromosome    start    end    band_name    band_type
chr1          0        1000000    p1.1        gneg
chr1          1000000  2000000    p1.2        gpos25
```

**Example creation:**
```bash
# Create cytoband file from reference
awk 'BEGIN{OFS="\t"} {print $1, 0, $2, "p1", "gneg"}' genome.fai > custom_cytoband.txt
```

## Single-Cell Analysis Customization

### Barcode Filtering Options

```bash
--min-reads-per-barcode N    # Minimum reads per cell barcode (default: 10000)
-sc, --singlecell_mode      # Force single-cell mode
```

**Examples:**
```bash
# Low-depth single-cell analysis
telomerehunter2-sc -ibt sc_sample.bam -o sc_results/ -p ScID \
  -b cytoband.txt --min-reads-per-barcode 1000

# High-quality cells only
telomerehunter2-sc -ibt sc_sample.bam -o sc_results/ -p ScID \
  -b cytoband.txt --min-reads-per-barcode 50000

# Force single-cell mode (even without CB tags)
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -sc
```

### Single-Cell Quality Control

```bash
# Combine quality filters for single-cell
telomerehunter2-sc -ibt scATAC.bam -o sc_results/ -p ScID \
  -b cytoband.txt \
  --min-reads-per-barcode 5000 \
  -mqt 15 \
  -d \
  -rt 6 \
  -c 8
```

## Plotting and Visualization Customization

### Plot Control Options

```bash
-pno, --plotNone            # Disable all plots
--plot_mode                 # Run only plotting on existing results
-pff, --plotFileFormat      # Output format: pdf, png, svg, all
```

### Individual Plot Types

```bash
-p1, --plotChr              # Chromosome distribution plots
-p2, --plotFractions        # Telomere fraction plots  
-p3, --plotTelContent       # Telomere content plots
-p4, --plotGC              # GC distribution plots
-p5, --plotRepeatFreq      # Repeat frequency histograms
-p6, --plotTVR             # TVR analysis plots
-p7, --plotSingleton       # Singleton TVR plots
-prc, --plotRevCompl       # Distinguish forward/reverse complement
```

**Examples:**

```bash
# Generate only essential plots
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt \
  -p3 -p4 -p6 -pff pdf

# All plots in multiple formats
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt \
  -p1 -p2 -p3 -p4 -p5 -p6 -p7 -prc -pff all

# Plots only mode (rerun plotting with different settings)
telomerehunter2 -ibt sample.bam -o existing_results/ -p ID -b cytoband.txt \
  --plot_mode -p6 -p7 -pff png

# No plots (fastest analysis)
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -pno
```

## Advanced Analysis Scenarios

### Tumor vs Control Comparisons

```bash
# Basic paired analysis
telomerehunter2 \
  -ibt tumor.bam \
  -ibc normal.bam \
  -o paired_results/ \
  -p PatientID \
  -b hg38_cytoBand.txt \
  -pl

# High-quality paired analysis with custom settings
telomerehunter2 \
  -ibt tumor.bam \
  -ibc normal.bam \
  -o paired_results/ \
  -p PatientID \
  -b hg38_cytoBand.txt \
  -pl \
  -rt 6 \
  -mqt 15 \
  -d \
  -c 12 \
  --plotFileFormat all
```

### Batch Processing

#### Process Multiple Samples
```bash
#!/bin/bash
# batch_analysis.sh

SAMPLES=("sample1" "sample2" "sample3")
CYTOBAND="/path/to/hg38_cytoBand.txt"
OUTPUT_DIR="/path/to/results"

for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."
    telomerehunter2 \
        -ibt "${sample}.bam" \
        -o "$OUTPUT_DIR/${sample}_results" \
        -p "$sample" \
        -b "$CYTOBAND" \
        -c 8 \
        -pno  # Skip plots for batch processing
done
```

#### Paired Sample Processing
```bash
#!/bin/bash
# batch_paired_analysis.sh

PATIENTS=("Patient01" "Patient02" "Patient03")
CYTOBAND="/path/to/hg38_cytoBand.txt"
OUTPUT_DIR="/path/to/results"

for patient in "${PATIENTS[@]}"; do
    echo "Processing $patient..."
    telomerehunter2 \
        -ibt "${patient}_tumor.bam" \
        -ibc "${patient}_normal.bam" \
        -o "$OUTPUT_DIR/${patient}_results" \
        -p "$patient" \
        -b "$CYTOBAND" \
        -pl \
        -c 16
done
```

### Time-Course and Longitudinal Analysis

```bash
# Process multiple time points
for timepoint in T0 T1 T2 T3; do
    telomerehunter2 \
        -ibt "patient_${timepoint}.bam" \
        -o "results/patient_${timepoint}" \
        -p "Patient_${timepoint}" \
        -b hg38_cytoBand.txt \
        -rt 6 \
        -c 8
done
```

### Multi-Species Comparative Analysis

```bash
# Human sample
telomerehunter2 -ibt human.bam -o human_results/ -p Human \
  -b hg38_cytoBand.txt

# Mouse sample with adjusted repeats
telomerehunter2 -ibt mouse.bam -o mouse_results/ -p Mouse \
  --repeats TTAGGG TCAGGG TGAGGG

# Plant sample with custom repeats
telomerehunter2 -ibt plant.bam -o plant_results/ -p Plant \
  --repeats TTTAGGG TTAGGG --repeatsContext TTTAGGG
```

## Container-Based Advanced Usage

### Docker Advanced Examples

#### Resource-Limited Environment
```bash
# Limit CPU and memory usage
docker run --rm -it \
  --cpus="4" \
  --memory="8g" \
  -v /data:/data \
  fpopp22/telomerehunter2 \
  telomerehunter2 -ibt /data/sample.bam -o /data/results -p ID \
  -b /data/cytoband.txt -c 4
```

#### Custom Analysis Pipeline
```bash
# Multi-stage analysis with intermediate results
docker run --rm -it \
  -v /data:/data \
  fpopp22/telomerehunter2 \
  bash -c "
    telomerehunter2 -ibt /data/sample.bam -o /data/stage1 -p S1 \
      -b /data/cytoband.txt -nf &&
    telomerehunter2 -ibt /data/sample.bam -o /data/stage2 -p S2 \
      -b /data/cytoband.txt --plot_mode -pff all
  "
```

### Apptainer/Singularity Advanced Usage

#### HPC Cluster Integration
```bash
# SLURM job script example
#!/bin/bash
#SBATCH --job-name=telomerehunter2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00

module load apptainer

apptainer run \
  -B /scratch:/scratch \
  -B /data:/data \
  telomerehunter2.sif \
  telomerehunter2 \
    -ibt /data/large_sample.bam \
    -o /scratch/results \
    -p LargeSample \
    -b /data/hg38_cytoBand.txt \
    -c 16 \
    -pl
```

## Performance Tuning

### Memory Optimization

```bash
# For large files (>50GB BAM)
telomerehunter2 -ibt large.bam -o results/ -p Large -b cytoband.txt \
  -c 8 \           # Moderate core count
  --subsample 0.8  # Reduce data if needed

# Memory-efficient single-cell analysis
telomerehunter2-sc -ibt sc_large.bam -o sc_results/ -p ScLarge \
  -b cytoband.txt \
  --min-reads-per-barcode 20000 \  # Higher threshold
  -c 4                             # Fewer cores
```

### Speed Optimization

```bash
# Fastest analysis (minimal plots, parallel processing)
telomerehunter2 -ibt tumor.bam -ibc control.bam -o results/ -p Fast \
  -b cytoband.txt \
  -pl \      # Parallel processing
  -pno \     # No plots
  -c 32      # Maximum cores

# Quick preview analysis  
telomerehunter2 -ibt sample.bam -o preview/ -p Preview -b cytoband.txt \
  --subsample 0.1 \  # 10% of reads
  -pno               # No plots
```

## Troubleshooting Advanced Scenarios

### Complex BAM Files

```bash
# Multiple read groups
telomerehunter2 -ibt multi_rg.bam -o results/ -p MultiRG -b cytoband.txt

# Mixed single/paired-end reads  
telomerehunter2 -ibt mixed.bam -o results/ -p Mixed -b cytoband.txt \
  -mqt 10  # Lower mapping quality threshold

# Very long reads (PacBio/Nanopore)
telomerehunter2 -ibt long_reads.bam -o results/ -p LongReads -b cytoband.txt \
  -rt 20   # Higher repeat threshold
```

### Custom Error Handling

```bash
# Robust analysis with error recovery
telomerehunter2 -ibt sample.bam -o results/ -p Robust -b cytoband.txt \
  --subsample 1.0 || \
telomerehunter2 -ibt sample.bam -o results_backup/ -p Robust -b cytoband.txt \
  --subsample 0.5 -pno  # Fallback with reduced data
```

## Integration with Other Tools

### Preprocessing Integration

```bash
# With samtools preprocessing
samtools view -h -q 10 -F 1024 input.bam | \
samtools sort -o filtered.bam - &&
samtools index filtered.bam &&
telomerehunter2 -ibt filtered.bam -o results/ -p Filtered -b cytoband.txt
```

### Post-processing Integration

```bash
# Combine multiple results
telomerehunter2 -ibt sample1.bam -o results1/ -p S1 -b cytoband.txt
telomerehunter2 -ibt sample2.bam -o results2/ -p S2 -b cytoband.txt

# Merge summary files (example)
cat results1/*/summary.tsv results2/*/summary.tsv > combined_summary.tsv
```

## Best Practices for Advanced Usage

1. **Test with subsamples** before full analysis on large datasets
2. **Monitor system resources** during parallel processing
3. **Use consistent parameters** across comparative studies
4. **Validate custom repeat sequences** with known positive controls
5. **Document parameter choices** for reproducible analysis
6. **Backup intermediate results** for long-running analyses

## Further Reading

- [Input and Output](Input-and-Output.md) - Detailed file format documentation
- [Installation and Usage](Installation-and-Usage.md) - Basic setup and usage
- [Contributing and Development](Contributing-and-Development.md) - Development guidelines
- [References and Citations](References-and-Citations.md) - Academic references