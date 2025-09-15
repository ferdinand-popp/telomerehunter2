# Installation and Usage

This page provides comprehensive installation instructions and usage guides for TelomereHunter2, supporting both bulk and single-cell telomere content analysis.

## System Requirements

### Operating System
- Linux (recommended)
- macOS 
- Windows (via WSL or Docker)

### Software Dependencies
- **Python**: ≥3.6 (Python 3.8+ recommended)
- **RAM**: 4GB minimum, 8GB+ recommended for large files
- **Storage**: 10GB free space for installation and temporary files

### Required Dependencies

The following Python packages are automatically installed:

- `pysam` - BAM/CRAM file handling
- `pandas` - Data manipulation and analysis
- `numpy` - Numerical computing
- `plotly` - Interactive visualizations  
- `PyPDF2` - PDF handling
- `kaleido` - Static image export (requires Chrome/Chromium)

## Installation Methods

### 1. PyPI Installation (Recommended)

Install the latest stable release from PyPI:

```bash
# Basic installation
pip install telomerehunter2

# Verify installation
telomerehunter2 --help
telomerehunter2-sc --help
```

### 2. Source Installation

#### With pip:
```bash
# Clone the repository
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/macOS
# venv\Scripts\activate   # Windows

# Install in development mode
pip install -e . --no-cache-dir

# Test installation
telomerehunter2 --help
```

#### With Poetry:
```bash
# Clone and setup
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2

# Install Poetry (if not already installed)
curl -sSL https://install.python-poetry.org | python3 -

# Setup environment and install
poetry env use python3
poetry install

# Activate environment
poetry shell

# Test installation
telomerehunter2 --help
```

#### With uv (Fast Python Package Manager):
```bash
# Clone repository
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2

# Install uv (if not already installed)
pip install uv

# Install TelomereHunter2
uv pip install -e . --no-cache-dir

# Test installation
telomerehunter2 --help
```

### 3. Docker Installation

#### Pull from Docker Hub (Recommended):
```bash
# Pull latest image
docker pull fpopp22/telomerehunter2

# Verify installation
docker run --rm fpopp22/telomerehunter2 telomerehunter2 --help
```

#### Build locally:
```bash
# Clone and build
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2

# Build Docker image  
docker build -t telomerehunter2 .

# Test installation
docker run --rm telomerehunter2 telomerehunter2 --help
```

### 4. Apptainer/Singularity Installation  

#### Pull from Docker Hub:
```bash
# Pull as Apptainer image
apptainer pull docker://fpopp22/telomerehunter2:latest

# Test installation
apptainer run telomerehunter2_latest.sif telomerehunter2 --help
```

#### Build from definition file:
```bash
# Clone repository
git clone https://github.com/ferdinand-popp/telomerehunter2.git
cd telomerehunter2

# Build from definition file
apptainer build telomerehunter2.sif Apptainer_TH2.def

# Test installation
apptainer run telomerehunter2.sif telomerehunter2 --help
```

## Quick Start Guides

### Bulk Analysis Quickstart

#### Single Sample Analysis

```bash
# Basic single sample run
telomerehunter2 \
  -ibt sample.bam \
  -o results/ \
  -p SampleID \
  -b /path/to/hg19_cytoBand.txt

# With custom settings
telomerehunter2 \
  -ibt sample.bam \
  -o results/ \
  -p SampleID \
  -b /path/to/hg19_cytoBand.txt \
  -rt 6 \
  -c 4 \
  --plotFileFormat pdf
```

#### Tumor vs Control Analysis

```bash
# Paired tumor-control analysis
telomerehunter2 \
  -ibt tumor.bam \
  -ibc control.bam \
  -o results/ \
  -p PatientID \
  -b /path/to/hg19_cytoBand.txt \
  -pl  # Enable parallel processing
```

#### Using Provided Reference Files

```bash
# Using included cytoband files
telomerehunter2 \
  -ibt sample.bam \
  -o results/ \
  -p SampleID \
  -b telomerehunter2/cytoband_files/hg19_cytoBand.txt

# For hg38
telomerehunter2 \
  -ibt sample.bam \
  -o results/ \
  -p SampleID \
  -b telomerehunter2/cytoband_files/hg38_cytoBand.txt
```

### Single-Cell Analysis Quickstart

```bash
# Basic single-cell analysis
telomerehunter2-sc \
  -ibt scATAC_sample.bam \
  -o sc_results/ \
  -p ScSampleID \
  -b telomerehunter2/cytoband_files/hg38_cytoBand.txt \
  --min-reads-per-barcode 10000

# With custom barcode threshold
telomerehunter2-sc \
  -ibt scATAC_sample.bam \
  -o sc_results/ \
  -p ScSampleID \
  -b telomerehunter2/cytoband_files/hg38_cytoBand.txt \
  --min-reads-per-barcode 5000 \
  -c 8
```

### Docker Usage Examples

#### Basic Docker run:
```bash
# Single sample with Docker
docker run --rm -it \
  -v /path/to/data:/data \
  fpopp22/telomerehunter2 \
  telomerehunter2 \
    -ibt /data/sample.bam \
    -o /data/results \
    -p SampleID \
    -b /data/hg19_cytoBand.txt
```

#### Docker with custom settings:
```bash
# Advanced Docker run with multiple options
docker run --rm -it \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  fpopp22/telomerehunter2 \
  telomerehunter2 \
    -ibt /data/tumor.bam \
    -ibc /data/control.bam \
    -o /results \
    -p PatientID \
    -b /data/hg38_cytoBand.txt \
    -rt 6 \
    -c 4 \
    -pl \
    --plotFileFormat all
```

### Apptainer Usage Examples

```bash
# Basic Apptainer run
apptainer run \
  -B /path/to/data:/data \
  telomerehunter2_latest.sif \
  telomerehunter2 \
    -ibt /data/sample.bam \
    -o /data/results \
    -p SampleID \
    -b /data/hg19_cytoBand.txt

# Single-cell with Apptainer
apptainer run \
  -B /path/to/data:/data \
  telomerehunter2_latest.sif \
  telomerehunter2-sc \
    -ibt /data/scATAC.bam \
    -o /data/sc_results \
    -p ScID \
    -b /data/hg38_cytoBand.txt \
    --min-reads-per-barcode 8000
```

## Common Usage Patterns

### Analysis Workflow

1. **Prepare input files**: Ensure BAM/CRAM files are indexed
2. **Choose cytoband file**: Select appropriate reference genome
3. **Run analysis**: Execute TelomereHunter2 with desired parameters
4. **Review results**: Check summary.tsv and plots
5. **Interpret findings**: Use output documentation for guidance

### Typical Command Structure

```bash
telomerehunter2 \
  -ibt <input_bam> \          # Required: input BAM/CRAM
  [-ibc <control_bam>] \      # Optional: control sample  
  -o <output_dir> \           # Required: output directory
  -p <sample_id> \            # Required: sample identifier
  -b <cytoband_file> \        # Recommended: cytoband file
  [additional_options]        # Optional: customization flags
```

## Performance Optimization

### CPU and Memory Settings

```bash
# Use multiple cores (default: all available)
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -c 8

# Enable parallel processing for tumor-control pairs
telomerehunter2 -ibt tumor.bam -ibc control.bam -o results/ -p ID -b cytoband.txt -pl

# Subsample for faster testing (10% of reads)
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt --subsample 0.1
```

### Large File Handling

```bash
# For large BAM files, consider:
# 1. Increase memory allocation (if using containers)
# 2. Use faster storage (SSD)
# 3. Enable parallel processing
# 4. Remove duplicates to reduce data size

telomerehunter2 \
  -ibt large_sample.bam \
  -o results/ \
  -p LargeSample \
  -b cytoband.txt \
  -d \     # Remove duplicates
  -c 16 \  # Use more cores
  -pl      # Enable parallel processing
```

## Troubleshooting

### Common Installation Issues

#### Missing Dependencies
```bash
# Error: No module named 'pysam'
pip install pysam pandas numpy plotly PyPDF2 kaleido

# Error: Chrome/Chromium not found (for static plots)
# Linux:
sudo apt-get install chromium-browser
# macOS:
brew install chromium
```

#### Permission Issues  
```bash
# Use user installation if root access unavailable
pip install --user telomerehunter2

# Or use virtual environment
python -m venv th2_env
source th2_env/bin/activate
pip install telomerehunter2
```

### Common Runtime Issues

#### Memory Errors
```bash
# Reduce parallel processes
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt -c 2

# Use subsampling for testing
telomerehunter2 -ibt sample.bam -o results/ -p ID -b cytoband.txt --subsample 0.5
```

#### BAM File Issues  
```bash
# Check BAM file integrity
samtools quickcheck sample.bam

# Index BAM file if needed
samtools index sample.bam

# Check for required tags (single-cell)
samtools view sample.bam | head -1000 | grep -c "CB:" 
```

#### Missing Cytoband File
```bash
# Analysis runs without cytoband but misses subtelomeric regions
# Warning: "No banding file specified"
# Solution: Download appropriate cytoband file or use provided ones
```

### Docker/Container Issues

#### Permission/Mount Issues
```bash
# Ensure proper volume mounting and permissions
docker run --rm -it \
  -v /path/to/data:/data:ro \     # Read-only data mount
  -v /path/to/results:/results \  # Write access for results
  --user $(id -u):$(id -g) \     # Use current user ID
  fpopp22/telomerehunter2 \
  telomerehunter2 -ibt /data/sample.bam -o /results -p ID -b /data/cytoband.txt
```

### Getting Help

#### Check Version and Help
```bash
# Version information
telomerehunter2 --version

# Detailed help
telomerehunter2 --help
telomerehunter2-sc --help

# Test with example data (if available)
cd tests/
python -m pytest test_telomerehunter2.py -v
```

#### Error Logs
Check log files in the output directory for detailed error information:
- `logs/error.log` - Error messages
- `logs/analysis.log` - Analysis progress
- `logs/filtering.log` - Read filtering details

## Best Practices

1. **Always index BAM/CRAM files** before analysis
2. **Use appropriate cytoband files** for your reference genome
3. **Start with default settings** before customizing
4. **Check available memory** before processing large files
5. **Use containers** for reproducible analysis environments
6. **Backup important results** before rerunning analysis

## Next Steps

After successful installation:

- Review [Input and Output](Input-and-Output.md) documentation
- Explore [Advanced Options and Customization](Advanced-Options-and-Customization.md)
- Check out [Contributing and Development](Contributing-and-Development.md) if interested in development

## Support

For technical support:
- [GitHub Issues](https://github.com/ferdinand-popp/telomerehunter2/issues)
- Email: f.popp@dkfz-heidelberg.de
- [Original TelomereHunter Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0)