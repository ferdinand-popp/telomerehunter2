# TelomereHunter2 Wiki

Welcome to the comprehensive documentation for **TelomereHunter2**, an advanced Python tool for estimating telomere content and analyzing telomeric variant repeats (TVRs) from genome sequencing data.

## Quick Navigation

### 📚 **Documentation Pages**

| Page | Description |
|------|-------------|
| [🔗 Installation and Usage](Installation-and-Usage.md) | Complete setup instructions and quickstart guides |
| [📊 Input and Output](Input-and-Output.md) | Detailed file format specifications and output explanations |
| [⚙️ Advanced Options and Customization](Advanced-Options-and-Customization.md) | Advanced configuration and specialized use cases |
| [👥 Contributing and Development](Contributing-and-Development.md) | Developer guidelines and project contribution |
| [📖 References and Citations](References-and-Citations.md) | Academic citations and acknowledgements |

### 🚀 **Quick Start**

New to TelomereHunter2? Start here:

1. **[Install TelomereHunter2](Installation-and-Usage.md#installation-methods)** - Choose your installation method
2. **[Run your first analysis](Installation-and-Usage.md#bulk-analysis-quickstart)** - Basic usage examples
3. **[Understand the output](Input-and-Output.md#main-output-files)** - Interpret your results
4. **[Explore advanced options](Advanced-Options-and-Customization.md)** - Customize your analysis

## What is TelomereHunter2?

TelomereHunter2 is a fast, container-friendly Python 3 implementation for analyzing telomere content from sequencing data. It provides:

### ✨ **Key Features**
- **Telomere Content Estimation**: GC-corrected telomere content from WGS data
- **TVR Analysis**: Comprehensive telomeric variant repeat analysis
- **Single-Cell Support**: Barcode-aware analysis for scATAC-seq and similar data
- **Multi-Species**: Support for human and non-human genomes
- **Containerized**: Docker and Apptainer/Singularity support
- **Interactive Visualizations**: HTML reports with Plotly integration

### 🔬 **Analysis Types**
- **Bulk Sequencing**: Standard genome-wide telomere analysis
- **Single-Cell**: Per-barcode telomere content estimation
- **Tumor vs Control**: Comparative analysis between paired samples
- **Custom Species**: Configurable repeat sequences for any organism

## Analysis Overview

### Input Requirements
- **BAM/CRAM files**: Aligned sequencing reads (indexed)
- **Cytoband files**: Chromosome banding information (optional but recommended)
- **Custom repeats**: Telomeric sequences for non-human species (optional)

### Output Products
- **summary.tsv**: Main results with telomere content estimates
- **Interactive plots**: HTML visualizations and static images  
- **TVR analysis**: Detailed telomeric variant repeat patterns
- **Single-cell results**: Per-barcode analysis (sc mode)
- **Quality logs**: Comprehensive analysis logs

## Usage Examples

### Basic Analysis
```bash
# Single sample analysis
telomerehunter2 -ibt sample.bam -o results/ -p SampleID \
  -b hg38_cytoBand.txt

# Tumor vs control
telomerehunter2 -ibt tumor.bam -ibc control.bam -o results/ -p PatientID \
  -b hg38_cytoBand.txt -pl
```

### Single-Cell Analysis
```bash
# Single-cell ATAC-seq analysis  
telomerehunter2-sc -ibt scATAC.bam -o sc_results/ -p ScSampleID \
  -b hg38_cytoBand.txt --min-reads-per-barcode 10000
```

### Container Usage
```bash
# Docker analysis
docker run --rm -v /data:/data fpopp22/telomerehunter2 \
  telomerehunter2 -ibt /data/sample.bam -o /data/results -p ID \
  -b /data/hg38_cytoBand.txt
```

## Scientific Background

### Telomere Biology
Telomeres are protective DNA-protein structures at chromosome ends, consisting of repetitive sequences (TTAGGG in humans). TelomereHunter2 analyzes:

- **Telomere content**: Overall abundance of telomeric DNA
- **Telomeric Variant Repeats (TVRs)**: Alternative repeat sequences
- **GC correction**: Accounts for sequence composition bias
- **Subtelomeric regions**: Chromosome-specific analysis

### Clinical Relevance
Telomere length analysis is important for:
- **Cancer research**: Telomere maintenance mechanisms
- **Aging studies**: Cellular senescence markers  
- **Genetic diseases**: Telomeropathies and related disorders
- **Therapeutic targets**: Telomerase inhibition strategies

## Community and Support

### Getting Help
- **📖 Documentation**: Comprehensive Wiki pages
- **🐛 Issues**: [GitHub Issues](https://github.com/ferdinand-popp/telomerehunter2/issues) for bug reports
- **💬 Discussions**: Community questions and feature requests  
- **📧 Direct Contact**: f.popp@dkfz-heidelberg.de

### Contributing
We welcome contributions! See our [Contributing Guidelines](Contributing-and-Development.md) for:
- Code contributions and bug fixes
- Documentation improvements
- Feature suggestions and testing
- Scientific applications and validation

## Installation Options

Choose the installation method that works best for your environment:

| Method | Best For | Command |
|--------|----------|---------|
| **PyPI** | Most users | `pip install telomerehunter2` |
| **Docker** | Containers, reproducibility | `docker pull fpopp22/telomerehunter2` |
| **Source** | Development, customization | `git clone` + `pip install -e .` |
| **Apptainer** | HPC environments | `apptainer pull docker://fpopp22/telomerehunter2` |

## Performance and Scalability

### System Requirements
- **Memory**: 4GB minimum, 8GB+ recommended
- **CPU**: Multi-core support with parallelization
- **Storage**: 10GB+ for installation and temp files
- **Python**: 3.6+ (3.8+ recommended)

### Optimization Tips
- Use parallel processing (`-pl`) for tumor-control pairs
- Adjust core count (`-c N`) based on available resources
- Consider subsampling (`--subsample`) for large datasets
- Use containers for consistent performance across systems

## Version Information

- **Current Version**: 1.0.0
- **License**: GNU General Public License v3.0
- **Python Support**: 3.6+
- **Container Images**: Available on Docker Hub

## Citation

If you use TelomereHunter2 in your research, please cite:

**Feuerbach, L., Sieverling, L., Deeg, K.I. et al.** TelomereHunter – in silico estimation of telomere content and composition from cancer genomes. *BMC Bioinformatics* **20**, 272 (2019). https://doi.org/10.1186/s12859-019-2851-0

*Additional citation for TelomereHunter2 application note will be provided when available.*

## Recent Updates

- **v1.0.0**: Major Python 3 rewrite with enhanced features
- **Single-cell support**: New `telomerehunter2-sc` command
- **Container optimization**: Improved Docker and Apptainer images
- **Documentation**: Comprehensive Wiki documentation
- **Performance**: Algorithm optimizations and parallelization

---

## Quick Links

### External Resources
- **🐙 [GitHub Repository](https://github.com/ferdinand-popp/telomerehunter2)**: Source code and releases
- **🐳 [Docker Hub](https://hub.docker.com/r/fpopp22/telomerehunter2)**: Container images
- **📦 [PyPI Package](https://pypi.org/project/telomerehunter2/)**: Python package distribution
- **📄 [Original Paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0)**: Scientific publication

### Internal Navigation
- **[Installation Instructions](Installation-and-Usage.md)**: Get started quickly
- **[Input/Output Guide](Input-and-Output.md)**: Understand file formats
- **[Advanced Usage](Advanced-Options-and-Customization.md)**: Power user features
- **[Developer Guide](Contributing-and-Development.md)**: Contribute to the project
- **[Citations](References-and-Citations.md)**: Academic references

---

*This Wiki is actively maintained and updated. For the latest information, please check the GitHub repository and official documentation.*

**Last updated**: 2024  
**Maintainers**: Ferdinand Popp, Lars Feuerbach (DKFZ Applied Bioinformatics)