# Contributing and Development

This page provides comprehensive guidelines for contributing to TelomereHunter2, setting up a development environment, and understanding the project structure.

## Contributing Guidelines

We welcome contributions to TelomereHunter2! Whether you're fixing bugs, adding features, improving documentation, or suggesting enhancements, your contributions are valuable.

### Types of Contributions

- **Bug Reports**: Report issues or unexpected behavior
- **Feature Requests**: Suggest new functionality or improvements
- **Code Contributions**: Submit bug fixes or new features
- **Documentation**: Improve or expand documentation
- **Testing**: Add test cases or improve test coverage

### Before Contributing

1. **Check existing issues**: Look for similar bug reports or feature requests
2. **Discuss major changes**: Open an issue to discuss significant modifications
3. **Read the code**: Familiarize yourself with the codebase structure
4. **Review coding standards**: Follow the project's coding conventions

## Development Workflow

### 1. Fork and Clone

```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/telomerehunter2.git
cd telomerehunter2

# Add upstream remote
git remote add upstream https://github.com/ferdinand-popp/telomerehunter2.git
```

### 2. Create a Development Environment

#### Using Virtual Environment (Recommended)
```bash
# Create and activate virtual environment
python -m venv dev_env
source dev_env/bin/activate  # Linux/macOS
# dev_env\Scripts\activate   # Windows

# Install in development mode
pip install -e . --no-cache-dir

# Install development dependencies
pip install pytest flake8 black isort
```

#### Using Poetry
```bash
# Install Poetry if not already available
curl -sSL https://install.python-poetry.org | python3 -

# Install dependencies
poetry env use python3
poetry install

# Activate environment
poetry shell
```

#### Using Conda
```bash
# Create conda environment
conda create -n telomerehunter2-dev python=3.9
conda activate telomerehunter2-dev

# Install dependencies
pip install -e . --no-cache-dir
pip install pytest flake8 black isort
```

### 3. Create Feature Branch

```bash
# Create and switch to feature branch
git checkout -b feature/your-feature-name
# or for bug fixes:
git checkout -b fix/issue-description
```

### 4. Make Changes

Follow these guidelines when making changes:

- **Write clear commit messages** describing your changes
- **Keep commits focused** - one logical change per commit
- **Update documentation** if your changes affect user-facing functionality
- **Add tests** for new features or bug fixes
- **Follow code style** guidelines (see below)

### 5. Test Your Changes

```bash
# Run existing tests
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_telomerehunter2.py -v

# Test installation
telomerehunter2 --help
telomerehunter2-sc --help

# Manual testing with example data
telomerehunter2 -ibt tests/example.bam -o test_output/ -p TestRun \
  -b src/telomerehunter2/cytoband_files/hg19_cytoBand.txt
```

### 6. Code Quality Checks

#### Linting and Formatting
```bash
# Check code style with flake8
flake8 src/

# Format code with black
black src/

# Sort imports with isort
isort src/

# Run all quality checks
flake8 src/ && black --check src/ && isort --check-only src/
```

#### Using tox (if configured)
```bash
# Run all tox environments
tox

# Run specific environment
tox -e py39
tox -e lint
tox -e black
```

### 7. Submit Pull Request

```bash
# Push your branch
git push origin feature/your-feature-name

# Create pull request on GitHub
# Include:
# - Clear description of changes
# - Reference to related issues
# - Screenshots for UI changes
# - Test results
```

## Project Structure

### Directory Layout

```
telomerehunter2/
├── src/telomerehunter2/           # Main package source
│   ├── __init__.py               # Package initialization
│   ├── telomerehunter2_main.py   # Main entry point
│   ├── telomerehunter2_sc.py     # Single-cell entry point
│   ├── TVR_context.py            # TVR context analysis
│   ├── TVR_screen.py             # TVR screening
│   ├── estimate_telomere_content.py  # Telomere content estimation
│   ├── filter_telomere_reads.py  # Read filtering
│   ├── plot_functions.py         # Visualization functions
│   ├── utils.py                  # Utility functions
│   └── cytoband_files/           # Reference cytoband files
├── tests/                        # Test files
│   ├── test_telomerehunter2.py   # Main test suite
│   ├── test_telomerehunter2_sc.py # Single-cell tests
│   └── case*_expected_result.tsv # Expected test outputs
├── wiki/                         # Documentation (Wiki pages)
├── pyproject.toml               # Project configuration
├── README.md                    # Project overview
├── LICENSE.txt                  # License information
├── Dockerfile                   # Docker configuration
└── Apptainer_TH2.def           # Apptainer/Singularity definition
```

### Key Modules

#### Core Analysis Modules

1. **`telomerehunter2_main.py`**: Main analysis pipeline entry point
2. **`filter_telomere_reads.py`**: BAM file processing and read filtering
3. **`estimate_telomere_content.py`**: Telomere content calculation
4. **`TVR_context.py`**: Telomeric variant repeat context analysis
5. **`plot_functions.py`**: Visualization and plotting functions

#### Supporting Modules

1. **`utils.py`**: Shared utility functions
2. **`sort_telomere_reads.py`**: Read sorting and organization
3. **`get_repeat_threshold.py`**: Repeat threshold calculation
4. **`normalize_TVR_counts.py`**: TVR normalization functions

### Configuration Files

- **`pyproject.toml`**: Project metadata, dependencies, and tool configuration
- **`requirements.txt`**: Runtime dependencies (if present)
- **`setup.py`**: Legacy setup configuration (if present)

## Coding Standards

### Python Style Guidelines

Follow [PEP 8](https://pep8.org/) Python style guidelines:

```python
# Good: Clear function names and docstrings
def calculate_telomere_content(intratelomeric_reads, total_reads, gc_correction=True):
    """
    Calculate GC-corrected telomere content.
    
    Args:
        intratelomeric_reads (int): Number of intratelomeric reads
        total_reads (int): Total number of reads
        gc_correction (bool): Apply GC content correction
    
    Returns:
        float: Telomere content estimate
    """
    if total_reads == 0:
        return 0.0
    
    raw_content = (intratelomeric_reads / total_reads) * 1000000
    
    if gc_correction:
        return apply_gc_correction(raw_content)
    return raw_content

# Bad: Unclear names and no documentation
def calc_tc(ir, tr, gc=True):
    if tr == 0:
        return 0.0
    rc = (ir / tr) * 1000000
    if gc:
        return apply_gc_correction(rc)
    return rc
```

### Documentation Standards

#### Docstring Format
Use Google-style docstrings:

```python
def process_bam_file(bam_path, output_dir, sample_id):
    """Process BAM file for telomere analysis.
    
    This function filters reads, identifies telomeric sequences,
    and prepares data for content estimation.
    
    Args:
        bam_path (str): Path to input BAM file
        output_dir (str): Directory for output files
        sample_id (str): Sample identifier
    
    Returns:
        dict: Processing statistics including read counts
        
    Raises:
        FileNotFoundError: If BAM file doesn't exist
        ValueError: If sample_id is invalid
        
    Example:
        >>> stats = process_bam_file('sample.bam', 'output/', 'S1')
        >>> print(stats['total_reads'])
        1000000
    """
```

#### Comment Guidelines
```python
# Use comments to explain WHY, not WHAT
def apply_gc_correction(raw_content, gc_bins):
    # GC correction is necessary because telomeric sequences
    # have different GC content than the genome average
    correction_factor = calculate_gc_bias(gc_bins)
    return raw_content * correction_factor
```

### Error Handling

Implement proper error handling and logging:

```python
import logging

logger = logging.getLogger(__name__)

def safe_file_processing(file_path):
    """Process file with proper error handling."""
    try:
        with open(file_path, 'r') as f:
            data = f.read()
            return process_data(data)
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except PermissionError:
        logger.error(f"Permission denied: {file_path}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error processing {file_path}: {e}")
        raise
```

## Testing Guidelines

### Test Structure

Tests are located in the `tests/` directory and use pytest framework:

```python
import unittest
import tempfile
from pathlib import Path
from telomerehunter2 import estimate_telomere_content

class TestTelomereContent(unittest.TestCase):
    """Test telomere content estimation functions."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.sample_data = self.create_test_data()
    
    def test_content_calculation(self):
        """Test basic telomere content calculation."""
        result = estimate_telomere_content.calculate_content(
            intratelomeric_reads=100,
            total_reads=1000000
        )
        self.assertAlmostEqual(result, 100.0, places=2)
    
    def test_invalid_input(self):
        """Test error handling for invalid input."""
        with self.assertRaises(ValueError):
            estimate_telomere_content.calculate_content(
                intratelomeric_reads=-1,
                total_reads=1000000
            )
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        shutil.rmtree(self.test_dir)
```

### Writing Good Tests

1. **Test one thing at a time**: Each test should verify a single behavior
2. **Use descriptive names**: Test names should clearly indicate what is being tested
3. **Include edge cases**: Test boundary conditions and error scenarios
4. **Use fixtures**: Set up test data consistently
5. **Clean up**: Ensure tests don't leave artifacts

### Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run tests with coverage
python -m pytest tests/ --cov=telomerehunter2 --cov-report=html

# Run specific test class
python -m pytest tests/test_telomerehunter2.py::TestTelomereHunter2 -v

# Run tests matching pattern
python -m pytest tests/ -k "test_content" -v
```

## Adding New Features

### Feature Development Checklist

1. **Plan the feature**: Write a clear specification
2. **Design the API**: Plan function signatures and interfaces
3. **Write tests first**: Test-driven development when possible
4. **Implement incrementally**: Small, testable changes
5. **Update documentation**: Keep docs current with code changes
6. **Test thoroughly**: Unit tests, integration tests, manual testing

### Example: Adding a New Analysis Function

```python
# 1. Add the function to appropriate module
def analyze_telomere_variants(bam_file, variant_sequences):
    """Analyze custom telomere variant sequences.
    
    Args:
        bam_file (str): Path to BAM file
        variant_sequences (list): List of variant sequences to analyze
        
    Returns:
        dict: Variant analysis results
    """
    # Implementation here
    pass

# 2. Add corresponding test
def test_analyze_telomere_variants(self):
    """Test custom variant analysis."""
    variants = ['TTAGGG', 'TCAGGG']
    results = analyze_telomere_variants(self.test_bam, variants)
    self.assertIn('TTAGGG', results)
    self.assertIn('TCAGGG', results)

# 3. Update command-line interface if needed
parser.add_argument('--custom-variants', nargs='+',
                   help='Custom variant sequences to analyze')

# 4. Update documentation
# Add to appropriate Wiki page or README
```

## Example Data and Testing

### Test Data Location

Test data is located in the `tests/` directory:

- `test_telomerehunter2.py` - Main test suite
- `test_telomerehunter2_sc.py` - Single-cell specific tests
- `case*_expected_result.tsv` - Expected output files for validation

### Creating Test BAM Files

For testing, the test suite creates synthetic BAM files:

```python
@classmethod
def create_test_bam(cls, output_dir):
    """Create synthetic BAM file for testing."""
    # Implementation in test files
    # Creates minimal BAM with telomeric reads
    pass
```

### Validating Results

Tests compare actual outputs with expected results:

```python
def validate_output(self, result_file, expected_file):
    """Compare test results with expected output."""
    result_df = pd.read_csv(result_file, sep='\t')
    expected_df = pd.read_csv(expected_file, sep='\t')
    
    # Compare key columns
    pd.testing.assert_frame_equal(
        result_df[['tel_content', 'TRPM']], 
        expected_df[['tel_content', 'TRPM']],
        rtol=1e-3
    )
```

## Release Process

### Version Management

TelomereHunter2 uses semantic versioning (MAJOR.MINOR.PATCH):

- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality, backward compatible
- **PATCH**: Bug fixes, backward compatible

### Release Checklist

1. **Update version** in `pyproject.toml`
2. **Update CHANGELOG** with new features and fixes
3. **Run full test suite** on multiple Python versions
4. **Update documentation** as needed
5. **Create release tag** and push to repository
6. **Build and upload** to PyPI (maintainers only)
7. **Update Docker images** (maintainers only)

## Getting Help

### Communication Channels

- **GitHub Issues**: For bug reports and feature requests
- **GitHub Discussions**: For general questions and discussions
- **Email**: f.popp@dkfz-heidelberg.de (maintainer contact)

### Reporting Bugs

When reporting bugs, include:

1. **TelomereHunter2 version**: `telomerehunter2 --version`
2. **Python version**: `python --version`
3. **Operating system**: Linux/macOS/Windows
4. **Complete command**: Command that caused the error
5. **Error message**: Full error output
6. **Input description**: BAM file size, genome reference, etc.
7. **Expected behavior**: What you expected to happen

### Feature Requests

When requesting features:

1. **Clear description**: What functionality do you need?
2. **Use case**: Why is this feature important?
3. **Proposed implementation**: Any ideas on how it could work?
4. **Alternatives considered**: Other approaches you've thought about

## Maintainer Information

### Current Maintainers

- **Ferdinand Popp** (f.popp@dkfz-heidelberg.de) - Lead developer
- **Lars Feuerbach** (l.feuerbach@dkfz-heidelberg.de) - Core developer

### Institutional Affiliation

Developed at the German Cancer Research Center (DKFZ) - Division Applied Bioinformatics

### Project History

TelomereHunter2 is the advanced Python implementation of the original TelomereHunter tool, with enhanced features for modern genomics workflows.

## License and Legal

### License Information

TelomereHunter2 is released under the GNU General Public License v3.0. See [LICENSE.txt](../LICENSE.txt) for full details.

### Contributor License Agreement

By contributing to TelomereHunter2, you agree that your contributions will be licensed under the same GNU GPL v3.0 license.

### Attribution

All contributors are recognized in the project documentation and releases. Significant contributions may warrant co-authorship in scientific publications.

## Further Reading

- [Installation and Usage](Installation-and-Usage.md) - Setup and basic usage
- [Input and Output](Input-and-Output.md) - File format documentation  
- [Advanced Options and Customization](Advanced-Options-and-Customization.md) - Advanced usage
- [References and Citations](References-and-Citations.md) - Academic references