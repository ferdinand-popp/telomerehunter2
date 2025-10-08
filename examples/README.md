# TelomereHunter2 Examples

This directory contains example scripts demonstrating various usage patterns for TelomereHunter2.

## Windows Multiprocessing Example

**File:** `windows_multiprocessing_example.py`

This example demonstrates the correct pattern for calling TelomereHunter2 functions that use multiprocessing, ensuring compatibility with Windows.

### Key Takeaways

1. **Always use `if __name__ == "__main__":` guard** when calling functions that use `ProcessPoolExecutor` or multiprocessing
2. **Why?** On Windows, Python's multiprocessing uses 'spawn' instead of 'fork', which re-imports the entire module in each subprocess
3. **Without the guard:** Each subprocess would try to start more subprocesses, causing infinite recursion and `BrokenProcessPool` errors

### Usage

```python
# ❌ WRONG - Would fail on Windows
from telomerehunter2.filter_telomere_reads import parallel_filter_telomere_reads
parallel_filter_telomere_reads(...)  # Not protected!

# ✅ CORRECT - Works on all platforms
if __name__ == "__main__":
    from telomerehunter2.filter_telomere_reads import parallel_filter_telomere_reads
    parallel_filter_telomere_reads(...)  # Protected by main guard
```

### For Production Use

The built-in CLI commands already include proper Windows compatibility:
- `telomerehunter2` - for bulk analysis
- `telomerehunter2-sc` - for single-cell analysis

These are the recommended way to use TelomereHunter2.

## Additional Resources

- See `README.md` in the project root for full documentation
- See `tests/test_windows_compatibility.py` for automated compatibility tests
- See the "Troubleshooting" section in the main README for more Windows-specific guidance
