#!/usr/bin/env python3
"""
Example script demonstrating Windows-compatible multiprocessing usage.

This script shows the correct pattern for calling TelomereHunter2 functions
that use multiprocessing, ensuring compatibility with Windows.
"""

def example_incorrect_usage():
    """
    ❌ INCORRECT: This would cause BrokenProcessPool errors on Windows!
    
    On Windows, when you run this script, it would be re-imported by each
    subprocess, causing infinite recursion of process creation.
    """
    pass
    # DON'T DO THIS:
    # from telomerehunter2.filter_telomere_reads import parallel_filter_telomere_reads
    # parallel_filter_telomere_reads(...)  # ❌ Not protected by main guard!


def example_correct_usage():
    """
    ✅ CORRECT: This is the proper way to use multiprocessing functions.
    
    The code is protected by `if __name__ == "__main__":`, which prevents
    subprocess re-execution on Windows.
    """
    from telomerehunter2.filter_telomere_reads import parallel_filter_telomere_reads
    
    # Example parameters
    bam_path = "sample.bam"
    out_dir = "./output"
    pid = "sample_id"
    
    # Call the multiprocessing function safely
    parallel_filter_telomere_reads(
        bam_path=bam_path,
        out_dir=out_dir,
        pid=pid,
        sample="test_sample",
        repeat_threshold_calc=4,
        mapq_threshold=0,
        repeats=["TTAGGG", "TCAGGG", "TGAGGG"],
        consecutive_flag=False,
        remove_duplicates=False,
        band_file=None,
        num_processes=None,
        singlecell_mode=False,
        fast_mode=False,
    )
    print(f"✅ Processing complete! Results in {out_dir}")


if __name__ == "__main__":
    # This guard is REQUIRED for Windows compatibility
    # 
    # Why? On Windows, Python's multiprocessing uses 'spawn' instead of 'fork',
    # which means each subprocess re-imports the entire script. Without this guard,
    # each subprocess would try to start more subprocesses infinitely!
    #
    # The built-in TelomereHunter2 CLI commands already include this protection.
    
    print("=" * 70)
    print("Windows-Compatible Multiprocessing Example")
    print("=" * 70)
    print()
    print("This script demonstrates the CORRECT way to use multiprocessing")
    print("functions from TelomereHunter2 in custom scripts.")
    print()
    print("Key points:")
    print("  1. Always protect multiprocessing code with if __name__ == '__main__':")
    print("  2. This prevents infinite subprocess creation on Windows")
    print("  3. The built-in CLI commands already include this protection")
    print()
    print("For production use, prefer the built-in commands:")
    print("  - telomerehunter2 (for bulk analysis)")
    print("  - telomerehunter2-sc (for single-cell analysis)")
    print()
    print("=" * 70)
    
    # Uncomment the next line to run the example (requires valid BAM file)
    # example_correct_usage()
