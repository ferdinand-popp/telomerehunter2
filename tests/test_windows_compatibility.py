#!/usr/bin/env python3
"""
Test Windows multiprocessing compatibility.

This test verifies that the multiprocessing components used in TelomereHunter2
are compatible with Windows' 'spawn' method for process creation.
"""

import unittest
import pickle
import sys
from pathlib import Path

# Add src to path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "src"))

from telomerehunter2.utils import measure_time


# Define a sample function at module level for pickle testing
@measure_time
def sample_function_for_pickle_test(x):
    """Sample function for testing."""
    return x * 2


class TestWindowsCompatibility(unittest.TestCase):
    """Test Windows multiprocessing compatibility."""

    def test_measure_time_decorator_is_pickleable(self):
        """Test that the measure_time decorator is pickle-safe."""
        # Use the module-level function
        # This should not raise an exception
        pickled = pickle.dumps(sample_function_for_pickle_test)
        unpickled = pickle.loads(pickled)
        
        # Function should work after unpickling
        result = unpickled(5)
        self.assertEqual(result, 10)

    def test_process_region_is_pickleable(self):
        """Test that process_region worker function is pickle-safe."""
        from telomerehunter2.filter_telomere_reads import process_region
        
        # This should not raise an exception
        pickled = pickle.dumps(process_region)
        unpickled = pickle.loads(pickled)
        
        # The function should be callable (though we won't call it without args)
        self.assertTrue(callable(unpickled))

    def test_process_unmapped_reads_is_pickleable(self):
        """Test that process_unmapped_reads worker function is pickle-safe."""
        from telomerehunter2.filter_telomere_reads import process_unmapped_reads
        
        # This should not raise an exception
        pickled = pickle.dumps(process_unmapped_reads)
        unpickled = pickle.loads(pickled)
        
        # The function should be callable
        self.assertTrue(callable(unpickled))

    def test_filter_telomere_reads_has_main_guard(self):
        """Test that filter_telomere_reads.py has proper main guard."""
        filter_telomere_reads_path = (
            project_root / "src" / "telomerehunter2" / "filter_telomere_reads.py"
        )
        
        with open(filter_telomere_reads_path, "r") as f:
            content = f.read()
        
        # Check for main guard
        self.assertIn('if __name__ == "__main__":', content)
        
        # Check for Windows compatibility documentation
        self.assertIn("WINDOWS COMPATIBILITY", content)
        self.assertIn("spawn", content.lower())

    def test_all_entry_points_have_main_guard(self):
        """Test that all main entry points have proper main guard."""
        entry_points = [
            "telomerehunter2_main.py",
            "telomerehunter2_sc.py",
            "sc_barcode_splitter_run.py",
            "filter_telomere_reads.py"
        ]
        
        for entry_point in entry_points:
            entry_point_path = project_root / "src" / "telomerehunter2" / entry_point
            
            if not entry_point_path.exists():
                continue
            
            with open(entry_point_path, "r") as f:
                content = f.read()
            
            # Check for main guard
            self.assertIn(
                'if __name__ == "__main__":',
                content,
                f"{entry_point} should have main guard for Windows compatibility"
            )

    def test_measure_time_preserves_function_metadata(self):
        """Test that measure_time decorator preserves function metadata."""
        # Use the module-level function
        # functools.wraps should preserve __name__ and __doc__
        self.assertEqual(sample_function_for_pickle_test.__name__, "sample_function_for_pickle_test")
        self.assertEqual(sample_function_for_pickle_test.__doc__, "Sample function for testing.")


if __name__ == "__main__":
    unittest.main()
