import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pysam


class TestTelomereHunter2(unittest.TestCase):
    """Test cases for TelomereHunter2 package."""

    @classmethod
    def setUpClass(cls):
        """Set up test environment once for all test cases."""
        cls.project_root = Path(__file__).resolve().parent.parent
        cls.test_dir = tempfile.mkdtemp(prefix="telomerehunter2_test_")
        cls.test_dir = Path(cls.test_dir)

        # Create the banding file path
        cls.banding_file = (
            cls.project_root
            / "src"
            / "telomerehunter2"
            / "cytoband_files"
            / "hg19_cytoBand.txt"
        )

        # Check if the banding file exists
        if not cls.banding_file.exists():
            raise FileNotFoundError(f"Banding file not found at: {cls.banding_file}")

        # Create test data
        cls.bam_file_path = cls.create_hg19_bam(cls.test_dir)
        cls.cram_file_path = cls.test_dir / "test_control.cram"
        cls.convert_bam_to_cram(cls.bam_file_path, cls.banding_file, cls.cram_file_path)

        # Define expected result paths
        cls.case1_expected = cls.project_root / "tests" / "case1_expected_result.tsv"
        cls.case2_expected = cls.project_root / "tests" / "case2_expected_result.tsv"

    @classmethod
    def tearDownClass(cls):
        """Clean up after all tests."""
        # Remove the temporary test directory
        shutil.rmtree(cls.test_dir)

    @classmethod
    def create_hg19_bam(cls, output_dir):
        """Create an artificial BAM file for testing."""
        bam_path = os.path.join(output_dir, "test.bam")

        # Define BAM header with reference sequences
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [
                {"LN": 248956422, "SN": "chr1"},  # Example chromosome from hg19
                {"LN": 242193529, "SN": "chr2"},  # Example chromosome from hg19
            ],
        }

        # Create a new BAM file
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
            for i in range(1, 11):  # Generate 10 reads
                read = pysam.AlignedSegment()
                read.query_name = f"read_{i:03d}"
                if i in [1, 2]:
                    read.query_sequence = "TCGACTTTAGGGTTAGGGTTAGGGTGAGGGTTAGGGTTAGGGTTAGGGTCGACT"  # singleton
                    read.flag = 4  # Unmapped flag
                    read.reference_id = -1
                    read.reference_start = 3
                elif i in [3, 4]:  # Add an unmapped read (every 3rd read)
                    read.query_sequence = "TTAGGGTCAGGG" * 3
                    read.flag = 4  # Unmapped flag
                    read.reference_id = -1
                    read.reference_start = 3  # unmapped
                elif i in [5, 6]:  # Add a read with 7 * TTAGGG
                    read.query_sequence = "TTAGGG" * 7  # 7 repeats of TTAGGG
                    read.flag = 0  # No flags set
                    read.reference_id = 0
                    read.reference_start = 100 + i * 10
                else:  # Add a read with mixed sequence
                    read.query_sequence = "ATCGTTAGGG" * 12
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = 100 + i * 10

                read.mapping_quality = 60
                read.cigar = (
                    [(0, len(read.query_sequence))] if read.flag == 0 else []
                )  # Entire sequence is a match
                read.next_reference_id = -1  # Unpaired
                read.query_qualities = pysam.qualitystring_to_array(
                    "E" * len(read.query_sequence)
                )  # Quality scores
                bam.write(read)

            # Add reads aligning to hg19 genome
            for i in range(12, 122):  # Generate 100 additional reads
                read = pysam.AlignedSegment()
                read.query_name = f"hg19_read_{i:03d}"
                read.query_sequence = (
                    "CAGTTCAGTTCAGTTCAGTTCAGTT"  # Random sequence for hg19
                )
                read.flag = 0  # No flags set
                read.reference_id = (
                    0 if i % 2 == 0 else 1
                )  # Alternate between chr1 and chr2
                read.reference_start = i * 100
                read.mapping_quality = 60
                read.cigar = [
                    (0, len(read.query_sequence))
                ]  # Entire sequence is a match
                read.next_reference_id = -1  # Unpaired
                read.query_qualities = pysam.qualitystring_to_array(
                    "E" * len(read.query_sequence)
                )  # Quality scores
                bam.write(read)

        # Sort and index BAM file
        sorted_bam_path = bam_path.replace(".bam", ".sorted.bam")
        pysam.sort("-o", sorted_bam_path, bam_path)
        pysam.index(sorted_bam_path)

        return Path(sorted_bam_path)

    @classmethod
    def convert_bam_to_cram(cls, bam_path, reference_path, output_cram_path):
        """Convert a BAM file to a CRAM file using the specified reference file."""
        try:
            with pysam.AlignmentFile(
                str(bam_path), "rb"
            ) as bam_file, pysam.AlignmentFile(
                str(output_cram_path),
                "wc",
                header=bam_file.header,
                reference_filename=None,
            ) as cram_file:
                for read in bam_file:
                    cram_file.write(read)

            pysam.index(str(output_cram_path))
            return output_cram_path
        except Exception as e:
            raise RuntimeError(f"Error converting BAM to CRAM: {e}")

    def run_telomerehunter_package(
        self, bam_file_path, results_path, patient_name, bam_file_path_control=None
    ):
        """Run the telomerehunter2 source script via command-line interface."""
        script_path = (
            self.project_root / "src" / "telomerehunter2" / "telomerehunter2_main.py"
        )
        command = [
            sys.executable,
            str(script_path),
            "-ibt",
            str(bam_file_path),
            "-o",
            str(results_path),
            "-p",
            patient_name,
            "-b",
            str(self.banding_file),
        ]

        if bam_file_path_control:
            command += ["-ibc", str(bam_file_path_control)]

        print(f"Running command: {' '.join(command)}")

        proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        stdout, stderr = proc.communicate()  # Capture both stdout and stderr
        if proc.returncode != 0:
            raise RuntimeError(
                f"Command failed with exit code {proc.returncode}.\n"
                f"Stdout: {stdout}\nStderr: {stderr}"
            )
        print(stdout)

    def validate_results(self, results_path, expected_result_path):
        """Validate the results of a test run."""
        patient_name = "TEST_PATIENT"
        summary_file = Path(results_path) / patient_name / f"{patient_name}_summary.tsv"
        merged_plots_file = (
            Path(results_path) / patient_name / f"{patient_name}_all_plots_merged.pdf"
        )

        self.assertTrue(
            summary_file.exists(), f"Summary file not found: {summary_file}"
        )
        self.assertTrue(
            merged_plots_file.exists(),
            f"Merged plots file not found: {merged_plots_file}",
        )
        self.assertGreater(
            merged_plots_file.stat().st_size,
            0,
            f"Merged plots file is empty: {merged_plots_file}",
        )

        # Read summary TSV into a DataFrame
        df = pd.read_csv(summary_file, sep="\t")
        self.assertIn(
            "tel_content",
            df.columns,
            f"'tel_content' column not found in summary file: {summary_file}",
        )

        # Read expected result file
        if not expected_result_path.exists():
            # Create a copy of current results as expected for future tests
            shutil.copy(summary_file, expected_result_path)
            df_expected = df
        else:
            df_expected = pd.read_csv(expected_result_path, sep="\t")

        # Compare the actual and expected values
        pd.testing.assert_frame_equal(
            df,
            df_expected,
            check_dtype=False,
            check_index_type=False,
            check_column_type=False,
        )
        print(f"✅ Validation passed for {results_path}")

    def test_case1_basic_run(self):
        """Test Case 1: Basic run with simple BAM file."""
        results_path = self.test_dir / "case1"
        results_path.mkdir(exist_ok=True)

        self.run_telomerehunter_package(
            bam_file_path=self.bam_file_path,
            results_path=results_path,
            patient_name="TEST_PATIENT",
        )

        self.validate_results(results_path, self.case1_expected)

    def test_case2_combined_run(self):
        """Test Case 2: Combined run with BAM/CRAM file tumor control."""
        results_path = self.test_dir / "case2"
        results_path.mkdir(exist_ok=True)

        self.run_telomerehunter_package(
            bam_file_path=self.cram_file_path,
            results_path=results_path,
            patient_name="TEST_PATIENT",
            bam_file_path_control=self.bam_file_path,
        )

        self.validate_results(results_path, self.case2_expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
