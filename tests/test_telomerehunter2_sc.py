import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pysam


class TestSCTelomereHunter2(unittest.TestCase):
    """Test cases for single-cell TelomereHunter2 summary script."""

    @classmethod
    def setUpClass(cls):
        cls.project_root = Path(__file__).resolve().parent.parent
        cls.test_dir = tempfile.mkdtemp(prefix="sc_telomerehunter2_test_")
        cls.test_dir = Path(cls.test_dir)
        cls.banding_file = (
            cls.project_root
            / "src"
            / "telomerehunter2"
            / "cytoband_files"
            / "hg19_cytoBand.txt"
        )
        if not cls.banding_file.exists():
            raise FileNotFoundError(f"Banding file not found at: {cls.banding_file}")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.test_dir)

    @classmethod
    def create_sc_bam_with_barcodes(cls, output_dir):
        bam_path = os.path.join(output_dir, "sc_test.bam")
        header = {
            "HD": {"VN": "1.6"},
            "SQ": [
                {"LN": 248956422, "SN": "chr1"},
                {"LN": 242193529, "SN": "chr2"},
            ],
        }
        with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
            for i in range(1, 6):
                read = pysam.AlignedSegment()
                read.query_name = f"sc_read_{i:03d}"
                read.query_sequence = "TTAGGG" * (i + 1)
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 100 + i * 10
                read.mapping_quality = 60
                read.cigar = [(0, len(read.query_sequence))]
                read.next_reference_id = -1
                read.query_qualities = pysam.qualitystring_to_array(
                    "E" * len(read.query_sequence)
                )
                read.set_tag("CB", f"CELL{i}")
                bam.write(read)
        pysam.sort("-o", bam_path.replace(".bam", ".sorted.bam"), bam_path)
        pysam.index(bam_path.replace(".bam", ".sorted.bam"))
        return Path(bam_path.replace(".bam", ".sorted.bam"))

    def test_sc_summary_script(self):
        results_path = self.test_dir / "sc_case"
        results_path.mkdir(exist_ok=True)
        sc_bam_path = self.create_sc_bam_with_barcodes(results_path)
        patient_name = "SC_TEST_PATIENT"
        sc_script = (
            self.project_root / "src" / "telomerehunter2" / "telomerehunter2_sc.py"
        )
        command = [
            sys.executable,
            str(sc_script),
            "-ibt",
            str(sc_bam_path),
            "-o",
            str(results_path),
            "-p",
            patient_name,
            "-b",
            str(self.banding_file),
            "--min-reads-per-barcode",
            "1",
        ]
        proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        for line in proc.stdout:
            print(line, end="")
        exit_code = proc.wait()
        if exit_code != 0:
            error_message = proc.stderr.read()
            raise RuntimeError(
                f"SC summary script failed with exit code {exit_code}. Error: {error_message}"
            )
        summary_file = results_path / patient_name / f"{patient_name}_sc_summary.tsv"
        expected_file = self.project_root / "tests" / "sc_case_expected_result.tsv"
        self.assertTrue(
            summary_file.exists(), f"SC summary file not found: {summary_file}"
        )
        df = pd.read_csv(summary_file, sep="\t")
        self.assertIn(
            "tel_content",
            df.columns,
            f"'tel_content' column not found in SC summary file: {summary_file}",
        )

        # Check TVR column order: all TVRs arbitrary, then all TVRs 100bp, then all TVRs singletons
        TVR_HEXAMERS = [
            "TCAGGG",
            "TGAGGG",
            "TTGGGG",
            "TTCGGG",
            "TTTGGG",
            "ATAGGG",
            "CATGGG",
            "CTAGGG",
            "GTAGGG",
            "TAAGGG",
        ]
        tvr_arbitrary = [
            f"{tvr}_arbitrary_context_norm_by_intratel_reads" for tvr in TVR_HEXAMERS
        ]
        tvr_100bp = [
            f"{tvr}_arbitrary_context_per_100bp_intratel_read" for tvr in TVR_HEXAMERS
        ]
        tvr_singletons = [f"{tvr}_singletons_norm_by_all_reads" for tvr in TVR_HEXAMERS]
        std_columns = [
            "PID",
            "sample",
            "tel_content",
            "total_reads",
            "read_lengths",
            "repeat_threshold_set",
            "repeat_threshold_used",
            "intratelomeric_reads",
            "junctionspanning_reads",
            "subtelomeric_reads",
            "intrachromosomal_reads",
            "tel_read_count",
            "gc_bins_for_correction",
            "total_reads_with_tel_gc",
        ]
        expected_columns = std_columns + tvr_arbitrary + tvr_100bp + tvr_singletons
        self.assertEqual(
            list(df.columns),
            expected_columns,
            f"SC summary columns are not in the expected order.\nExpected: {expected_columns}\nFound: {list(df.columns)}",
        )

        if not expected_file.exists():
            shutil.copy(summary_file, expected_file)
            df_expected = df
        else:
            df_expected = pd.read_csv(expected_file, sep="\t")
        pd.testing.assert_frame_equal(
            df,
            df_expected,
            check_dtype=False,
            check_index_type=False,
            check_column_type=False,
        )
        print(f"✅ SC summary validation passed for {results_path}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
