import os
import tempfile
import unittest

from telomerehunter2.estimate_telomere_content import get_telomere_counts

HEADER = "chr\tband\treads_with_pattern\tTTAGGG\tother\n"


class TestGetTelomereCounts(unittest.TestCase):
    """get_telomere_counts must find the 'unmapped' row by name, not by position."""

    def _counts(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, "P_spectrum.tsv"), "w") as f:
                f.write(HEADER + "".join(rows))
            return get_telomere_counts(tmp, "P")

    def test_unmapped_row_last(self):
        # banded layout: unmapped is always written last
        counts = self._counts(
            [
                "1\tp36.33\t7\t0\t0\n",
                "2\tp25.3\t3\t0\t0\n",
                "unmapped\tunmapped\t100\t0\t0\n",
            ]
        )
        self.assertEqual(counts, (110, 100))

    def test_unmapped_row_first(self):
        # no-banding layout: rows follow dict insertion order (issue #27)
        counts = self._counts(
            [
                "unmapped\tunmapped\t100\t0\t0\n",
                "1\t_\t7\t0\t0\n",
                "chrUn_x\t_\t3\t0\t0\n",
            ]
        )
        self.assertEqual(counts, (110, 100))

    def test_unmapped_row_in_middle(self):
        counts = self._counts(
            [
                "1\t_\t7\t0\t0\n",
                "unmapped\tunmapped\t100\t0\t0\n",
                "chrUn_x\t_\t3\t0\t0\n",
            ]
        )
        self.assertEqual(counts, (110, 100))

    def test_no_unmapped_row(self):
        # must not fall back to the last chromosome's count
        counts = self._counts(["1\t_\t7\t0\t0\n", "chrUn_x\t_\t3\t0\t0\n"])
        self.assertEqual(counts, (10, 0))

    def test_only_unmapped_row(self):
        self.assertEqual(self._counts(["unmapped\tunmapped\t5\t0\t0\n"]), (5, 5))

    def test_header_only(self):
        self.assertEqual(self._counts([]), (0, 0))

    def test_missing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(get_telomere_counts(tmp, "P"), (0, 0))


if __name__ == "__main__":
    unittest.main()
