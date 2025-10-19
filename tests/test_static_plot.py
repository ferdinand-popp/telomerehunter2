import os
import tempfile
import unittest

import plotly.express as px


class TestPlotlyKaleidoStaticExport(unittest.TestCase):
    def test_static_export(self):
        fig = px.scatter(x=[1, 2, 3], y=[4, 5, 6], title="Kaleido Static Export Test")
        formats = ["png", "pdf", "svg"]
        with tempfile.TemporaryDirectory() as tmpdir:
            for ext in formats:
                out_file = os.path.join(tmpdir, f"test_plot.{ext}")
                fig.write_image(out_file)
                self.assertTrue(os.path.exists(out_file), f"{ext} file was not created")
                self.assertGreater(os.path.getsize(out_file), 0, f"{ext} file is empty")
                os.remove(out_file)
