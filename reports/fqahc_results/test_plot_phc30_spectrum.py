from pathlib import Path
import sys
import unittest

import numpy as np


REPORT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(REPORT_DIR))

from plot_phc30_spectrum import (
    format_summary_text,
    load_spectrum,
    summarize_low_energy,
)


DATA_PATH = (
    REPORT_DIR.parents[1] / "results" / "phc30_np12_v1_100_cdw" /
    "data" / "spectrum.dat"
)


class SpectrumPlotTests(unittest.TestCase):
    def test_loads_eight_levels_in_each_of_fifteen_sectors(self):
        rows = load_spectrum(DATA_PATH)

        self.assertEqual(rows.shape, (120, 2))
        self.assertEqual(set(rows[:, 0].astype(int)), set(range(15)))
        self.assertTrue(all(
            np.count_nonzero(rows[:, 0] == k) == 8 for k in range(15)
        ))

    def test_summarizes_three_state_cdw_ground_manifold(self):
        summary = summarize_low_energy(load_spectrum(DATA_PATH))

        self.assertEqual(summary["manifold_size"], 3)
        self.assertEqual(summary["ground_sectors"], (0, 5, 10))
        self.assertAlmostEqual(summary["manifold_width"], 0.0001245629)
        self.assertAlmostEqual(summary["fourth_energy"], 0.0213937117)
        self.assertAlmostEqual(summary["separation"], 0.0212691488)

    def test_summary_annotation_uses_real_line_breaks(self):
        summary = summarize_low_energy(load_spectrum(DATA_PATH))

        text = format_summary_text(summary)

        self.assertNotIn("\\n", text)
        self.assertEqual(text.count("\n"), 2)
        self.assertIn("CDW ground-state triplet", text)


if __name__ == "__main__":
    unittest.main()
