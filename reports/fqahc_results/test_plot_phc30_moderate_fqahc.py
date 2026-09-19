from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


REPORT_DIR = Path(__file__).resolve().parent
PLOT_SOURCE = REPORT_DIR / "plot_phc30_moderate_fqahc.py"
sys.path.insert(0, str(REPORT_DIR))

from plot_phc30_moderate_fqahc import (  # noqa: E402
    generate_figures,
    load_optical_curve,
    load_structure_factor,
    summarize_manifold,
    write_optical_manifold_average,
)


class ModerateFqahcPlotTests(unittest.TestCase):
    def test_plot_labels_moderate_manifold_as_fqahc_candidate(self):
        text = PLOT_SOURCE.read_text()
        self.assertIn("3x5 候选低能流形放大", text)
        self.assertIn("moderate-coupling FQAHC candidate", text)
        self.assertNotIn("CDW 低能流形放大", text)

    def write_synthetic_result(self, directory):
        spectrum_path = directory / "spectrum.dat"
        spectrum_lines = ["# k E-E0"]
        for sector in range(15):
            spectrum_lines.append(f"{sector} {0.002 * sector}")
            for level in range(1, 8):
                spectrum_lines.append(
                    f"{sector} {0.8 + 0.1 * level + 0.001 * sector}"
                )
        spectrum_path.write_text("\n".join(spectrum_lines) + "\n")

        structure_path = directory / "structure.dat"
        structure_lines = ["# state_sector q N(q) qx qy"]
        for sector in range(15):
            for q in range(15):
                peak = 1.0 if q in (5, 10) else 0.05
                structure_lines.append(
                    f"{sector} {q} {peak + sector / 1000} {q} 0"
                )
        for q in range(15):
            peak = 1.007 if q in (5, 10) else 0.057
            structure_lines.append(f"average {q} {peak} {q} 0")
        structure_path.write_text("\n".join(structure_lines) + "\n")

        optical_paths = {}
        for sector in range(15):
            path = directory / f"sector_{sector}_optical_response.dat"
            lines = ["# omega Re_total Im_total Re_regular Im_regular Re_drude Im_drude"]
            for index in range(101):
                omega = index / 20
                regular_real = (1 + sector / 100) / (1 + (omega - 2) ** 2)
                regular_imag = (omega - 2) * regular_real
                lines.append(
                    f"{omega} {regular_real - 0.1} {regular_imag} "
                    f"{regular_real} {regular_imag} -0.1 0"
                )
            path.write_text("\n".join(lines) + "\n")
            optical_paths[sector] = path
        return spectrum_path, structure_path, optical_paths

    def test_summarizes_lowest_fifteen_levels(self):
        rows = []
        for sector in range(15):
            rows.append((sector, 0.002 * sector))
            for level in range(1, 8):
                rows.append((sector, 0.8 + 0.1 * level + 0.001 * sector))

        summary = summarize_manifold(np.asarray(rows), manifold_size=15)

        self.assertEqual(summary["sectors"], tuple(range(15)))
        self.assertAlmostEqual(summary["width"], 0.028)
        self.assertAlmostEqual(summary["next_energy"], 0.9)
        self.assertAlmostEqual(summary["separation"], 0.872)
        self.assertTrue(summary["isolated"])

    def test_loads_fifteen_structure_factor_curves_and_average(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "structure.dat"
            lines = ["# state_sector q N(q) qx qy"]
            for sector in range(15):
                for q in range(15):
                    lines.append(f"{sector} {q} {sector + q / 10} {q} 0")
            for q in range(15):
                lines.append(f"average {q} {q / 10} {q} 0")
            path.write_text("\n".join(lines) + "\n")

            curves, average = load_structure_factor(path)

        self.assertEqual(set(curves), set(range(15)))
        self.assertTrue(all(curve.shape == (15, 4) for curve in curves.values()))
        self.assertEqual(average.shape, (15, 4))
        np.testing.assert_allclose(average[:, 0], np.arange(15))

    def test_loads_seven_column_optical_curve(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sector_5_optical_response.dat"
            path.write_text(
                "# omega Re_total Im_total Re_regular Im_regular "
                "Re_drude Im_drude\n"
                "0.0 1 2 3 4 5 6\n"
                "0.1 2 3 4 5 6 7\n"
            )

            curve = load_optical_curve(path)

        self.assertEqual(curve.shape, (2, 7))
        np.testing.assert_allclose(curve[:, 0], [0.0, 0.1])

    def test_generates_spectrum_structure_and_four_optical_figures(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            spectrum, structure, optical = self.write_synthetic_result(directory)
            output_dir = directory / "figures"

            outputs = generate_figures(
                spectrum, structure, optical, output_dir
            )

            self.assertEqual(
                set(outputs),
                {
                    "spectrum",
                    "structure_factor",
                    "regular_real",
                    "regular_imag",
                    "total_real",
                    "total_imag",
                },
            )
            self.assertTrue(all(path.is_file() for path in outputs.values()))
            self.assertTrue(all(path.stat().st_size > 10_000 for path in outputs.values()))

    def test_writes_fifteen_state_optical_statistics(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            _, _, optical = self.write_synthetic_result(directory)
            output_path = directory / "optical_manifold_average.dat"

            statistics = write_optical_manifold_average(optical, output_path)

            self.assertEqual(statistics.shape, (101, 25))
            self.assertTrue(output_path.is_file())
            curves = np.stack([load_optical_curve(optical[k]) for k in range(15)])
            np.testing.assert_allclose(statistics[:, 1:7], curves[:, :, 1:7].mean(axis=0))
            np.testing.assert_allclose(statistics[:, 7:13], curves[:, :, 1:7].std(axis=0))
            np.testing.assert_allclose(statistics[:, 13:19], curves[:, :, 1:7].min(axis=0))
            np.testing.assert_allclose(statistics[:, 19:25], curves[:, :, 1:7].max(axis=0))


if __name__ == "__main__":
    unittest.main()
