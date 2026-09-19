from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / "results" / "phc30_np12_v1_10_v2_2_v3_2_fqahc"


class ModerateFqahcResultPackageTests(unittest.TestCase):
    def require_file(self, relative_path):
        path = RESULT / relative_path
        self.assertTrue(path.is_file(), str(path))
        return path

    def test_result_package_has_required_small_files(self):
        required = (
            RESULT / "parameters.toml",
            RESULT / "manifest.toml",
            RESULT / "data" / "spectrum.dat",
            RESULT / "data" / "manifold_diagnostics.txt",
            RESULT / "data" / "structure_factor_ground_manifold.dat",
        )
        for path in required:
            self.assertTrue(path.is_file(), str(path))
            self.assertGreater(path.stat().st_size, 0)
        self.assertEqual(list(RESULT.rglob("*.jld2")), [])

    def test_spectrum_has_isolated_fifteen_state_manifold(self):
        path = self.require_file(Path("data") / "spectrum.dat")
        rows = np.loadtxt(path, comments="#")
        self.assertEqual(rows.shape, (120, 2))
        sectors = rows[:, 0].astype(int)
        self.assertEqual(set(sectors), set(range(15)))
        self.assertTrue(all(np.count_nonzero(sectors == k) == 8 for k in range(15)))

        order = np.argsort(rows[:, 1], kind="stable")
        lowest = rows[order[:15]]
        self.assertEqual(set(lowest[:, 0].astype(int)), set(range(15)))
        width = lowest[-1, 1] - lowest[0, 1]
        separation = rows[order[15], 1] - lowest[-1, 1]
        self.assertLess(width, separation)

    def test_structure_factor_has_period_three_peaks(self):
        path = self.require_file(
            Path("data") / "structure_factor_ground_manifold.dat"
        )
        rows = [
            line.split()
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        self.assertEqual(len(rows), 240)
        labels = {row[0] for row in rows}
        self.assertEqual(labels, {str(k) for k in range(15)} | {"average"})
        average = np.asarray(
            [[float(value) for value in row[1:]] for row in rows if row[0] == "average"]
        )
        self.assertEqual(average.shape, (15, 4))
        peak_q = set(average[np.argsort(average[:, 1])[-2:], 0].astype(int))
        self.assertEqual(peak_q, {5, 10})

    def test_manifest_records_all_remote_wavefunctions(self):
        manifest = self.require_file("manifest.toml").read_text()
        self.assertIn('storage_host = "W003"', manifest)
        self.assertEqual(manifest.count("[[partial_files]]"), 15)
        self.assertEqual(manifest.count("sha256 = "), 15)
        parameters = self.require_file("parameters.toml").read_text()
        for entry in ("Np = 12", "V1 = 10.0", "V2 = 2.0", "V3 = 2.0"):
            self.assertIn(entry, parameters)


if __name__ == "__main__":
    unittest.main()
