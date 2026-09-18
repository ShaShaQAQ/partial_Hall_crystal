from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / "results" / "phc30_np12_v1_100_cdw"


def load_numeric_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        rows.append(stripped.split())
    return rows


class CdwResultPackageTests(unittest.TestCase):
    def test_case3_plot_marks_only_the_three_lowest_levels(self):
        case_dir = ROOT / "ed" / "cases" / "case3_30sites_Np12"
        merge_source = (case_dir / "merge.jl").read_text()
        plot_source = (case_dir / "plot.jl").read_text()

        self.assertIn("n_gs = 3", merge_source)
        self.assertNotIn("n_gs = 15", merge_source)
        self.assertIn("sortperm(es)[1:min(n_gs, length(es))]", plot_source)
        self.assertNotIn("gs1_set = Set(v[1]", plot_source)

    def test_spectrum_identifies_ground_state_triplet(self):
        spectrum = RESULT / "data" / "spectrum.dat"
        self.assertTrue(spectrum.is_file())
        rows = load_numeric_rows(spectrum)
        levels = sorted((float(row[1]), int(row[0])) for row in rows)

        self.assertEqual({sector for _, sector in levels[:3]}, {0, 5, 10})
        self.assertAlmostEqual(levels[2][0] - levels[0][0], 0.0001245629)
        self.assertAlmostEqual(levels[3][0] - levels[2][0], 0.0212691488)

    def test_structure_factor_contains_three_states_and_average(self):
        structure_factor = (
            RESULT / "data" / "structure_factor_ground_states.dat"
        )
        self.assertTrue(structure_factor.is_file())
        rows = load_numeric_rows(structure_factor)
        by_state = {}
        for row in rows:
            state, momentum = row[0], int(row[1])
            by_state.setdefault(state, {})[momentum] = float(row[2])

        self.assertEqual(set(by_state), {"0", "5", "10", "average"})
        for values in by_state.values():
            self.assertEqual(set(values), set(range(15)))
            peak_momenta = sorted(values, key=values.get, reverse=True)[:2]
            self.assertEqual(set(peak_momenta), {5, 10})
            self.assertGreater(values[5], 1.4)
            self.assertGreater(values[10], 1.4)

    def test_parameters_and_wavefunction_manifest_are_traceable(self):
        parameters_path = RESULT / "parameters.toml"
        manifest_path = RESULT / "manifest.toml"
        self.assertTrue(parameters_path.is_file())
        self.assertTrue(manifest_path.is_file())
        parameters = parameters_path.read_text()
        for assignment in (
            "Ns = 30", "Nuc = 15", "Np = 12", "t1 = 1.0",
            "t3 = 0.2", "V1 = 100.0", "V2 = 0.0", "V3 = 0.0",
        ):
            self.assertIn(assignment, parameters)

        manifest = manifest_path.read_text()
        self.assertEqual(manifest.count("[[wavefunctions]]"), 3)
        self.assertEqual(
            set(re.findall(r"sector\s*=\s*(\d+)", manifest)),
            {"0", "5", "10"},
        )
        hashes = re.findall(r'sha256\s*=\s*"([0-9a-f]{64})"', manifest)
        sizes = [int(value) for value in re.findall(
            r"size_bytes\s*=\s*(\d+)", manifest
        )]
        self.assertEqual(len(hashes), 3)
        self.assertEqual(len(sizes), 3)
        self.assertTrue(all(size > 400_000_000 for size in sizes))


if __name__ == "__main__":
    unittest.main()
