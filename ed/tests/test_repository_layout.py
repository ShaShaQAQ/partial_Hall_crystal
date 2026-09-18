from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]


class RepositoryLayoutTests(unittest.TestCase):
    def test_long_lived_modules_have_explicit_top_level_owners(self):
        expected = [
            ROOT / "ed" / "shared",
            ROOT / "ed" / "cases",
            ROOT / "ed" / "projected",
            ROOT / "dmrg",
            ROOT / "reports" / "fqahc_results",
            ROOT / "results",
        ]

        self.assertTrue(all(path.is_dir() for path in expected))
        self.assertFalse((ROOT / "shared").exists())
        self.assertFalse((ROOT / "projected_ed").exists())
        self.assertFalse(any(ROOT.glob("case[0-9]*")))
        self.assertFalse((ROOT / "dmrg" / "report").exists())

    def test_ed_sources_do_not_reference_legacy_root_shared_directory(self):
        offenders = []
        for path in (ROOT / "ed").rglob("*.jl"):
            if 'include("../shared/' in path.read_text():
                offenders.append(path.relative_to(ROOT).as_posix())

        self.assertEqual(offenders, [])

    def test_large_runtime_artifacts_are_ignored(self):
        probes = [
            "results/example/wavefunction.jld2",
            "ed/cases/example/output/run.out",
            "ed/cases/example/output/run.err",
        ]
        process = subprocess.run(
            ["git", "check-ignore", "--stdin"],
            cwd=ROOT,
            input="\n".join(probes) + "\n",
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

        self.assertEqual(process.returncode, 0)
        ignored = set(process.stdout.splitlines())
        self.assertEqual(ignored, set(probes))


if __name__ == "__main__":
    unittest.main()
