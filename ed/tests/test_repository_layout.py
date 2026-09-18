from pathlib import Path
import re
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

    def test_shared_includes_resolve_inside_ed_shared(self):
        shared_root = (ROOT / "ed" / "shared").resolve()
        offenders = []
        for owner in (ROOT / "ed", ROOT / "dmrg"):
            for path in owner.rglob("*.jl"):
                for include_path in re.findall(
                    r'include\("([^"]*shared/[^"]+)"\)', path.read_text()
                ):
                    resolved = (path.parent / include_path).resolve()
                    if shared_root not in resolved.parents or not resolved.is_file():
                        offenders.append(
                            "{} -> {}".format(
                                path.relative_to(ROOT).as_posix(), include_path
                            )
                        )

        self.assertEqual(offenders, [])

    def test_literal_julia_includes_resolve(self):
        offenders = []
        for owner in (ROOT / "ed", ROOT / "dmrg", ROOT / "reports"):
            for path in owner.rglob("*.jl"):
                for include_path in re.findall(
                    r'include\("([^"]+)"\)', path.read_text()
                ):
                    resolved = (path.parent / include_path).resolve()
                    if not resolved.is_file():
                        offenders.append(
                            "{} -> {}".format(
                                path.relative_to(ROOT).as_posix(), include_path
                            )
                        )

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

    def test_w003_optical_job_uses_available_node_shape(self):
        optical_job = (
            ROOT / "ed" / "cases" / "case3_30sites_Np12" /
            "submit" / "run_optical_response_w003.pbs"
        )
        structure_job = (
            ROOT / "ed" / "cases" / "case3_30sites_Np12" /
            "submit" / "run_cdw_structure_factor_w003.pbs"
        )
        for job in (optical_job, structure_job):
            self.assertTrue(job.is_file())
            text = job.read_text()
            self.assertIn("#PBS -q short", text)
            self.assertIn("select=1:ncpus=24:mem=90gb", text)
            self.assertIn("--threads=24", text)
        self.assertIn("run_optical_response.jl", optical_job.read_text())
        self.assertIn(
            "compute_cdw_structure_factor.jl", structure_job.read_text()
        )


if __name__ == "__main__":
    unittest.main()
