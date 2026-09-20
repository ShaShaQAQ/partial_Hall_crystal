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

    def test_w003_moderate_fqahc_workflow_is_isolated(self):
        submit_dir = (
            ROOT / "ed" / "cases" / "case3_30sites_Np12" / "submit"
        )
        jobs = {
            "spectrum": submit_dir / "run_spectrum_w003.pbs",
            "analysis": submit_dir / "analyze_moderate_fqahc_w003.pbs",
            "optical": submit_dir / "run_moderate_optical_w003.pbs",
            "optical_all": submit_dir / "run_moderate_optical_all_w003.pbs",
        }
        for path in jobs.values():
            self.assertTrue(path.is_file())
            text = path.read_text()
            self.assertIn("#PBS -q short", text)
            self.assertIn("select=1:ncpus=24:mem=90gb", text)
            self.assertIn("phc30_np12_v1_10_v2_2_v3_2_fqahc", text)

        spectrum = jobs["spectrum"].read_text()
        self.assertIn("run_spectrum.jl", spectrum)
        self.assertIn("PBS_ARRAY_INDEX", spectrum)
        self.assertIn("PHC_SECTOR_START=${PBS_ARRAY_INDEX}", spectrum)
        self.assertIn("PHC_SECTOR_END=${PHC_SECTOR_START}", spectrum)
        self.assertNotIn("5 * PBS_ARRAY_INDEX", spectrum)
        self.assertIn("--lattice corrected", spectrum)
        self.assertIn("--V1 10.0", spectrum)
        self.assertIn("--V2 2.0", spectrum)
        self.assertIn("--V3 2.0", spectrum)

        analysis = jobs["analysis"].read_text()
        self.assertIn("analyze_ground_manifold.jl", analysis)

        optical = jobs["optical"].read_text()
        self.assertIn("run_optical_response.jl", optical)
        self.assertIn("PBS_ARRAY_INDEX", optical)
        self.assertIn("array_sectors=(5 10 0)", optical)
        self.assertIn("--lattice corrected", optical)
        self.assertIn("--V1 10.0", optical)
        self.assertIn("--V2 2.0", optical)
        self.assertIn("--V3 2.0", optical)

        optical_all = jobs["optical_all"].read_text()
        self.assertIn("run_optical_response.jl", optical_all)
        self.assertIn("PHC_SECTOR=${PBS_ARRAY_INDEX}", optical_all)
        self.assertNotIn("PHC_SECTOR=${PHC_SECTOR:-", optical_all)
        self.assertIn("--lattice corrected", optical_all)
        self.assertIn("--V1 10.0", optical_all)
        self.assertIn("--V2 2.0", optical_all)
        self.assertIn("--V3 2.0", optical_all)

    def test_w003_np13_workflow_is_parameterized_and_isolated(self):
        submit_dir = (
            ROOT / "ed" / "cases" / "case3_30sites_Np12" / "submit"
        )
        jobs = {
            "spectrum": submit_dir / "run_np13_spectrum_w003.pbs",
            "analysis": submit_dir / "run_np13_analyze_w003.pbs",
            "optical": submit_dir / "run_np13_optical_w003.pbs",
        }
        for path in jobs.values():
            self.assertTrue(path.is_file(), str(path))
            text = path.read_text()
            self.assertIn("#PBS -q short", text)
            self.assertIn("select=1:ncpus=24:mem=90gb", text)
            self.assertIn("--threads=24", text)
            self.assertIn("PHC_RESULT_ID", text)
            self.assertIn("PHC_V1", text)
            self.assertIn("PHC_V2", text)
            self.assertIn("PHC_V3", text)
            self.assertIn("--Np 13", text)
            self.assertIn("--lattice corrected", text)
            self.assertIn("PBS_JOBNAME", text)

        spectrum = jobs["spectrum"].read_text()
        self.assertIn("PBS_ARRAY_INDEX", spectrum)
        self.assertIn("--nev 8", spectrum)
        self.assertIn("--krylovdim 60", spectrum)
        self.assertIn("run_spectrum.jl", spectrum)
        self.assertIn("phc13s_spec", spectrum)
        self.assertIn("phc13m_spec", spectrum)

        analysis = jobs["analysis"].read_text()
        self.assertIn("PHC_MANIFOLD_SIZE", analysis)
        self.assertIn("PHC_PHASE_LABEL", analysis)
        self.assertIn("analyze_ground_manifold.jl", analysis)
        self.assertIn("phc13s_m", analysis)
        self.assertIn("phc13m_m", analysis)

        optical = jobs["optical"].read_text()
        self.assertIn("PHC_SECTORS", optical)
        self.assertIn("PBS_ARRAY_INDEX", optical)
        self.assertIn("--mmax 600", optical)
        self.assertIn("--eta 0.065", optical)
        self.assertIn("run_optical_response.jl", optical)
        self.assertIn("phc13s_opt", optical)
        self.assertIn("phc13m_opt", optical)
        self.assertIn("manifold_diagnostics.txt", optical)
        self.assertIn("candidate_sectors", optical)


if __name__ == "__main__":
    unittest.main()
