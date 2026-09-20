from pathlib import Path
import sys
import tempfile
import unittest


REPORT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(REPORT_DIR))

from package_phc30_moderate_optics import package_optical_results  # noqa: E402


class ModerateOpticalPackagingTests(unittest.TestCase):
    def write_source_data(
        self, source_dir, sectors=range(15), nparticles=12,
        interactions=(10.0, 2.0, 2.0),
    ):
        v1, v2, v3 = interactions
        for sector in sectors:
            data_path = source_dir / f"sector_{sector}_optical_response.dat"
            data_path.write_text(
                "# sector={0} Eg=1.0 eta=0.065 area=25.0 Kexp=8.0 "
                "source_norm2=31.0 lattice=corrected Np={1} t1=1.0 t3=0.2 "
                "V1={2} V2={3} V3={4} requested_mmax=600 "
                "lanczos_steps=600 lanczos_breakdown=false residual=1e-12 "
                "drude_weight=-1.0\n"
                "# omega Re_total Im_total Re_regular Im_regular Re_drude Im_drude\n"
                "0.0 1 2 3 4 5 6\n"
                "0.1 2 3 4 5 6 7\n".format(
                    sector, nparticles, v1, v2, v3
                )
            )
            (source_dir / f"sector_{sector}_optical_response.jld2").write_bytes(
                ("synthetic-jld2-sector-%d" % sector).encode()
            )

    def test_packages_small_curves_and_hashes_remote_outputs(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            source_dir = directory / "source"
            result_dir = directory / "result"
            source_dir.mkdir()
            self.write_source_data(source_dir)

            summary = package_optical_results(source_dir, result_dir)

            optical_dir = result_dir / "data" / "optical"
            self.assertEqual(len(list(optical_dir.glob("sector_*.dat"))), 15)
            self.assertEqual(list(result_dir.rglob("*.jld2")), [])
            self.assertTrue((optical_dir / "optical_manifold_average.dat").is_file())
            manifest = (result_dir / "optical_manifest.toml").read_text()
            self.assertEqual(manifest.count("[[responses]]"), 15)
            self.assertEqual(manifest.count("data_sha256 = "), 15)
            self.assertEqual(manifest.count("jld2_sha256 = "), 15)
            self.assertEqual(summary["sectors"], tuple(range(15)))
            self.assertEqual(summary["maximum_residual"], 1e-12)

    def test_packages_selected_np13_manifold(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            source_dir = directory / "source"
            result_dir = directory / "result"
            source_dir.mkdir()
            sectors = (2, 7, 12)
            self.write_source_data(
                source_dir, sectors=sectors, nparticles=13,
                interactions=(100.0, 0.0, 0.0),
            )

            summary = package_optical_results(
                source_dir, result_dir,
                sectors=sectors,
                result_id="phc30_np13_v1_100_v2_0_v3_0",
                expected_metadata={
                    "Np": "13", "V1": "100.0",
                    "V2": "0.0", "V3": "0.0",
                },
            )

            optical_dir = result_dir / "data" / "optical"
            self.assertEqual(len(list(optical_dir.glob("sector_*.dat"))), 3)
            self.assertEqual(summary["sectors"], sectors)
            manifest = (result_dir / "optical_manifest.toml").read_text()
            self.assertIn(
                'result_id = "phc30_np13_v1_100_v2_0_v3_0"', manifest
            )
            self.assertEqual(manifest.count("[[responses]]"), 3)
            average = (optical_dir / "optical_manifold_average.dat").read_text()
            self.assertIn("3 态等权光电导统计量", average)


if __name__ == "__main__":
    unittest.main()
