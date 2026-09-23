from pathlib import Path
import re
import subprocess
import unittest


REPORT_DIR = Path(__file__).resolve().parent
SOURCE = REPORT_DIR / "dmrg_summary.tex"


class ReportSourceTests(unittest.TestCase):
    def test_compiled_note_pdf_is_present(self):
        pdf = REPORT_DIR / "dmrg_summary.pdf"
        self.assertTrue(pdf.is_file())
        self.assertGreater(pdf.stat().st_size, 1_000_000)

    def test_compiled_note_pdf_is_not_ignored(self):
        repository = REPORT_DIR.parents[1]
        relative_pdf = (REPORT_DIR / "dmrg_summary.pdf").relative_to(repository)
        process = subprocess.run(
            ["git", "check-ignore", "--quiet", str(relative_pdf)],
            cwd=str(repository),
            check=False,
        )
        self.assertEqual(process.returncode, 1)

    def test_cdw_section_uses_three_ground_states(self):
        text = SOURCE.read_text()
        self.assertIn(r"\(k=0,5,10\) 的三重 CDW 基态", text)
        self.assertIn(r"E_3-E_1=0.0001245629", text)
        self.assertIn(r"E_4-E_3=0.0212691488", text)
        self.assertNotIn("最低流形实际包含 15 个近简并态", text)
        self.assertNotIn("30-site FQAHC/PHC", text)

    def test_cdw_figures_exist_and_are_referenced(self):
        text = SOURCE.read_text()
        figures = (
            "phc30_manybody_spectrum.png",
            "phc30_cdw_structure_factor.png",
        )
        for name in figures:
            path = REPORT_DIR / "figures" / "optical_response" / name
            self.assertTrue(path.is_file())
            self.assertGreater(path.stat().st_size, 10_000)
            self.assertIn(name, text)

    def test_moderate_fqahc_candidate_is_explained_in_chinese(self):
        text = SOURCE.read_text()
        for phrase in (
            "温和参数的 30-site FQAHC 候选态",
            "E_{15}-E_1=0.034430068375",
            "E_{16}-E_{15}=0.560257695317",
            r"15=3\times5",
            "隔离的 15 态准简并低能流形",
            "五重拓扑部分",
            "flux insertion",
            "many-body Chern number",
            "15 态等权平均",
            "不是 DMRG 动力学计算",
            "response-Lanczos",
        ):
            self.assertIn(phrase, text)
        for rejected in (
            "温和参数的 30-site CDW 态",
            "本参数点归类为 CDW",
            "不再将该参数点解释为 FQAHC",
        ):
            self.assertNotIn(rejected, text)

    def test_moderate_fqahc_figures_exist_and_are_referenced(self):
        text = SOURCE.read_text()
        figures = (
            "phc30_moderate_manybody_spectrum.png",
            "phc30_moderate_structure_factor.png",
            "phc30_moderate_sigma_xx_regular_real.png",
            "phc30_moderate_sigma_xx_regular_imag.png",
            "phc30_moderate_sigma_xx_total_real.png",
            "phc30_moderate_sigma_xx_total_imag.png",
        )
        for name in figures:
            path = REPORT_DIR / "figures" / "optical_response" / name
            self.assertTrue(path.is_file(), str(path))
            self.assertGreater(path.stat().st_size, 10_000)
            self.assertIn(name, text)

    def test_np13_full_ed_parameters_and_spectrum_are_reported(self):
        text = SOURCE.read_text()
        for phrase in (
            "Np=13（13/30 filling）",
            "30-site 倾斜环面",
            "corrected lattice convention",
            "强耦合 Np=13",
            "E_0=892.7693923475454",
            "0.11898243171037848",
            "0.20949350939918077",
            "4.93366\\times10^{-12}",
            r"\overline{N(q=5)}=\overline{N(q=10)}=1.321613404399",
            "1.081443970873\\text{--}1.579102752688",
            "15 态流形",
            "扇区 0:14 各 1 态",
            "每扇区保留 8 能级",
            "候选流形",
        ):
            self.assertIn(phrase, text)

    def test_np13_optical_algorithms_and_numbers_are_reported(self):
        text = SOURCE.read_text()
        for phrase in (
            "response-Lanczos",
            "eta=0.065",
            "M=600",
            "15 态等权平均",
            "1.200",
            "3.773389760031857",
            "1.416",
            "2.05082509",
            "0.216",
            "0.902612709",
            "1.160--1.428",
            "3.44686--6.96221",
            "Drude",
            "-1.923004722924851",
            "0.8082451494173644",
            "-2.919037372848053",
            "-1.092715887328596",
            "不是 DMRG",
            "same-sector",
            "current-active",
            "有限 torus 固定 twist 曲率",
        ):
            self.assertIn(phrase, text)
        self.assertNotIn("shift current", text.lower())

    def test_np13_moderate_parameters_and_observables_are_reported(self):
        text = SOURCE.read_text()
        for phrase in (
            "温和参数 Np=13",
            "E_0=157.19949322891767",
            "0.1780098138966082",
            "0.25265895492933055",
            "6.36322\\times10^{-13}",
            r"\overline{N(q=5)}=\overline{N(q=10)}=1.076334465756",
            "0.819697952425",
            "1.376882017093",
            "q=6,9",
            "0.388442756221",
            "8.287294881631763\\times10^{-13}",
            "1.562",
            "8.69526118214673",
            "1.392",
            "6.16009960",
            "2.412",
            "0.715658572",
            "1.384--1.598",
            "8.93081--15.13362",
            "-4.38000670699518",
            "0.4861527438891664",
            "-4.888351330592577",
            "-3.370946817502301",
        ):
            self.assertIn(phrase, text)

    def test_np13_all_twelve_figures_are_referenced(self):
        text = SOURCE.read_text()
        for regime in ("strong", "moderate"):
            for suffix in (
                "manybody_spectrum.png",
                "structure_factor.png",
                "sigma_xx_regular_real.png",
                "sigma_xx_regular_imag.png",
                "sigma_xx_total_real.png",
                "sigma_xx_total_imag.png",
            ):
                name = f"phc30_np13_{regime}_{suffix}"
                path = REPORT_DIR / "figures" / "optical_response" / name
                self.assertTrue(path.is_file(), str(path))
                self.assertGreater(path.stat().st_size, 10_000)
                self.assertIn(name, text)

    def test_latex_environments_are_balanced(self):
        text = SOURCE.read_text()
        tokens = re.findall(r"\\(begin|end)\{([^}]+)\}", text)
        stack = []
        for action, environment in tokens:
            if action == "begin":
                stack.append(environment)
            else:
                self.assertTrue(stack, "end without begin: " + environment)
                self.assertEqual(stack.pop(), environment)
        self.assertEqual(stack, [])


if __name__ == "__main__":
    unittest.main()
