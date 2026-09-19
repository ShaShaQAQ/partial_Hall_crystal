from pathlib import Path
import re
import unittest


REPORT_DIR = Path(__file__).resolve().parent
SOURCE = REPORT_DIR / "dmrg_summary.tex"


class ReportSourceTests(unittest.TestCase):
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

    def test_moderate_cdw_results_are_explained_in_chinese(self):
        text = SOURCE.read_text()
        for phrase in (
            "温和参数的 30-site CDW 态",
            "E_{15}-E_1=0.034430068375",
            "E_{16}-E_{15}=0.560257695317",
            "15 态等权平均",
            "归类为 CDW",
            "不是 DMRG 动力学计算",
            "response-Lanczos",
        ):
            self.assertIn(phrase, text)
        for rejected in (
            "FQAHC 候选态",
            "五重分数拓扑简并",
            r"15=3\times5",
        ):
            self.assertNotIn(rejected, text)

    def test_moderate_cdw_figures_exist_and_are_referenced(self):
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
