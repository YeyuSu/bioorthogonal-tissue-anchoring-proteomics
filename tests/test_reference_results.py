import json
import tempfile
import unittest
from pathlib import Path

from reproduce import EXPECTED_INPUT_SHA256, run_reproduction
from scripts.publication.build_publication_tables import CONTRAST_ORDER
from scripts.validate_release import (
    compare_reference_results,
    promote_reference_results,
    verify_reference_results,
)

ROOT = Path(__file__).resolve().parents[1]


class ReferenceResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.base = Path(cls.temp.name)
        cls.generated = run_reproduction(ROOT, cls.base / "generated")
        (cls.generated / "publication_figures" / "forbidden.png").write_bytes(b"png")
        (cls.generated / "publication_tables" / "forbidden.xlsx").write_bytes(b"xlsx")
        cls.repo = cls.base / "repo"
        cls.reference = cls.repo / "results" / "reference"
        cls.paths = promote_reference_results(cls.generated, cls.reference)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_promotion_copies_only_whitelisted_artifacts(self):
        copied = [path.resolve().relative_to(self.reference.resolve()).as_posix() for path in self.paths]
        self.assertTrue(copied)
        self.assertFalse(any(name.endswith((".png", ".pdf", ".docx", ".xlsx", ".pyc")) for name in copied))
        self.assertFalse(any("forbidden" in name for name in copied))
        self.assertEqual(len(list((self.reference / "figures").glob("Volcano_*.svg"))), 6)
        self.assertTrue((self.reference / "figures" / "three_tissue_treatment_effect_heatmap.svg").is_file())

    def test_validation_report_records_reproduction_contract(self):
        report = json.loads((self.reference / "validation_report.json").read_text(encoding="utf-8"))
        self.assertTrue(report["passed"])
        self.assertEqual(report["input_sha256"], EXPECTED_INPUT_SHA256)
        self.assertEqual(report["contrasts"], CONTRAST_ORDER)
        self.assertEqual(report["primary_excluded_samples"], ["17"])
        self.assertEqual(report["sensitivity_restored_samples"], ["17"])
        self.assertIn("python", report)
        self.assertIn("packages", report)
        self.assertEqual(len(report["volcano_svgs"]), 6)
        for spec in report["volcano_svgs"]:
            self.assertEqual(spec["width"], "44mm")
            self.assertEqual(spec["height"], "38mm")
            self.assertTrue(spec["editable_arial_text"])

    def test_reference_verifier_and_comparison_detect_changes(self):
        self.assertEqual(verify_reference_results(self.repo), [])
        self.assertEqual(compare_reference_results(self.generated, self.reference), [])
        candidate = self.reference / "tables" / "candidate_summary.csv"
        original = candidate.read_text(encoding="utf-8")
        candidate.write_text(original + "tampered\n", encoding="utf-8")
        try:
            errors = compare_reference_results(self.generated, self.reference)
            self.assertTrue(any("candidate_summary.csv" in error for error in errors))
        finally:
            candidate.write_text(original, encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
