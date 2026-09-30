import tempfile
import unittest
from pathlib import Path

from scripts.validate_release import (
    audit_repository,
    verify_license_coverage,
    verify_reference_hashes,
)


class ReleaseScopeTests(unittest.TestCase):
    def test_release_scope_rejects_out_of_scope_material(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data").mkdir()
            (root / "configs").mkdir()
            (root / "scripts").mkdir()
            (root / "results").mkdir()
            (root / "data" / "B16F10_input.tsv").write_text("protein\n", encoding="utf-8")
            (root / ".DS_Store").write_bytes(b"metadata")
            (root / "scripts" / "cache.pyc").write_bytes(b"bytecode")
            (root / "configs" / "analysis.json").write_text(
                '{"contrast": "NHS_vs_NCS"}', encoding="utf-8"
            )

            errors = audit_repository(root)

            joined = "\n".join(errors)
            self.assertIn("B16F10", joined)
            self.assertIn(".DS_Store", joined)
            self.assertIn(".pyc", joined)
            self.assertIn("NHS_vs_NCS", joined)

    def test_release_scope_rejects_machine_specific_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "scripts").mkdir()
            machine_path = "/" + "Users" + "/" + "syysmac"
            (root / "scripts" / "example.py").write_text(
                f'INPUT = "{machine_path}/Documents/input.tsv"\n', encoding="utf-8"
            )

            errors = audit_repository(root)

            self.assertTrue(any(machine_path in error for error in errors))

    def test_readme_may_state_that_b16f10_is_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "README.md").write_text(
                "B16F10 data are not included in this repository.\n", encoding="utf-8"
            )

            self.assertEqual(audit_repository(root), [])


class ReleaseDocumentationTests(unittest.TestCase):
    def test_reference_hashes_match_frozen_release(self):
        root = Path(__file__).resolve().parents[1]
        self.assertEqual(verify_reference_hashes(root), [])

    def test_license_coverage_is_explicit(self):
        root = Path(__file__).resolve().parents[1]
        self.assertEqual(verify_license_coverage(root), [])
        self.assertIn("MIT License", (root / "LICENSE-CODE").read_text(encoding="utf-8"))
        data_license = (root / "LICENSE-DATA").read_text(encoding="utf-8")
        self.assertIn("CC BY 4.0", data_license)
        self.assertIn("https://creativecommons.org/licenses/by/4.0/", data_license)

    def test_readme_documents_reproduction_and_statistical_limits(self):
        root = Path(__file__).resolve().parents[1]
        readme = (root / "README.md").read_text(encoding="utf-8")
        for required in [
            "python reproduce.py",
            "27 mice",
            "mouse 17",
            "Welch",
            "Benjamini–Hochberg",
            "streptavidin",
            "results/reference",
            "B16F10 data are not included",
        ]:
            self.assertIn(required, readme)


if __name__ == "__main__":
    unittest.main()
