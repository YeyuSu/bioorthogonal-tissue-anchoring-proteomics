import tempfile
import unittest
from pathlib import Path

from scripts.validate_release import audit_repository


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
            (root / "scripts" / "example.py").write_text(
                'INPUT = "/Users/syysmac/Documents/input.tsv"\n', encoding="utf-8"
            )

            errors = audit_repository(root)

            self.assertTrue(any("/Users/syysmac" in error for error in errors))

    def test_readme_may_state_that_b16f10_is_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "README.md").write_text(
                "B16F10 data are not included in this repository.\n", encoding="utf-8"
            )

            self.assertEqual(audit_repository(root), [])


if __name__ == "__main__":
    unittest.main()
