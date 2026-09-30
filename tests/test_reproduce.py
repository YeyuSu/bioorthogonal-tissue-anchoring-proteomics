import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from reproduce import run_reproduction

ROOT = Path(__file__).resolve().parents[1]


class ReproduceTests(unittest.TestCase):
    def test_reproduce_rejects_wrong_input_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "repo"
            (root / "data").mkdir(parents=True)
            (root / "configs").mkdir()
            (root / "data" / "27_mice_protein_groups.tsv").write_text("corrupt\n", encoding="utf-8")
            shutil.copy(ROOT / "configs" / "27_mice.json", root / "configs" / "27_mice.json")
            output = Path(tmp) / "output"

            with self.assertRaisesRegex(ValueError, "SHA256"):
                run_reproduction(root, output)

            self.assertFalse(output.exists())

    def test_reproduce_rejects_nonempty_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "output"
            output.mkdir()
            (output / "keep.txt").write_text("do not replace", encoding="utf-8")

            with self.assertRaisesRegex(FileExistsError, "not empty"):
                run_reproduction(ROOT, output)

            self.assertEqual((output / "keep.txt").read_text(encoding="utf-8"), "do not replace")

    def test_reproduce_is_independent_of_cwd(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output = tmp_path / "generated"
            original = Path.cwd()
            try:
                os.chdir(tmp_path)
                result = run_reproduction(ROOT, output)
            finally:
                os.chdir(original)

            self.assertEqual(result, output.resolve())
            self.assertTrue((output / "tables" / "sample_metadata.csv").is_file())
            self.assertTrue((output / "tables" / "differential_proteins.csv").is_file())
            self.assertTrue((output / "validation_report.json").is_file())


if __name__ == "__main__":
    unittest.main()
