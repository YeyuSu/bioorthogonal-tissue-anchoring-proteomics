import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.publication.build_publication_tables import (
    CONTRAST_ORDER,
    build_publication_tables,
)


class PublicationTableTests(unittest.TestCase):
    def make_core_output(self, root: Path, include_extra: bool = False) -> Path:
        tables = root / "tables"
        tables.mkdir(parents=True)
        rows = []
        for index, contrast in enumerate(CONTRAST_ORDER):
            tissue, treated = contrast.split("_", 1)[0], "NHS" if "_NHS_" in contrast else "NCS"
            for protein_id, description, log2fc, q in [
                ("P1", "Protein one", 3.0 - index * 0.1, 0.01 if index == 0 else 0.20),
                ("P2", "Protein two", 1.5 + index * 0.1, 0.20),
                ("P3", "Protein three", -2.0, 0.01 if index == 1 else 0.20),
                ("P4", "Protein four", -0.5, 0.90),
            ]:
                rows.append(
                    {
                        "analysis": "primary",
                        "contrast": contrast,
                        "tissue": tissue,
                        "treated": treated,
                        "control": "biotin",
                        "protein_id": protein_id,
                        "description": description,
                        "n_treated": 3,
                        "n_control": 3,
                        "log2FC": log2fc,
                        "p_Welch": min(q, 0.04),
                        "q_BH": q,
                        "significant_enriched": q < 0.05 and log2fc >= 1,
                        "significant_depleted": q < 0.05 and log2fc <= -1,
                    }
                )
        if include_extra:
            rows.append({**rows[0], "contrast": "subcutaneous_NHS_vs_NCS", "control": "NCS"})
        differential = pd.DataFrame(rows)
        differential.to_csv(tables / "differential_proteins.csv", index=False)
        sensitivity = differential.copy()
        sensitivity["analysis"] = "sensitivity"
        sensitivity.to_csv(tables / "sensitivity_differential_proteins.csv", index=False)
        pd.DataFrame({"sample_id": ["1"], "excluded_from_primary": [False]}).to_csv(
            tables / "sample_metadata.csv", index=False
        )
        pd.DataFrame({"sample_id": ["1"], "reference_abundance": [10.0]}).to_csv(
            tables / "sample_qc.csv", index=False
        )
        return root

    def test_builds_ranked_top10_and_six_column_heatmap_union(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            core = self.make_core_output(tmp_path / "core")
            destination = tmp_path / "publication"

            paths = build_publication_tables(core, destination)

            self.assertTrue(paths)
            top10 = pd.read_csv(destination / "top10_enriched_proteins.csv")
            self.assertTrue((top10["log2FC"] > 0).all())
            first = top10[top10["contrast"] == CONTRAST_ORDER[0]]
            self.assertEqual(first.iloc[0]["protein_id"], "P1")
            self.assertEqual(first.iloc[0]["rank"], 1)

            log2fc = pd.read_csv(destination / "three_tissue_effect_heatmap_log2FC.csv")
            q_bh = pd.read_csv(destination / "three_tissue_effect_heatmap_q_BH.csv")
            self.assertEqual(list(log2fc.columns[1:]), CONTRAST_ORDER)
            self.assertEqual(list(q_bh.columns[1:]), CONTRAST_ORDER)
            self.assertEqual(set(log2fc["protein_key"]), {"P1|Protein one", "P3|Protein three"})

    def test_rejects_non_treatment_vs_biotin_contrast(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            core = self.make_core_output(tmp_path / "core", include_extra=True)

            with self.assertRaisesRegex(ValueError, "six treatment-versus-biotin"):
                build_publication_tables(core, tmp_path / "publication")


if __name__ == "__main__":
    unittest.main()
