import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.core import run_proteomics_workflow as workflow

ROOT = Path(__file__).resolve().parents[1]


class CoreWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "configs" / "27_mice.json").read_text(encoding="utf-8"))

    def test_27_samples_have_expected_tissue_and_treatment_mapping(self):
        samples = self.config["samples"]
        self.assertEqual(len(samples), 27)
        expected = {
            **{str(i): ("subcutaneous", "biotin") for i in range(1, 4)},
            **{str(i): ("subcutaneous", "NHS") for i in range(4, 7)},
            **{str(i): ("subcutaneous", "NCS") for i in range(7, 10)},
            **{str(i): ("muscle", "biotin") for i in range(10, 13)},
            **{str(i): ("muscle", "NHS") for i in range(13, 16)},
            **{str(i): ("muscle", "NCS") for i in range(16, 19)},
            **{str(i): ("periodontal", "biotin") for i in range(19, 22)},
            **{str(i): ("periodontal", "NHS") for i in range(22, 25)},
            **{str(i): ("periodontal", "NCS") for i in range(25, 28)},
        }
        observed = {item["sample_id"]: (item["tissue"], item["treatment"]) for item in samples}
        self.assertEqual(observed, expected)

    def test_primary_excludes_17_and_sensitivity_restores_it(self):
        self.assertEqual(self.config["analysis"]["exclude_samples"], ["17"])
        self.assertTrue(self.config["sensitivity"]["restore_excluded_samples"])

    def test_reference_and_six_treatment_vs_biotin_contrasts_are_fixed(self):
        self.assertEqual(
            self.config["normalization"],
            {
                "reference_match_column": "PG.ProteinAccessions",
                "reference_match_value": "CON__Streptavidin",
            },
        )
        contrasts = self.config["contrasts"]
        self.assertEqual(len(contrasts), 6)
        self.assertTrue(all(item["control"] == "biotin" for item in contrasts))
        self.assertFalse(any(item["treated"] == "NCS" and item["control"] == "NHS" for item in contrasts))
        self.assertFalse(any("NHS_vs_NCS" in item["name"] for item in contrasts))

    def test_bh_adjustment_matches_known_values(self):
        adjusted = workflow.bh_adjust(np.array([0.01, 0.04, 0.03, np.nan]))
        np.testing.assert_allclose(adjusted[:3], [0.03, 0.04, 0.04])
        self.assertTrue(np.isnan(adjusted[3]))

    def test_reference_normalization_preserves_missing_values(self):
        config = {
            "columns": {"protein_id": "ProteinID", "description": "Description"},
            "normalization": {
                "reference_match_column": "ProteinID",
                "reference_match_value": "CON__Streptavidin",
            },
        }
        samples = [
            {"sample_id": "1", "column": "S1"},
            {"sample_id": "2", "column": "S2"},
        ]
        table = pd.DataFrame(
            {
                "ProteinID": ["CON__Streptavidin", "P1"],
                "Description": ["Reference", "Protein"],
                "S1": [10.0, 20.0],
                "S2": [20.0, np.nan],
            }
        )
        reference = workflow.extract_reference_values(table, config, samples)
        normalized = workflow.normalize_table(table, reference, samples)
        self.assertEqual(normalized.loc[1, "S1"], 2.0)
        self.assertTrue(np.isnan(normalized.loc[1, "S2"]))


if __name__ == "__main__":
    unittest.main()
