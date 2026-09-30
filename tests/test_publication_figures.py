import tempfile
import unittest
import warnings
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from scripts.core.run_proteomics_workflow import run_workflow
from scripts.publication.build_publication_tables import CONTRAST_ORDER, build_publication_tables
from scripts.publication.render_publication_figures import render_publication_figures

ROOT = Path(__file__).resolve().parents[1]


class PublicationFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.core = cls.root / "core"
        cls.tables = cls.root / "publication_tables"
        cls.figures = cls.root / "publication_figures"
        run_workflow(
            ROOT / "data" / "27_mice_protein_groups.tsv",
            ROOT / "configs" / "27_mice.json",
            cls.core,
        )
        build_publication_tables(cls.core, cls.tables)
        log2_path = cls.tables / "three_tissue_effect_heatmap_log2FC.csv"
        q_path = cls.tables / "three_tissue_effect_heatmap_q_BH.csv"
        log2fc = pd.read_csv(log2_path)
        q_bh = pd.read_csv(q_path)
        if not log2fc.empty:
            log2fc.iloc[0, 1] = float("nan")
            q_bh.iloc[0, 1] = float("nan")
            log2fc.to_csv(log2_path, index=False)
            q_bh.to_csv(q_path, index=False)
        with warnings.catch_warnings():
            warnings.simplefilter("error", PendingDeprecationWarning)
            cls.paths = render_publication_figures(cls.core, cls.tables, cls.figures)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_six_volcanoes_follow_svg_contract(self):
        volcanoes = sorted((self.figures / "standalone_volcanoes").glob("*.svg"))
        observed = [path.stem.removeprefix("Volcano_") for path in volcanoes]
        self.assertEqual(observed, sorted(CONTRAST_ORDER))
        self.assertEqual(len(volcanoes), 6)
        for path in volcanoes:
            root = ET.parse(path).getroot()
            self.assertEqual(root.attrib["width"], "44mm")
            self.assertEqual(root.attrib["height"], "38mm")
            svg = path.read_text(encoding="utf-8")
            self.assertIn("Arial", svg)
            self.assertIn("<text", svg)
            self.assertIn("font-size: 5px", svg)
            self.assertIn("font-size: 6px", svg)
            self.assertNotIn("NHS_vs_NCS", path.name)
            self.assertNotIn('id="Arial', svg)

    def test_heatmap_has_six_columns_exact_stars_and_missing_dash(self):
        heatmap = self.figures / "three_tissue_treatment_effect_heatmap.svg"
        self.assertTrue(heatmap.is_file())
        texts = [text.strip() for text in ET.parse(heatmap).getroot().itertext() if text.strip()]
        for contrast in CONTRAST_ORDER:
            self.assertIn(contrast, texts)

        log2fc = pd.read_csv(self.tables / "three_tissue_effect_heatmap_log2FC.csv")
        q_bh = pd.read_csv(self.tables / "three_tissue_effect_heatmap_q_BH.csv")
        expected_stars = 0
        expected_dashes = 0
        for contrast in CONTRAST_ORDER:
            expected_stars += int(((q_bh[contrast] < 0.05) & (log2fc[contrast].abs() >= 1)).sum())
            expected_dashes += int(log2fc[contrast].isna().sum())
        self.assertEqual(texts.count("*"), expected_stars)
        self.assertEqual(texts.count("–"), expected_dashes)


if __name__ == "__main__":
    unittest.main()
