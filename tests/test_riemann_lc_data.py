"""Checks the scientific reductions and rejects incomplete paired evidence."""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('riemann_lc', ROOT/'scripts/build_riemann_lc_data.py')
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)

class RiemannLCTables(unittest.TestCase):
    def test_shock_width_and_mass_drift_are_distinct(self):
        with (ROOT/'generated/riemann_lc_summary.csv').open() as stream:
            rows = {(r['scheme'],r['s']):r for r in csv.DictReader(stream)}
        self.assertAlmostEqual(float(rows['NN1','4']['width_standard']),1.4276033,places=5)
        self.assertAlmostEqual(float(rows['BAR2','2']['width_standard']),2.20509,places=5)
        self.assertGreater(float(rows['CV2','2']['mass_standard']),5e-6)
        self.assertLess(float(rows['CV2','2']['mass_lc']),5e-13)

    def test_incomplete_pair_set_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            shutil.copytree(ROOT/'data/riemann_lc',root/'data/riemann_lc')
            (root/'provenance').mkdir()
            manifest=json.loads((ROOT/'provenance/riemann_lc_inputs.json').read_text())
            source=root/'data/riemann_lc/paired_lc_vs_off.csv'
            lines=source.read_bytes().splitlines(keepends=True)
            source.write_bytes(b''.join(lines[:-1]))
            manifest['sha256'][source.name]=hashlib.sha256(source.read_bytes()).hexdigest()
            (root/'provenance/riemann_lc_inputs.json').write_text(json.dumps(manifest))
            old_root,old_source=BUILDER.ROOT,BUILDER.SOURCE
            try:
                BUILDER.ROOT,BUILDER.SOURCE=root,root/'data/riemann_lc'
                with self.assertRaisesRegex(ValueError,'36 unique matched'):
                    BUILDER.main()
            finally:
                BUILDER.ROOT,BUILDER.SOURCE=old_root,old_source

if __name__=='__main__':
    unittest.main()
