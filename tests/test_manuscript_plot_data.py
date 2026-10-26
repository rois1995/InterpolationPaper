"""Check pairing and sampling choices that could silently mislead the new figures."""
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('plots', ROOT/'scripts/build_manuscript_plot_data.py')
plots = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plots)

export_spec = importlib.util.spec_from_file_location('export', ROOT/'scripts/prepare_cf_submission.py')
export = importlib.util.module_from_spec(export_spec)
export_spec.loader.exec_module(export)


class FigureDataTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        for folder in ('data','provenance','RANS_PostProcessing'):
            shutil.copytree(ROOT/folder,self.root/folder)
        (self.root/'generated').mkdir()

    def test_error_endpoint_and_seed_envelopes(self):
        plots.build(self.root)
        f = pd.read_csv(self.root/'generated/riemann_error_plot.csv')
        pairs = pd.read_csv(self.root/'data/riemann_lc/paired_lc_vs_off.csv')
        self.assertEqual(len(f),42)
        for scheme in plots.SCHEMES:
            expected = pairs[(pairs.scheme == scheme) & (pairs.dIters == 2) & (pairs.SF == 2)].L1_post399_change_pct.mean()
            self.assertAlmostEqual(f.iloc[-1][scheme+'_mean'],expected,delta=.001)
            self.assertTrue((f[scheme+'_min'] <= f[scheme+'_mean']+1e-10).all())
            self.assertTrue((f[scheme+'_max'] >= f[scheme+'_mean']-1e-10).all())

    def test_flat_export_resolves_nested_macro_paths(self):
        source = r'{\DataDir/riemann_lc/density_profiles.csv} {\GeneratedDir/rans_drag_plot.csv} {data/riemann_lc/density_profiles.csv}'
        mapping = {'data/riemann_lc/density_profiles.csv':'density_profiles.csv',
                   'generated/rans_drag_plot.csv':'rans_drag_plot.csv'}
        self.assertEqual(export.flatten_source(source,mapping),
                         '{density_profiles.csv} {rans_drag_plot.csv} {density_profiles.csv}')

    def test_missing_paired_time_is_rejected(self):
        path = self.root/'data/riemann_lc/timeseries.csv'
        f = pd.read_csv(path)
        drop = f.index[(f['case'] == 'WD_0.1_NN_Every_2Steps_SF2_First_LC') & (f.step == 20)]
        self.assertEqual(len(drop),1)
        f.drop(drop).to_csv(path,index=False)
        with self.assertRaisesRegex(ValueError,'Unaligned paired error histories'):
            plots.build(self.root)


if __name__ == '__main__':
    unittest.main()
