import unittest,importlib.util,tempfile,json,subprocess,sys,csv
from pathlib import Path
P=Path(__file__).resolve().parents[1]
def mod(p):
 s=importlib.util.spec_from_file_location('unit_module',P/p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
class ReleaseTests(unittest.TestCase):
 def test_ties_and_plus_one(self):
  m=mod('04_gene_inference/empirical.py');self.assertEqual(m.empirical_p(2,[-2,0,1,2]),.6);self.assertEqual(m.empirical_p(9,[0,1]),1/3)
 def test_bh_and_bad_input(self):
  m=mod('04_gene_inference/empirical.py');self.assertEqual(m.bh([.01,.04,.03]),[.03,.04,.04]);self.assertRaises(ValueError,m.bh,[float('nan')]);self.assertRaises(ValueError,m.empirical_p,1,[])
 def test_frozen_data_counts(self):
  with (P/'data/CANDIDATE_EVIDENCE_MATRIX.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
  self.assertEqual(len(r),1492);self.assertEqual(sum(x['tier']=='TIER_A' for x in r),0);self.assertEqual(sum(x['residual']=='True' for x in r),1354)
 def test_demo_is_explicitly_noninferential(self):
  # Temp files stay under the requested release workspace.
  with tempfile.TemporaryDirectory(dir=P/'tests') as d:
   out=Path(d)/'run';subprocess.run([sys.executable,str(P/'run.py'),'demo','--output',str(out)],check=True,stdout=subprocess.DEVNULL);self.assertFalse(json.loads((out/'demo.json').read_text())['scientific_inference'])
 def test_reject_fractional_counts(self):
  with tempfile.TemporaryDirectory(dir=P/'tests') as d:
   d=Path(d);(d/'bad.tsv').write_text('sgRNA\tGene\tA\tB\tC\tD\tE\tF\ng1\tG\t0.2\t1\t1\t1\t1\t1\n');cfg=dict(contrast='mTKO_B16_T17',counts='bad.tsv',samples=list('ABCDEF'),labels=[0,0,0,1,1,1],replicate_provenance='test');(d/'cfg.json').write_text(json.dumps(cfg));self.assertRaises(ValueError,mod('01_prepare_inputs/prepare.py').prepare,d/'cfg.json',d/'out',P)
 def test_reconstruct_frozen_classification(self):
  import pandas as pd
  g=pd.read_csv(P/'data/B16_gene_inference.tsv',sep='\t')
  with tempfile.TemporaryDirectory(dir=P/'tests') as d:
   d=Path(d);g[['gene','sample_robust','guide_robust']].to_csv(d/'robust.tsv',sep='\t',index=False)
   x=mod('06_functional_reclassification/classify.py').classify(P/'data/B16_gene_inference.tsv',d/'robust.tsv',P/'data/B16_gene_inference.tsv',P/'data/CANDIDATE_EVIDENCE_MATRIX.tsv',d/'reconstructed.tsv')
   self.assertEqual(int(x.HIGH_CONFIDENCE.sum()),1492);self.assertEqual(int(x.residual.sum()),1354);self.assertEqual(int(x.tier.eq('TIER_A').sum()),0);self.assertEqual(int(x.tier.eq('TIER_B').sum()),179);self.assertEqual(int(x.tier.eq('TIER_C').sum()),1313)
if __name__=='__main__':unittest.main()
