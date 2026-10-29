#!/usr/bin/env python3
"""Portable release entry point; no writes to the frozen analysis archive."""
from pathlib import Path
import argparse,importlib.util,json,subprocess,csv,shutil,sys
P=Path(__file__).resolve().parent
def load(path):
 s=importlib.util.spec_from_file_location(Path(path).stem,P/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def command(argv,log):
 with log.open('w') as f:subprocess.run(argv,stdout=f,stderr=subprocess.STDOUT,check=True)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['demo','full','audit']);ap.add_argument('--config');ap.add_argument('--output',required=True);ap.add_argument('--preflight-only',action='store_true');a=ap.parse_args();w=Path(a.output).resolve()
 if w.exists():raise SystemExit('Output must be a new directory; refusing to overwrite an existing run')
 w.mkdir(parents=True)
 if a.mode=='demo':
  m=load('04_gene_inference/empirical.py');p=[m.empirical_p(2,[-2,0,1,2]),m.empirical_p(0,[-2,0,1,2])]
  out={'mode':'SOFTWARE_DEMONSTRATION_ONLY','seed':'not stochastic','P':p,'BH':m.bh(p),'scientific_inference':False,'note':'Toy finite-reference arithmetic; not a simulated benchmark or biological analysis'}
  (w/'demo.json').write_text(json.dumps(out,indent=2));print(json.dumps(out));return
 if a.mode=='audit':
  print(json.dumps(load('07_external_benchmarks/audit.py').audit(P,w/'external_audit.json')));return
 if not a.config:raise SystemExit('--config is required for full reproduction')
 cfg=load('01_prepare_inputs/prepare.py').prepare(a.config,w,P);rs=shutil.which(cfg.get('rscript','Rscript'))
 if not rs:raise SystemExit('Rscript unavailable; install R 4.4.3 and the recorded dependencies')
 command([rs,'-e',"library(edgeR);library(data.table);library(jsonlite);stopifnot(as.character(packageVersion('edgeR'))=='4.4.0');sessionInfo()"],w/'11_logs/environment.log')
 if a.preflight_only:print('PREFLIGHT_PASS');return
 cid=cfg['contrast']; command([rs,str(w/'10_code/phase1rr_fidelity.R'),str(w),cid],w/'11_logs/fidelity.log')
 with (w/('06_qc/phase1rr/'+cid+'.NULL_FIDELITY.tsv')).open() as f:checks=list(csv.DictReader(f,delimiter='\t'))
 good=all(r['PASS']=='TRUE' for r in checks if r['scope']=='MANDATORY')
 if not good:
  (w/'STATUS.json').write_text(json.dumps({'SCREEN_FIDELITY_STATUS':'FAIL','CALIBRATION_STATUS':'NOT_RUN','formal_P_FDR':'WITHHELD'}));print('FIDELITY_FAIL: no formal empirical inference');return
 command([rs,str(w/'10_code/phase1rr_simulate.R'),str(w),cid,str(cfg.get('cores',3))],w/'11_logs/simulation.log')
 with (w/('06_qc/phase1rr/'+cid+'.CALIBRATION_SUMMARY.tsv')).open() as f:cal=next(csv.DictReader(f,delimiter='\t'))
 with (w/('06_qc/phase1rr/'+cid+'.SPIKE_SUMMARY.tsv')).open() as f:sp=next(csv.DictReader(f,delimiter='\t'))
 if cal['CALIBRATION_PASS']!='TRUE' or sp['SPIKE_STATUS']!='PASS':
  (w/'STATUS.json').write_text(json.dumps({'SCREEN_FIDELITY_STATUS':'PASS','CALIBRATION_STATUS':'NOT_QUALIFIED','formal_P_FDR':'PROVISIONAL_NOT_FOR_INTERPRETATION'}));return
 command([rs,str(P/'05_robustness/robustness.R'),str(w),cid],w/'11_logs/robustness.log')
 if cid=='mTKO_B16_T17':
  t=w/'07_tables/phase1rr';load('06_functional_reclassification/classify.py').classify(t/(cid+'.EMPIRICAL_INFERENCE.tsv'),t/(cid+'.RELEASE_ROBUSTNESS.tsv'),P/'data/B16_gene_inference.tsv',P/'data/CANDIDATE_EVIDENCE_MATRIX.tsv',w/'B16_RECONSTRUCTED_CLASSIFICATION.tsv')
 (w/'STATUS.json').write_text(json.dumps({'SCREEN_FIDELITY_STATUS':'PASS','CALIBRATION_STATUS':'PASS','downstream':'B16 classifications use explicitly frozen ancillary and biological evidence inputs'}));print('FULL_STATISTICAL_RUN_COMPLETE')
if __name__=='__main__':main()
