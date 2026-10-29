"""Validate the six-culture case-study contract; preserve all raw-library rows."""
import csv, json, shutil
from pathlib import Path
def prepare(config,work,package):
 cfg=json.loads(Path(config).read_text()); p=Path(cfg['counts'])
 if not p.is_absolute():p=Path(config).resolve().parent/p
 if cfg['contrast'] not in ['mTKO_B16_T17','mTKO_MC38_T15','mTKO_4T1_T25']:raise ValueError('This release runner implements the three frozen Lawson primary contrasts')
 if cfg['labels']!=[0,0,0,1,1,1] or len(cfg['samples'])!=6:raise ValueError('Six unpaired cultures, comparator first, are required')
 if not cfg.get('replicate_provenance'):raise ValueError('Replicate provenance required')
 seen=set(); n=0
 with p.open() as f:
  r=csv.reader(f,delimiter='\t'); head=next(r)
  if head!=['sgRNA','Gene']+cfg['samples']:raise ValueError('Expected sgRNA, Gene, then exactly the configured six samples')
  for row in r:
   if len(row)!=8 or not row[0] or not row[1] or row[0] in seen:raise ValueError('Missing or duplicate guide, gene, or malformed row')
   seen.add(row[0]); n+=1
   if any(not v.isdecimal() for v in row[2:]):raise ValueError('Counts must be nonnegative integers; CPM is not accepted')
 if n<3:raise ValueError('Insufficient guides')
 for d in ['10_code','01_analysis_freeze','04_screen_reconstruction/phase1','04_screen_reconstruction/phase1rr','06_qc/phase1rr','07_tables/phase1rr','11_logs']:(work/d).mkdir(parents=True,exist_ok=True)
 for script in (package/'src/10_code').iterdir():
  if not script.name.startswith('._'):shutil.copyfile(script,work/'10_code'/script.name)
 target='04_screen_reconstruction/phase1/'+cfg['contrast']+'.raw.tsv';shutil.copyfile(p,work/target)
 spec=[dict(contrast=cfg['contrast'],raw=target,samples=cfg['samples'],role='PRIMARY')]
 (work/'01_analysis_freeze/CONTRAST_SPECIFICATIONS.json').write_text(json.dumps(spec,indent=2))
 (work/'INPUT_VALIDATION.json').write_text(json.dumps(dict(status='PASS',guides=n,contrast=cfg['contrast'],raw_input=str(p.name)),indent=2))
 return cfg
