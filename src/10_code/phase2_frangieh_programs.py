from pathlib import Path
import json,hashlib,pandas as pd
R=Path(__file__).resolve().parents[1];M=R/'03_metadata/phase2';f=json.loads((R/'01_analysis_freeze/PHASE2_METHODS_DECISIONS.json').read_text());rows=[]
for n,gg in f['Frangieh']['programs'].items():rows.append(dict(name=n,source='Prespecified Frangieh mechanistic axes',version='2026-10-20 before phenotype inspection',genes=';'.join(gg),known_or_exploratory='KNOWN',citation='https://doi.org/10.1038/s41588-021-00779-1',native_condition='ALL'))
s=pd.read_csv(M/'FRANGIEH_S5_ORIGINAL_ROWS.tsv',sep='\t');orig=[json.loads(a) for a in s['values']];head=orig[37]
for j,n in enumerate(head):
 if not n:continue
 gg=sorted({a[j] for a in orig[38:] if isinstance(a[j],str) and a[j].strip()});rows.append(dict(name='AUTHOR_JACKSTRAW_'+n,source='Frangieh TableS5 original program members',version='2021 original publication',genes=';'.join(gg),known_or_exploratory='KNOWN',citation='https://doi.org/10.1038/s41588-021-00779-1',native_condition='Co-culture' if n.startswith('Post_CocX') else 'IFNγ' if n.startswith('Treated') else 'Control'))
s=pd.read_csv(M/'FRANGIEH_S6.tsv',sep='\t');s['Feature Cluster']=pd.to_numeric(s['Feature Cluster'],errors='coerce')
for k,z in s.dropna(subset=['Feature Cluster']).groupby('Feature Cluster'):
 gg=sorted(set(z.Feature.dropna().astype(str)));rows.append(dict(name=f'AUTHOR_COFUNCTIONAL_FEATURE_CLUSTER_{int(k)}',source='Frangieh TableS6 original co-regulated features; exploratory phenotype axes, no novel module naming',version='2021 original publication',genes=';'.join(gg),known_or_exploratory='KNOWN',citation='https://doi.org/10.1038/s41588-021-00779-1',native_condition='ALL'))
p=M/'FRANGIEH_PREDEFINED_PROGRAMS.tsv';pd.DataFrame(rows).to_csv(p,sep='\t',index=False);p.with_suffix('.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n');print(len(rows),'programs frozen')
