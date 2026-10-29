from pathlib import Path
from collections import defaultdict
import pandas as pd,numpy as np,json,openpyxl,hashlib
from scipy.stats import rankdata
R=Path(__file__).resolve().parents[1];T=R/'07_tables';M=R/'03_metadata';E=T/'phase2';Q=R/'06_qc/phase2';D=R/'02_raw_data/phase2'
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def out(d,p):pd.DataFrame(d).to_csv(p,sep='\t',index=False)
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
# HGNC-approved aliases are an identifier resolution layer, never name-similarity mapping.
hg=rd(D/'HGNC_complete_set.tsv').fillna('');alias=defaultdict(set);approved=set(hg.symbol)
for z in hg.to_dict('records'):
 alias[z['symbol']].add(z['symbol'])
 for k in ['alias_symbol','prev_symbol']:
  for a in str(z[k]).split('|'):
   if a:alias[a].add(z['symbol'])
def canon(s):
 if s in approved:return s
 v=alias.get(s,set());return next(iter(v)) if len(v)==1 else s
out([dict(input_symbol=k,approved_symbol=canon(k),n_possible=len(v),source='HGNC complete set retrieved2026-10-20') for k,v in sorted(alias.items())],M/'phase2/HUMAN_SYMBOL_RESOLUTION.tsv')
# All primary mouse->human mappings remain the frozen MGI one-to-one mappings.
o=rd(M/'MOUSE_HUMAN_ORTHOLOGY.tsv');oo=o[o.one_to_one].copy();oo['human_gene']=oo.human_gene.map(canon);oo=oo[['mouse_gene','human_gene']];assert not oo.mouse_gene.duplicated().any() and not oo.human_gene.duplicated().any()
g=rd(T/'phase1rr/PHASE1RR_GENE_INFERENCE.tsv').query("contrast=='mTKO_B16_T17'");g=g.merge(oo,left_on='gene',right_on='mouse_gene',how='left');h=g[g.HIGH_CONFIDENCE].copy();h['arm']=np.where(h.median_log2FC>0,'ESCAPE_ARM','SENSITIZER_ARM')
# RefSeq->gene supplied by Joung authors. Drop ambiguous RefSeq assignments.
ann=rd(M/'phase2/JOUNG_AUTHOR_GENE_REFSEQ.tsv')[['Gene','id']].dropna();ann.Gene=ann.Gene.map(canon);counts=ann.groupby('id').Gene.nunique();ann=ann[ann.id.map(counts)==1].drop_duplicates('id').set_index('id').Gene
cache=D/'JOUNG_CPM_WITH_AUTHOR_GENE.tsv.gz'
if cache.exists():j=rd(cache)
else:
 w=openpyxl.load_workbook(R/'02_raw_data/phase0r/Joung_MOESM5.xlsx',read_only=True,data_only=True);v=w['sgRNA CPM'].values;cols=next(v);j=pd.DataFrame(v,columns=cols);w.close();j['human_gene']=j['RefSeq ID'].map(ann);j.to_csv(cache,sep='\t',index=False,compression='gzip')
valid=j.human_gene.notna();j=j[valid].copy();rows=[]
contrasts={'ANTIGEN_SPECIFIC_INCREMENT':('ESOTcells','Tcells'),'TOTAL_T_CELL_CONTEXT':('ESOTcells','noTcells'),'UNMODIFIED_T_CELL_CONTEXT':('Tcells','noTcells')}
for name,(tr,ct) in contrasts.items():
 for rep in [1,2]:
  a=j[f'Acute {tr} S{rep}'].to_numpy(float);b=j[f'Acute {ct} S{rep}'].to_numpy(float);j[f'{name}_S{rep}']=np.log2((a+1)/(b+1))
 cols=[name+'_S1',name+'_S2'];z=j.groupby('human_gene')[cols].median();z['effect']=z.median(axis=1);z['n_guides']=j.groupby('human_gene').size();z['percentile']=z.effect.rank(pct=True);z['donor_direction_consistent']=np.sign(z[cols[0]])==np.sign(z[cols[1]]);z['contrast']=name;z['modality']='CRISPRa';z['gene']=z.index;z['human_gene']=z.index;out(z.reset_index(drop=True),E/f'JOUNG_{name}_ALL_GENES.tsv')
 x=h.merge(z.reset_index(drop=True),on='human_gene',suffixes=('_B16','_Joung'));x['same_direction_descriptive']=np.sign(x.median_log2FC)==np.sign(x.effect);x['opposite_direction_descriptive']=np.sign(x.median_log2FC)==-np.sign(x.effect);x['tail_same_coding']=np.where(x.median_log2FC>0,x.percentile,1-x.percentile);x['tail_inverse_coding']=np.where(x.median_log2FC>0,1-x.percentile,x.percentile);x['support_status']='NOT_DIRECTIONALLY_COMPARABLE_CRISPRA_VS_KO';x['independent_KO_directional_support']=False;rows.append(x)
out(pd.concat(rows,ignore_index=True).rename(columns={'contrast_Joung':'external_contrast','contrast_B16':'B16_contrast'}),T/'JOUNG_DIRECTIONAL_VALIDATION.tsv')
# Set-level matched-universe permutation. Same/inverse coding both shown as cross-modality diagnostics.
nullrows=[];rng=np.random.default_rng(20261020)
for name in contrasts:
 z=rd(E/f'JOUNG_{name}_ALL_GENES.tsv');matched=g.merge(z,on='human_gene',suffixes=('_B16','_Joung'));vals=matched.percentile.to_numpy(float);zgenes=set(matched.human_gene)
 for arm,sign in [('ESCAPE_ARM',1),('SENSITIZER_ARM',-1)]:
  a=matched[matched.HIGH_CONFIDENCE & (np.sign(matched.median_log2FC)==sign)];obs=float(np.mean(sign*(a.percentile-.5)));n=len(a);null=np.array([np.mean(sign*(rng.choice(vals,n,replace=False)-.5)) for _ in range(10000)]);psame=(1+sum(null>=obs))/10001;pinv=(1+sum(null<=obs))/10001
  nullrows.append(dict(contrast=name,arm=arm,n_HC_overlap=n,n_matched_background=len(matched),mean_directional_percentile=.5+obs,same_coding_P=psame,inverse_coding_P=pinv,interpretation='CRISPRa_vs_KO directional-coding sensitivity only; neither is automatic replication'))
x=pd.DataFrame(nullrows);x['same_coding_FDR']=bh(x.same_coding_P);x['inverse_coding_FDR']=bh(x.inverse_coding_P);out(x,E/'JOUNG_SET_LEVEL_PERMUTATION.tsv')
# Yeh targeted author-score support; keep discordance with highlighted validation rather than changing sign per gene.
y=rd(D/'YEH_AUTHOR_Z_SCORES.tsv');y.gene=y.gene.map(canon);y=y.rename(columns={'gene':'human_gene'});z=h.merge(y,on='human_gene');z['T_cell_direction']=np.sign(z['Tcell.z']);z['NK_direction']=np.sign(z['NK.z']);z['T_support']=(np.sign(z.median_log2FC)==z.T_cell_direction)&(z['Tcell.z'].abs()>=1);z['NK_support']=(np.sign(z.median_log2FC)==z.NK_direction)&(z['NK.z'].abs()>=1);z['category']=np.select([z.T_support&z.NK_support,z.T_support&~z.NK_support,z.NK_support&~z.T_support],['T_NK_SHARED','T_CELL_SHARED','NK_DOMINANT'],default='DISCORDANT');z['formal_global_validation']=False;z['source']='Yeh Fig7b 74-gene author z, positive resistance; targeted |z|>=1 descriptive, not new P';out(z,T/'YEH_DIRECTIONAL_SUPPORT.tsv')
# Akana gene activation context. Do not infer KO concordance by inverting the signs.
a=rd(M/'phase2/AKANA_Table_S2b_CRISPRa_TCR.tsv');a['human_gene']=a.Genes.map(canon);z=h.merge(a,on='human_gene');z['activation_effect']=z[['S1.log2FC','S2.log2FC']].median(axis=1);z['directional_status']='NOT_DIRECTIONALLY_COMPARABLE';z['independent_KO_replication']=False;z['source']='Akana TableS2b; CRISPRa not KO';out(z,T/'AKANA_CROSS_MODALITY_CONTEXT.tsv')
# Frangieh actual target list used for overlap; harmonized assignment list verified when available.
f=rd(M/'phase2/FRANGIEH_S1.tsv');f['target']=f['Guide Name'].str.rsplit('_',n=1).str[0].map(canon);targets=set(f.target.dropna());fo=h.merge(pd.DataFrame({'human_gene':sorted(targets)}),on='human_gene',how='inner');fo['overlap_definition']='frozen strict1to1 + exact official-alias-resolved author targeted gene';out(fo,T/'FRANGIEH_OVERLAP.tsv')
lf=[]
for ratio,sfx in [(1,'e'),(2,'i'),(4,'m')]:
 z=rd(M/f'phase2/FRANGIEH_S2{sfx}.tsv');z['human_gene']=z.id.map(canon);z['LFC']=pd.to_numeric(z['pos|lfc'],errors='coerce');lf.append(z[['human_gene','LFC']].rename(columns={'LFC':f'TIL{ratio}_LFC'}).set_index('human_gene'))
fv=pd.concat(lf,axis=1).reset_index();fv=fo.merge(fv,on='human_gene');v=fv[['TIL1_LFC','TIL2_LFC','TIL4_LFC']].to_numpy(float);sign=np.sign(fv.median_log2FC.to_numpy(float));fv['all_ratios_direction_consistent']=(np.sign(v)==sign[:,None]).all(axis=1);fv['pressure_monotonic_descriptive']=((np.diff(v,axis=1)*sign[:,None])>=0).all(axis=1)&fv.all_ratios_direction_consistent;fv['support_category']=np.where(fv.pressure_monotonic_descriptive,'PRESSURE_MONOTONIC_DESCRIPTIVE',np.where(fv.all_ratios_direction_consistent,'CONSISTENT_DIRECTION','DISCORDANT'));fv['inference']='n=1 perratio, descriptive only, no new P/FDR';out(fv,E/'FRANGIEH_BULK_VIABILITY_DESCRIPTIVE.tsv')
# Resolve published mask aliases for a coverage audit; original source rows unchanged.
pub=rd(M/'PUBLISHED_IMMUNE_CRISPR_HITS.tsv');pub['canonical_human_gene']=pub.gene.map(canon);out(pub,M/'phase2/PUBLISHED_MASK_HGNC_RESOLUTION.tsv');audit=[]
pre=rd(E/'PRE_EXTERNAL_EVIDENCE_MATRIX.tsv')
for z in pre.to_dict('records'):
 human=canon(z.get('human_gene',''));prior=pub[(pub.study!='Lawson2020')&pub.published_hit & ((pub.canonical_human_gene==human)|(pub.gene==z['gene']))] if human else pub[(pub.study!='Lawson2020')&pub.published_hit&(pub.gene==z['gene'])]
 if bool(len(prior))!=bool(z['other_published_hit']):audit.append(dict(mouse_gene=z['gene'],human_gene=human,old_other_hit=z['other_published_hit'],resolved_other_hit=bool(len(prior)),studies=';'.join(sorted(set(prior.study)))))
out(audit if audit else pd.DataFrame(columns=['mouse_gene','human_gene','old_other_hit','resolved_other_hit','studies']),Q/'HGNC_MASK_DIFFERENCES.tsv')
summary={'n_HC':1492,'n_one_to_one_HC':int(h.human_gene.notna().sum()),'Frangieh_overlap_author_target':len(fo),'Joung_best_contrast_overlap':int(pd.concat(rows).query("contrast_Joung=='ANTIGEN_SPECIFIC_INCREMENT'").shape[0]),'Joung_modality':'CRISPRa; no direct KO replication','Yeh_overlap':len(z) if False else len(rd(T/'YEH_DIRECTIONAL_SUPPORT.tsv')),'Akana_overlap':len(rd(T/'AKANA_CROSS_MODALITY_CONTEXT.tsv')),'HGNC_prior_hit_differences':len(audit)};(Q/'EXTERNAL_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(summary)
