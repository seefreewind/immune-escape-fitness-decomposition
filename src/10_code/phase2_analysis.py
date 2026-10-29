from pathlib import Path
from collections import defaultdict
import os
for x in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[x]='1'
import numpy as np,pandas as pd,json,hashlib
from scipy.stats import hypergeom
from scipy.spatial.distance import cdist
import networkx as nx
import gseapy as gp
R=Path(__file__).resolve().parents[1];T=R/'07_tables';Q=R/'06_qc/phase2';M=R/'03_metadata';E=T/'phase2'
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def out(x,p):
 if '--stage-reconciled' in __import__('sys').argv:
  p=T/'phase2/reconciled'/p.relative_to(T);p.parent.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(x).to_csv(p,sep='\t',index=False)
def bh(p):
 p=np.asarray(p,float);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
for z in rd(R/'01_analysis_freeze/PHASE2_ANNOTATION_FREEZE.tsv').to_dict('records'):
 p=R/z['file']
 if z['file']=='03_metadata/PUBLISHED_IMMUNE_CRISPR_HITS.tsv':p=R/'01_analysis_freeze/phase2_initial/PUBLISHED_IMMUNE_CRISPR_HITS.tsv'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256']
for z in rd(R/'01_analysis_freeze/PHASE2_PUBLISHED_MASK_RECONCILIATION_FREEZE.tsv').to_dict('records'):assert hashlib.sha256((R/z['file']).read_bytes()).hexdigest()==z['sha256']
g=rd(T/'phase1rr/PHASE1RR_GENE_INFERENCE.tsv').query("contrast=='mTKO_B16_T17'").sort_values('gene').reset_index(drop=True);assert len(g)==18710 and g.HIGH_CONFIDENCE.sum()==1492;U=set(g.gene);h=g[g.HIGH_CONFIDENCE].copy();h['arm']=np.where(h.median_log2FC>0,'ESCAPE_ARM','SENSITIZER_ARM')
s=rd(M/'MECHANISM_GENESETS.tsv').fillna('');sg={z['name']:set(z['tested_genes'].split(';'))-set(['']) for z in s.to_dict('records')};known=defaultdict(set);broad=defaultdict(set)
for z in s.to_dict('records'):
 for gene in sg[z['name']]:
  if z['dominant_core']:known[gene].add(z['known_mechanism_class'])
  elif z['known_mechanism_class']!='OTHER':broad[gene].update(z['known_mechanism_class'].split(';'))
o=rd(M/'MOUSE_HUMAN_ORTHOLOGY.tsv').fillna('');oo=o[o.one_to_one].set_index('mouse_gene');hits=rd(M/'PUBLISHED_IMMUNE_CRISPR_HITS.tsv').fillna('');hc=rd(M/'phase2/HUMAN_SYMBOL_RESOLUTION.tsv');cm=dict(zip(hc.input_symbol,hc.approved_symbol));canon=lambda z:cm.get(z,z);oo['human_gene']=oo.human_gene.map(canon);hits['gene']=hits.gene.map(canon);hs=hits[hits.published_hit & (hits.study!='Lawson2020')];screened=defaultdict(set);hitby=defaultdict(set)
for z in hits.to_dict('records'):screened[z['gene']].add(z['study'])
for z in hs.to_dict('records'):hitby[z['gene']].add(z['study'])
for i,z in h.iterrows():
 gene=z.gene;human=oo.loc[gene,'human_gene'] if gene in oo.index else '';prior=hitby[gene]|hitby[human];studies=screened[human]|screened[gene]
 h.loc[i,'human_gene']=human;h.loc[i,'one_to_one']=bool(human);h.loc[i,'known_core']=bool(z.published_core_gene);h.loc[i,'other_published_hit']=bool(prior);h.loc[i,'published_hit_studies']=';'.join(sorted(prior));h.loc[i,'known_mechanism_classes']=';'.join(sorted(known[gene]));h.loc[i,'broad_known_programs']=';'.join(sorted(broad[gene]));h.loc[i,'dominant_known_member']=bool(known[gene]);h.loc[i,'previously_screened']=';'.join(sorted(studies));h.loc[i,'residual']=not z.published_core_gene and not prior and not known[gene]
 h.loc[i,'novelty_mask_status']='KNOWN_PUBLISHED_HIT' if prior or z.published_core_gene else 'KNOWN_PATHWAY_NEW_MEMBER' if known[gene] or broad[gene] else 'PREVIOUSLY_SCREENED_NOT_HIGHLIGHTED' if studies else 'NO_CLEAR_PRIOR_CRISPR_HIT'
 h.loc[i,'residualization_class']='KNOWN_CORE' if z.published_core_gene else 'KNOWN_PUBLISHED_HIT' if prior else 'KNOWN_PATHWAY' if known[gene] else 'RESIDUAL'
 for st in ['Frangieh2021','Joung2022','Yeh2024','Akana2026']:h.loc[i,st+'_overlap']=st in studies
out(h,T/'KNOWN_MECHANISM_CLASSIFICATION.tsv');out(h,E/'PRE_EXTERNAL_EVIDENCE_MATRIX.tsv')
# Arm-specific ORA with fixed eligible background.
el=s[(s.n_tested>=15)&(s.n_tested<=500)&(~s.dominant_core)].copy();ora=[]
for arm,a in h.groupby('arm'):
 A=set(a.gene)
 for z in el.to_dict('records'):
  G=sg[z['name']];over=sorted(A&G);ora.append(dict(pathway=z['name'],direction=arm,k=len(over),n_HC=len(A),n_pathway=len(G),background_n=len(U),P=float(hypergeom.sf(len(over)-1,len(U),len(G),len(A))),genes=';'.join(over),source=z['source'],known_mechanism_class=z['known_mechanism_class']))
ora=pd.DataFrame(ora);ora['FDR']=bh(ora.P);out(ora,T/'ORA_PATHWAY_ENRICHMENT.tsv');print('ORA',len(ora),'significant',sum(ora.FDR<=.05),flush=True)
if '--reuse-enrichment' not in __import__('sys').argv:
 # Full ranked analysis, unique complete universe, fixed score, no HC selection.
 rnk=g[['gene','PRIMARY_GENE_SCORE']].sort_values(['PRIMARY_GENE_SCORE','gene'],ascending=[False,True]);res=gp.prerank(rnk=rnk,gene_sets={n:sg[n] for n in el.name},outdir=None,permutation_num=10000,threads=8,min_size=15,max_size=500,weight=1,seed=20261020,no_plot=True,verbose=False)
 x=res.res2d.copy();x=x.rename(columns={'Term':'pathway','NOM p-val':'P','Lead_genes':'leading_edge_genes','NES':'NES_equivalent','ES':'ES','FDR q-val':'GSEA_FDR','FWER p-val':'FWER_P'});x['P']=pd.to_numeric(x.P);x['NES_equivalent']=pd.to_numeric(x.NES_equivalent);x['FDR']=bh(x.P);x['direction']=np.where(x.NES_equivalent>0,'positive-tail','negative-tail');x['leading_edge_n']=x.leading_edge_genes.fillna('').map(lambda z:len(z.split(';')) if z else 0);x['known_mechanism_class']=x.pathway.map(s.set_index('name').known_mechanism_class);x['background_n']=18710;out(x,T/'RANKED_PATHWAY_ENRICHMENT.tsv');print('GSEA',len(x),'BH05',sum(x.FDR<=.05),flush=True)
# Evidence features: primary/early/legacy/rank/baseline. Different axes remain explicit in saved matrix.
feat=['median_log2FC','early_effect','author_NormZ','RRA_signed_strength','independent_count_rank','B_g'];base=h.set_index('gene');raw=base[feat].astype(float);med=raw.median();scale=(raw-med).abs().median().replace(0,1);standard=((raw-med)/scale).clip(-4,4).fillna(0);out(standard.reset_index(),E/'MODULE_EFFECT_FEATURES.tsv')
los=[]
for p in sorted((T/'phase1rr').glob('mTKO_B16_T17.OMIT_*.GENE_SCORES.tsv')):
 if not p.name.startswith('._'):los.append(rd(p).set_index('gene')['median_log2FC'])
assert len(los)==6
nonGO=el[~el.source.str.contains('GO_BP')];ngo=nonGO.name.tolist();go=el[el.source.str.contains('GO_BP')];rng=np.random.default_rng(20261020)
allassign=[];allmodules=[];bootrows=[];moduleGO=[]
def similarity_graph(A,F):
 weights=1/np.maximum(1,A.sum(axis=0));inter=(A*weights)@A.T;denom=(A*weights).sum(axis=1);j=inter/np.maximum(denom[:,None]+denom[None,:]-inter,1e-12);sim=.8*j+.2*np.exp(-cdist(F,F,metric='cityblock')/F.shape[1]);sim[j<.1]=0;np.fill_diagonal(sim,0);G=nx.Graph();G.add_nodes_from(range(len(A)))
 for i in range(len(A)):
  inds=np.argsort(-sim[i],kind='stable')[:15]
  for k in inds:
   if sim[i,k]>0:G.add_edge(i,int(k),weight=float(sim[i,k]))
 return G
for scope in ['ALL_HC','RESIDUAL_ONLY']:
 for arm,aa in h.groupby('arm'):
  aa=aa[aa.residual.astype(bool)] if scope=='RESIDUAL_ONLY' else aa;genes=aa.gene.tolist()
  if not genes:continue
  A=np.array([[float(gene in sg[n]) for n in ngo] for gene in genes]);active=A.sum(axis=0)>=2;A=A[:,active];F=standard.loc[genes].to_numpy(float);G=similarity_graph(A,F);com=list(nx.community.louvain_communities(G,weight='weight',resolution=1,seed=20261020)) if G.number_of_edges() else [{i} for i in range(len(genes))];com=sorted(com,key=lambda z:(-len(z),min(genes[k] for k in z)));labs=np.zeros(len(genes),int)
  for i,c in enumerate(com):labs[list(c)]=i
  matches=np.zeros((100,len(com)));cohesions=np.zeros((100,len(com)))
  for b in range(100):
   ids=rng.choice(A.shape[1],max(1,int(.8*A.shape[1])),replace=False);FB=F.copy();le=los[int(rng.integers(6))].reindex(genes).to_numpy(float);FB[:,0]=np.clip((le-med.iloc[0])/scale.iloc[0],-4,4);GB=similarity_graph(A[:,ids],FB);cb=list(nx.community.louvain_communities(GB,weight='weight',resolution=1,seed=20261020+b)) if GB.number_of_edges() else [{i} for i in range(len(genes))]
   lb=np.zeros(len(genes),int)
   for k,c in enumerate(cb):lb[list(c)]=k
   for k,c in enumerate(com):
    counts=np.bincount(lb[list(c)],minlength=len(cb));unions=np.array([len(c)+len(d) for d in cb])-counts;matches[b,k]=np.max(counts/np.maximum(unions,1));cohesions[b,k]=sum(v*(v-1) for v in counts)/max(len(c)*(len(c)-1),1)
  for k,c in enumerate(com):
   mid=f'{scope}_{arm}_{k+1:03d}';gg={genes[j] for j in c};stable=len(c)>=5 and np.median(matches[:,k])>=.65 and np.mean(matches[:,k]>=.5)>=.8
   allmodules.append(dict(module=mid,scope=scope,arm=arm,n_genes=len(c),genes=';'.join(sorted(gg)),median_bootstrap_Jaccard=np.median(matches[:,k]),fraction_bootstrap_Jaccard_ge_050=np.mean(matches[:,k]>=.5),mean_pair_coassignment=cohesions[:,k].mean(),stability_PASS=stable,graph_edges=G.subgraph(c).number_of_edges(),interpretation='PENDING_COHERENCE'))
   for gene in sorted(gg):allassign.append(dict(gene=gene,module=mid,scope=scope,arm=arm,module_stability=stable,degree=G.degree[genes.index(gene)]))
   for b in range(100):bootrows.append(dict(module=mid,resample=b+1,best_Jaccard=matches[b,k],pair_coassignment=cohesions[b,k]))
   for z in (go.to_dict('records') if len(c)>=5 else []):
    overlap=gg&sg[z['name']];moduleGO.append(dict(module=mid,pathway=z['name'],n_module=len(gg),k=len(overlap),P=float(hypergeom.sf(len(overlap)-1,18710,len(sg[z['name']]),len(gg))),genes=';'.join(sorted(overlap))))
  print('Modules',scope,arm,len(genes),len(com),'stable',sum(z['stability_PASS'] for z in allmodules if z['scope']==scope and z['arm']==arm),flush=True)
mg=pd.DataFrame(moduleGO);mg['FDR']=bh(mg.P);out(mg,E/'MODULE_HELDOUT_GO_ENRICHMENT.tsv');mods=pd.DataFrame(allmodules);best=mg.sort_values(['FDR','P']).groupby('module',sort=False).first();mods['coherent_GO_FDR']=mods.module.map(best.FDR);mods['coherent_GO_term']=mods.module.map(best.pathway);mods['coherence_PASS']=mods.coherent_GO_FDR<=.05;mods['stable_coherent']=mods.stability_PASS & mods.coherence_PASS;mods['interpretation']=np.where(mods.stable_coherent,'STABLE_ANNOTATION_COHERENT','UNRESOLVED_CLUSTER');out(mods,T/'FUNCTIONAL_MODULES.tsv');out(allassign,E/'MODULE_GENE_ASSIGNMENTS.tsv');out(bootrows,E/'MODULE_STABILITY_RESAMPLES.tsv');print('ANALYSIS_COMPLETE',flush=True)
