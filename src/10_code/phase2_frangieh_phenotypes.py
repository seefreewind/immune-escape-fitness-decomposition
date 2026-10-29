from pathlib import Path
from collections import defaultdict
import os
for s in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS']:os.environ[s]='1'
import pandas as pd,numpy as np,json,hashlib,anndata as ad
from scipy import sparse
from sklearn.linear_model import ElasticNet
R=Path(__file__).resolve().parents[1];D=R/'02_raw_data/phase2';M=R/'03_metadata/phase2';T=R/'07_tables';E=T/'phase2';Q=R/'06_qc/phase2'
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def out(x,p):pd.DataFrame(x).to_csv(p,sep='\t',index=False)
resolution=rd(M/'HUMAN_SYMBOL_RESOLUTION.tsv');cm=dict(zip(resolution.input_symbol,resolution.approved_symbol));canon=lambda x:cm.get(x,x)
programs=rd(M/'FRANGIEH_PREDEFINED_PROGRAMS.tsv');ph=Path(M/'FRANGIEH_PREDEFINED_PROGRAMS.sha256');assert hashlib.sha256((M/'FRANGIEH_PREDEFINED_PROGRAMS.tsv').read_bytes()).hexdigest()==ph.read_text().split()[0]
fo=rd(T/'FRANGIEH_OVERLAP.tsv');targets=set(fo.human_gene)|{'CD58'}
a=ad.read_h5ad(D/'Frangieh_2021_PERTPY_PROCESSED.h5ad',backed='r');obs=a.obs.copy();obs['sgRNA']=obs.sgRNA.astype(str);obs['condition']=obs.condition.astype(str)
# The mirror's perturbation_name strips gene digits (IFNGR1 -> IFNGR); rebuild from exact sgRNA IDs.
obs['target']=obs.sgRNA.str.rsplit('_',n=1).str[0].map(canon);obs['NTC']=obs.sgRNA.str.startswith('NO_SITE_');obs['INTERGENIC']=obs.sgRNA.str.startswith('ONE_NON-GENE_SITE_');obs['target']=np.where(obs.NTC,'NTC',obs.target)
keep=obs.target.isin(targets)|obs.NTC|obs.INTERGENIC;obs=obs[keep].copy();ids=np.flatnonzero(keep.to_numpy());assert pd.to_numeric(obs.MOI).eq(1).all()
# Source already log1p CPM (10^6 per cell), verified inverse-transform sums; no second normalization.
varcanon=np.array([canon(str(x)) for x in a.var_names]);pg={z['name']:{canon(x) for x in z['genes'].split(';') if x} for z in programs.to_dict('records')};needed=set().union(*pg.values());colidx=np.flatnonzero(np.isin(varcanon,list(needed)));names=varcanon[colidx];unique,counts=np.unique(names,return_counts=True);ambig=set(unique[counts>1]);colidx=np.array([i for i in colidx if varcanon[i] not in ambig]);names=varcanon[colidx]
# Read sparse RNA in cell chunks, retain only exact perturbations and reference cells and predefined genes.
chunks=[]
for start in range(0,len(a.obs),4096):
 stop=min(len(a.obs),start+4096);sel=np.flatnonzero(keep.iloc[start:stop].to_numpy())
 if len(sel):chunks.append(a.X[start:stop,:][sel,:][:,colidx].tocsr())
X=sparse.vstack(chunks,format='csr');prot=np.asarray(a.obsm['protein_expression'][ids],float);pnames=[str(z).replace(' protein','') for z in a.uns['protein_names']];assert prot.shape[1]==20
# Exact row alignment and author isotype correction verified against independently deposited raw ADT.
b=ad.read_h5ad(D/'FrangiehIzar2021_protein.h5ad');b=b[obs.index,:].copy();B=b.X.toarray();bindex={n:i for i,n in enumerate(b.var_names)};expected=np.column_stack([np.maximum(0,np.log((B[:,bindex[p]]+1)/(B[:,bindex[str(b.var.loc[p,'Isotype_control'])]]+1))) for p in pnames]);maxdiff=float(np.max(np.abs(expected-prot)));assert maxdiff<2e-5,maxdiff
del b,B,expected
targets_obs=set(obs.loc[~obs.NTC & ~obs.INTERGENIC,'target']);out(fo.assign(in_actual_single_guide_assignments=fo.human_gene.isin(targets_obs)),T/'FRANGIEH_OVERLAP.tsv')
# Save feature membership/coverage before fitting. This is an anchor analysis, no sample-level statistical replication.
pn=programs.name.tolist();W=np.zeros((len(names),len(pn)))
coverage=[]
for j,n in enumerate(pn):
 z=np.isin(names,list(pg[n]));W[z,j]=1/max(z.sum(),1);coverage.append(dict(program=n,n_defined=len(pg[n]),n_measured=int(z.sum()),measured_genes=';'.join(names[z]),source=programs.iloc[j].source))
out(coverage,E/'FRANGIEH_PROGRAM_COVERAGE.tsv')
guide_rows=[];gene_rows=[];sizes=[];mapping_issues=[]
for condition in ['Control','IFNγ','Co-culture']:
 idx=np.flatnonzero(obs.condition.to_numpy()==condition);ob=obs.iloc[idx].copy();XX=X[idx];P=prot[idx];ntc=ob.NTC.to_numpy(bool);inter=ob.INTERGENIC.to_numpy(bool);train=~inter
 assert ntc.sum()>=100
 mu=np.asarray(XX[ntc].mean(axis=0)).ravel();sd=np.sqrt(np.maximum(0,np.asarray(XX[ntc].power(2).mean(axis=0)).ravel()-mu**2));sd[sd<1e-6]=1
 rscore=np.asarray(XX@(W/sd[:,None]))-(mu/sd)@W
 rmu=rscore[ntc].mean(axis=0);rsd=rscore[ntc].std(axis=0);rsd[rsd<1e-6]=1;rscore=(rscore-rmu)/rsd
 pmu=P[ntc].mean(axis=0);psd=P[ntc].std(axis=0);psd[psd<1e-6]=1;P=(P-pmu)/psd
 Y=np.column_stack([rscore,P]);labels=['RNA:'+z for z in pn]+['PROTEIN:'+z for z in pnames]
 guides=sorted(ob.loc[train&~ntc,'sgRNA'].unique());guide_index={s:i for i,s in enumerate(guides)};design=np.zeros((train.sum(),len(guides)+2));obtr=ob.iloc[np.flatnonzero(train)]
 for i,sg in enumerate(obtr.sgRNA):
  if sg in guide_index:design[i,guide_index[sg]]=1
 quality=np.log1p(pd.to_numeric(obtr.UMI_count).to_numpy(float));design[:,-2]=(quality-quality.mean())/max(quality.std(),1e-8);cyc=rscore[train,pn.index('CELL_STATE_CYCLE')];design[:,-1]=(cyc-cyc.mean())/max(cyc.std(),1e-8)
 y=Y[train];model=ElasticNet(alpha=.0005,l1_ratio=.5,max_iter=10000,precompute=True,tol=1e-4);model.fit(design,y);coefs=model.coef_[:,:len(guides)].copy()
 # A cycle score cannot be both response and nuisance for its own effect.
 jcycle=pn.index('CELL_STATE_CYCLE');m2=ElasticNet(alpha=.0005,l1_ratio=.5,max_iter=10000,precompute=True,tol=1e-4);m2.fit(design[:,:-1],y[:,jcycle]);coefs[jcycle]=m2.coef_[:len(guides)]
 for j,label in enumerate(labels):
  for sg,k in guide_index.items():
   select=ob.sgRNA.to_numpy()==sg;n=int(select.sum());target=ob.loc[select,'target'].iloc[0];guide_rows.append(dict(human_gene=target,sgRNA=sg,condition=condition,feature=label,n_cells=n,model_effect=float(coefs[j,k]),NTC_standardized_difference=float(Y[select,j].mean()),intergenic_reference_difference=float(Y[select,j].mean()-Y[inter,j].mean()) if inter.any() else np.nan,n_NTC=int(ntc.sum()),n_intergenic=int(inter.sum())))
 for gene,z in ob[~ob.NTC & ~ob.INTERGENIC].groupby('target',observed=True):
  sglist=z.sgRNA.unique();valid=[sg for sg in sglist if (z.sgRNA==sg).sum()>=20];jj=[guide_index[sg] for sg in valid];n=len(z)
  for j,label in enumerate(labels):
   eff=float(np.median(coefs[j,jj])) if jj else np.nan;arr=coefs[j,jj];same=int(np.sum(np.sign(arr)==np.sign(eff)));qualified=n>=30 and len(valid)>=2 and same>=2 and same/len(valid)>=2/3 and abs(eff)>=.5
   gene_rows.append(dict(human_gene=gene,condition=condition,feature=label,median_guide_model_effect=eff,unadjusted_NTC_effect=float(Y[ob.target.to_numpy()==gene,j].mean()),n_cells=n,n_NTC=int(ntc.sum()),n_valid_guides=len(valid),guide_direction_concordance=same/max(len(valid),1),phenotype_supported=qualified,reference='verified NO_SITE non-targeting, MOI1 only',model='author ElasticNet alpha0.0005 l1_ratio0.5; complexity+cycle covariates; no EM needed after single-guide restriction; no new P/FDR'))
 sizes.append(dict(condition=condition,n_cells=len(ob),n_NTC=int(ntc.sum()),n_intergenic=int(inter.sum()),n_target_guides=len(guides),model_max_iter_used=int(np.max(model.n_iter_))))
print('RNA and protein effects fitted',flush=True);out(gene_rows,E/'FRANGIEH_PROGRAM_PROTEIN_EFFECTS.tsv');out(guide_rows,E/'FRANGIEH_GUIDE_FEATURE_EFFECTS.tsv.gz');out(sizes,Q/'FRANGIEH_MODEL_SIZES.tsv')
x=pd.DataFrame(gene_rows);phen=[]
for (gene,cond),z in x.groupby(['human_gene','condition'],sort=True):
 supported=z[z.phenotype_supported];features=set(supported.feature);ap=bool(features&{'RNA:ANTIGEN_PRESENTATION','PROTEIN:HLA_A','PROTEIN:HLA_E'});ifng=bool(features&{'RNA:IFNG_RESPONSE','PROTEIN:CD119'});other=bool(features&{'RNA:PDL1','RNA:CELL_STATE_CYCLE','RNA:ADHESION','RNA:STRESS_DEATH','PROTEIN:CD274','PROTEIN:CD58','PROTEIN:CD29','PROTEIN:CD44','PROTEIN:CD49f'})
 cl='AP_IFNG_COMBINED' if ap and ifng else 'AP_AXIS' if ap else 'IFNG_AXIS' if ifng else 'ORTHOGONAL_TO_AP_IFNG' if other else 'UNRESOLVED'
 dic={'human_gene':gene,'condition':cond,'axis_classification':cl,'n_cells':int(z.n_cells.iloc[0]),'n_NTC':int(z.n_NTC.iloc[0]),'n_valid_guides':int(z.n_valid_guides.iloc[0]),'mechanistic_phenotype_supported':cl!='UNRESOLVED','supported_features':';'.join(sorted(features)),'interpretation':'single-cell mechanistic anchor; no independent B16 statistical replication'}
 for n in ['ANTIGEN_PRESENTATION','IFNG_RESPONSE','PDL1','CELL_STATE_CYCLE','ADHESION','STRESS_DEATH']:
  v=z[z.feature=='RNA:'+n];dic[n+'_RNA_effect']=float(v.median_guide_model_effect.iloc[0]) if len(v) else np.nan
 for n in ['HLA_A','CD119','CD274','CD58']:
  v=z[z.feature=='PROTEIN:'+n];dic[n+'_protein_effect']=float(v.median_guide_model_effect.iloc[0]) if len(v) else np.nan
 phen.append(dic)
ph=pd.DataFrame(phen);out(ph[ph.human_gene!='CD58'].merge(fo[['gene','human_gene']],on='human_gene'),T/'FRANGIEH_MECHANISTIC_PHENOTYPES.tsv');out(ph[ph.human_gene=='CD58'].assign(control_status='KNOWN_ORTHOGONAL_POSITIVE_CONTROL'),E/'FRANGIEH_CD58_POSITIVE_CONTROL.tsv')
qc={'processed_RNA_source':'pertpy/sc-pert curation, not claimed exact original authors regulatory matrix','RNA_log1p_CPM_scale':1000000,'source_curation':'https://nbviewer.org/github/theislab/sc-pert/blob/main/datasets/Frangieh_2021_curation.ipynb','n_raw_protein_cells':218331,'n_processed_single_guide_RNA_cells':126966,'n_analyzed_with_NTC_and_intergenic':len(obs),'n_exact_HC_overlap_actual_assignment':len(targets_obs&set(fo.human_gene)),'n_HC_with_any_mechanistic_phenotype':int(ph[(ph.human_gene!='CD58')&ph.mechanistic_phenotype_supported].human_gene.nunique()),'protein_features':len(pnames),'protein_isotype_formula_verified_maxabs':maxdiff,'RNA_selected_features':len(names),'programs':len(pn),'processed_target_annotation_issue':'perturbation_name removes gene digits; all targets rebuilt from exact sgRNA IDs and HGNC approved aliases','no_P_or_FDR_generated':True,'no_cell_as_independent_screen_replication':True};(Q/'FRANGIEH_PHENOTYPE_QC.json').write_text(json.dumps(qc,indent=2));a.file.close();print(qc,flush=True)
