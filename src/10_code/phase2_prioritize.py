from pathlib import Path
import pandas as pd,numpy as np,json
R=Path(__file__).resolve().parents[1];T=R/'07_tables';S=T/'phase2/reconciled';E=T/'phase2'
def rd(p):return pd.read_csv(p,sep='\t',low_memory=False)
def out(x,p):x.to_csv(p,sep='\t',index=False)
h=rd(S/'KNOWN_MECHANISM_CLASSIFICATION.tsv');m=rd(S/'FUNCTIONAL_MODULES.tsv');a=rd(S/'phase2/MODULE_GENE_ASSIGNMENTS.tsv');ph=rd(T/'FRANGIEH_MECHANISTIC_PHENOTYPES.tsv');y=rd(T/'YEH_DIRECTIONAL_SUPPORT.tsv');b=rd(E/'FRANGIEH_BULK_VIABILITY_DESCRIPTIVE.tsv')
f=ph.groupby('gene').mechanistic_phenotype_supported.any();h['Frangieh_mechanistic_support']=h.gene.map(f).fillna(False).astype(bool);h['Yeh_T_directional_support']=h.gene.map(y.set_index('gene').T_support).fillna(False).astype(bool);h['Frangieh_bulk_directional_descriptive']=h.gene.map(b.set_index('gene').all_ratios_direction_consistent).fillna(False).astype(bool)
h['independent_directional_sources']=h.Yeh_T_directional_support.astype(int)+h.Frangieh_bulk_directional_descriptive.astype(int);h['evidence_types']=h.independent_directional_sources+h.Frangieh_mechanistic_support.astype(int)
a=a.merge(m[['module','stable_coherent','median_bootstrap_Jaccard','coherent_GO_term','coherent_GO_FDR']],on='module',validate='many_to_one')
ra=a.query("scope=='RESIDUAL_ONLY'").set_index('gene');ha=a.query("scope=='ALL_HC'").set_index('gene')
for prefix,d in [('residual_module',ra),('allHC_module',ha)]:
 for col in ['module','stable_coherent','degree','median_bootstrap_Jaccard','coherent_GO_term','coherent_GO_FDR']:
  h[prefix+'_'+col]=h.gene.map(d[col])
h['stable_coherent_residual_module']=h.residual_module_stable_coherent.fillna(False).astype(bool);h['mechanistic_coherence']=h.stable_coherent_residual_module|h.allHC_module_stable_coherent.fillna(False).astype(bool)|h.Frangieh_mechanistic_support
h['residual_level_pre_audit']=np.where(h.residual,np.where((h.independent_directional_sources>0)|h.Frangieh_mechanistic_support|h.stable_coherent_residual_module,'RESIDUAL_LEVEL_2','RESIDUAL_LEVEL_1'),'NOT_RESIDUAL');h['selected_shortlist']=False;h['selection_reason']=''
selected=[]
for mid,z in h[h.residual&h.stable_coherent_residual_module].groupby('residual_module_module'):
 z=z.sort_values(['evidence_types','residual_module_median_bootstrap_Jaccard','residual_module_degree','gene'],ascending=[False,False,False,True]).head(2)
 for g in z.gene:selected.append(g);h.loc[h.gene==g,'selection_reason']='two representatives per stable coherent residual module, evidence types then stability then degree then symbol'
extras=h[h.residual&(h.independent_directional_sources>0)&~h.gene.isin(selected)].sort_values(['evidence_types','gene'],ascending=[False,True])
for g in extras.gene:
 if len(selected)>=30:break
 selected.append(g);h.loc[h.gene==g,'selection_reason']='additional residual with independent experimental directional support'
h.loc[h.gene.isin(selected),'selected_shortlist']=True
assert len(selected)<=30;out(h,E/'PHASE2_PRE_AUDIT_EVIDENCE_MATRIX.tsv');out(h[h.selected_shortlist],E/'PHASE2_NOVELTY_SHORTLIST.tsv');print(h[h.selected_shortlist][['gene','human_gene','arm','independent_directional_sources','residual_module_module','residual_module_degree']].to_string(index=False));print('shortlist',len(selected))
