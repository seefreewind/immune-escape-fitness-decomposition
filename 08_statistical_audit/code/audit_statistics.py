from pathlib import Path
from fractions import Fraction
import json,math
import numpy as np,pandas as pd
from scipy.stats import binom
P=Path(__file__).resolve().parents[2];R=P.parent;Q=P/'08_QC_Reports';D=Q/'audit_data';m=18710
save=lambda x,n: pd.DataFrame(x).to_csv(Q/n,sep='\t',index=False)
def bh(p):
 p=np.asarray(p);o=np.argsort(p,kind='stable');v=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[o]=np.minimum(v,1);return q
def matrix(n,B):return np.fromfile(D/(n+'.bin'),dtype='<f8').reshape((m,B),order='F')
genes=pd.read_csv(D/'TRAIN_genes.tsv',sep='\t').gene.to_numpy()
tr=matrix('TRAIN_score',1000);ho=matrix('HELDOUT_score',100);sp=matrix('SPIKE_score',60)
assert np.array_equal(genes,pd.read_csv(D/'HELDOUT_genes.tsv',sep='\t').gene)
ss={x:set(pd.read_csv(D/(x+'_seeds.tsv'),sep='\t').seed) for x in ['TRAIN','HELDOUT','SPIKE']}
assert not(ss['TRAIN']&ss['HELDOUT'] or ss['TRAIN']&ss['SPIKE'] or ss['HELDOUT']&ss['SPIKE'])
# Exact integer arithmetic for all 1001 attainable P values.
rows=[]
for j in range(1,1002):
 k=(m*10*j+1000)//1001
 rows.append({'grid_numerator':j,'P_exact':f'{j}/1001','P':j/1001,'first_BH_rank_at_q010':k,'feasible_rank':k<=m,'expected_count_at_or_below_under_discrete_uniform':m*j/1001,'expected_count_in_single_bin':m/1001})
save(rows,'P_RESOLUTION_BH_ANALYTICAL_AUDIT.tsv')
obs=pd.read_csv(R/'07_tables/phase1rr/mTKO_B16_T17.EMPIRICAL_INFERENCE.tsv',sep='\t').set_index('gene').loc[genes]
inf=pd.read_csv(R/'07_tables/phase1rr/PHASE1RR_GENE_INFERENCE.tsv',sep='\t');inf=inf[inf.contrast=='mTKO_B16_T17'].set_index('gene').loc[genes]
pcalc=(1+(np.abs(tr)>=np.abs(obs.PRIMARY_GENE_SCORE.to_numpy())[:,None]).sum(axis=1))/1001
qcalc=bh(pcalc)
assert np.max(np.abs(pcalc-obs.empirical_P.to_numpy()))<1e-12
assert np.max(np.abs(qcalc-obs.empirical_FDR.to_numpy()))<1e-12
assert np.sum(qcalc<=.1)==3999
hc=inf.HIGH_CONFIDENCE.astype(str).str.lower().eq('true').to_numpy();effect=(qcalc<=.1)&(np.abs(obs.median_log2FC.to_numpy())>=.5)
assert hc.sum()==1492 and effect.sum()==2675
bins=np.rint(pcalc*1001).astype(int)
rr=[]
for j in range(1,1002):
 z=bins==j
 rr.append({'grid_numerator':j,'P':j/1001,'n_eligible':int(z.sum()),'n_FDR':int((z&(qcalc<=.1)).sum()),'n_FDR_effect':int((z&effect).sum()),'n_HC':int((z&hc).sum()),'cumulative_n':int((bins<=j).sum()),'BH_critical_rank':(m*10*j+1000)//1001,'abs_score_min':float(np.abs(obs.PRIMARY_GENE_SCORE.to_numpy())[z].min()) if z.any() else np.nan,'abs_score_max':float(np.abs(obs.PRIMARY_GENE_SCORE.to_numpy())[z].max()) if z.any() else np.nan})
save(rr,'OBSERVED_P_GRID_DIAGNOSTICS.tsv')
o=np.argsort(pcalc,kind='stable');ks=np.arange(1,m+1);valid=pcalc[o]<=ks*.1/m;k=int(ks[valid].max()) if valid.any() else 0
save([{'rank':int(i+1),'gene':genes[o[i]],'P':pcalc[o[i]],'BH_critical_P':(i+1)*.1/m,'BH_adjusted':qcalc[o[i]],'rejected':qcalc[o[i]]<=.1,'grid_numerator':bins[o[i]]} for i in range(max(0,k-20),min(m,k+20))],'BH_REJECTION_BOUNDARY.tsv')
# Bounded first-B stored-reference sensitivity, no fitting or new simulated outcomes.
sens=[];nullrows=[];spikerows=[];common=[]
spgenes=pd.read_csv(R/'06_qc/phase1rr/mTKO_B16_T17.PRESPECIFIED_SPIKE_GENES.tsv',sep='\t');lookup=dict(zip(genes,range(m)));idx=np.array([lookup[g] for g in spgenes.gene]);isnull=np.ones(m,bool);isnull[idx]=False
Pmat={}
for B in [250,500,1000]:
 sortedtr=np.sort(np.abs(tr[:,:B]),axis=1)
 for name,mat in [('OBSERVED',obs.PRIMARY_GENE_SCORE.to_numpy()[:,None]),('HELDOUT',ho),('SPIKE',sp)]:
  pp=np.empty(mat.shape)
  for g in range(m):pp[g]=(1+B-np.searchsorted(sortedtr[g],np.abs(mat[g]),side='left'))/(B+1)
  Pmat[(B,name)]=pp
  ns=[int((bh(pp[:,b])<=.1).sum()) for b in range(pp.shape[1])]
  sens.append({'B':B,'stage':name,'screens':pp.shape[1],'mean_BH_discoveries':np.mean(ns),'screens_with_BH_discovery':sum(n>0 for n in ns),'min_BH_discoveries':min(ns),'max_BH_discoveries':max(ns),'scope':'sensitivity only; first B archived training draws; not replacement inference'})
  if B==1000 and name=='HELDOUT':
   saved=matrix('HELDOUT_P_A',100);assert np.max(np.abs(pp-saved))<1e-12
   for b in range(100):nullrows.append({'simulation':b+1,'p05':np.mean(pp[:,b]<=.05),'p01':np.mean(pp[:,b]<=.01),'BH_n':ns[b],'min_P_n':int((pp[:,b]==1/1001).sum()),'first_BH_rank_needed':187})
  if B==1000 and name=='SPIKE':
   saved=matrix('SPIKE_P_A',60);assert np.max(np.abs(pp-saved))<1e-12
   for b in range(60):
    q=bh(pp[:,b]);sel=q<=.1
    spikerows.append({'simulation':b+1,'FDP':sum(sel&isnull)/max(1,sum(sel)),'null_nominal_P05':np.mean(pp[isnull,b]<=.05),'null_nominal_P01':np.mean(pp[isnull,b]<=.01),'null_minP_n':sum((pp[:,b]==1/1001)&isnull),'rejections':int(sel.sum()),'false_rejections':int(sum(sel&isnull))})
# Round comparator B upward onto the identical P grid; no randomization or mid-P.
for stage,N in [('HELDOUT',100),('SPIKE',60)]:
 pb=matrix(stage+'_P_B',N)
 for b in range(N):
  pg=np.maximum(1,np.ceil(pb[:,b]*1001-1e-12))/1001
  for pipe,pp in [('A_original',Pmat[(1000,stage)][:,b]),('B_original',pb[:,b]),('B_rounded_up_1001_grid',pg)]:
   sel=bh(pp)<=.1
   common.append({'stage':stage,'simulation':b+1,'pipeline':pipe,'BH_n':int(sel.sum()),'FDP':sum(sel&isnull)/max(1,sum(sel)) if stage=='SPIKE' else np.nan,'P05':np.mean(pp<=.05)})
save(sens,'FINITE_REFERENCE_SENSITIVITY.tsv');save(nullrows,'GLOBAL_NULL_RECOMPUTED.tsv');save(spikerows,'SPIKE_NULL_RECOMPUTED.tsv');save(common,'COMMON_P_GRID_DIAGNOSTIC.tsv')
# Cross-gene signed score dependence over independent held-out simulation draws.
# A deterministic, evenly spaced 1000-gene diagnostic subset avoids an enormous 18,710-square matrix.
ix=np.linspace(0,m-1,1000,dtype=int);corr=np.corrcoef(ho[ix]);cv=corr[np.triu_indices(1000,1)]
dep={'genes':1000,'simulation_draws':100,'selection':'evenly spaced indices of frozen alphabetical gene order','mean_correlation':float(cv.mean()),'median_correlation':float(np.median(cv)),'q01':float(np.quantile(cv,.01)),'q99':float(np.quantile(cv,.99)),'fraction_abs_corr_gt_02':float(np.mean(np.abs(cv)>.2)),'interpretation':'descriptive signed-score dependence diagnostic; does not establish BH dependence assumptions'}
(Q/'DEPENDENCE_DIAGNOSTIC.json').write_text(json.dumps(dep,indent=2))
smalln=int((pcalc==1/1001).sum());boundary=pcalc[o[k-1]]
summary={'m':m,'training_B':1000,'P_min':1/1001,'first_BH_rank':187,'null_expected_minP_count':m/1001,'independent_binomial_minP_ge187_probability':float(binom.sf(186,m,1/1001)),'observed_minP_n':smalln,'observed_minP_fraction':smalln/m,'FDR_minP_n':int(sum((qcalc<=.1)&(pcalc==1/1001))),'observed_BH_n':k,'boundary_P':boundary,'boundary_grid_numerator':round(boundary*1001),'training_heldout_spike_seeds_disjoint':True,'observed_P_max_reproduction_error':float(np.max(np.abs(pcalc-obs.empirical_P.to_numpy()))),'observed_Q_max_reproduction_error':float(np.max(np.abs(qcalc-obs.empirical_FDR.to_numpy()))),'HC_n':int(hc.sum()),'effect_FDR_n':int(effect.sum()),'null_minP_count_min':min(r['min_P_n'] for r in nullrows),'null_minP_count_max':max(r['min_P_n'] for r in nullrows),'all_reference_draws_reused_across_heldout':'shared frozen training bank; heldout random seeds disjoint; generator shared','common_grid_summary':pd.DataFrame(common).groupby(['stage','pipeline']).agg(BH_any_fraction=('BH_n',lambda s:float(np.mean(s>0))),mean_BH_n=('BH_n','mean')).reset_index().to_dict('records'),'sensitivity':sens,'null_summary':pd.DataFrame(nullrows).describe().to_dict(),'spike_summary':pd.DataFrame(spikerows).describe().to_dict()}
(Q/'STATISTICAL_COMPUTATION_SUMMARY.json').write_text(json.dumps(summary,indent=2))
# Explicitly requested vector diagnostic PDF, with source PNG for Word embedding.
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(6.69,6),layout='constrained')
axs[0,0].bar(np.arange(1,1002),np.bincount(bins,minlength=1002)[1:],color='#397CA1');axs[0,0].set(xlabel='Empirical P grid numerator j (P=j/1001)',ylabel='Observed genes')
axs[0,1].plot(ks,pcalc[o],color='#397CA1');axs[0,1].plot(ks,ks*.1/m,'--',color='#C17B4E');axs[0,1].axvline(k,color='gray',ls=':');axs[0,1].set(xlabel='Sorted gene rank',ylabel='P and BH critical P',xlim=(1,6000),ylim=(0,.05))
axs[1,0].hist([r['min_P_n'] for r in nullrows],color='#368D87');axs[1,0].axvline(187,color='#C17B4E',ls='--');axs[1,0].set(xlabel='Minimum-P genes per held-out null screen',ylabel='Simulated screens')
a=pd.DataFrame(common);g=a[a.stage=='HELDOUT'].groupby('pipeline',sort=False).BH_n.apply(lambda s:np.mean(s>0));axs[1,1].bar(range(len(g)),g,color=['#397CA1','#C17B4E','#AEB6BA']);axs[1,1].set_xticks(range(len(g)),['A','B','B rounded\nto 1001 grid']);axs[1,1].set(ylabel='Null screens with ≥1 BH discovery',ylim=(0,1));
for ax,c in zip(axs.flat,'abcd'):ax.text(-.15,1.04,c,transform=ax.transAxes,weight='bold');ax.spines[['top','right']].set_visible(False)
fig.savefig(Q/'P_GRID_DIAGNOSTIC_FIGURE.pdf');fig.savefig(Q/'P_GRID_DIAGNOSTIC_FIGURE.png',dpi=600);plt.close(fig)
print(json.dumps({k:summary[k] for k in ['observed_minP_n','observed_BH_n','boundary_P','common_grid_summary','sensitivity']},indent=2))
