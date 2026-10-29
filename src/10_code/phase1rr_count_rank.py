from pathlib import Path
import pandas as pd,numpy as np,json
from scipy import stats
R=Path(__file__).resolve().parents[1];T=R/'07_tables/phase1rr';Q=R/'06_qc/phase1rr';out=[]
for s in json.loads((R/'01_analysis_freeze/CONTRAST_SPECIFICATIONS.json').read_text()):
 if s['role']!='PRIMARY':continue
 cid=s['contrast'];x=pd.read_csv(R/s['raw'],sep='\t');v=x[s['samples']].to_numpy();keep=(v.sum(1)>=60)&((v>=10).sum(1)>=3)&(v[:,:3].sum(1)>=30);n=pd.read_csv(Q/(cid+'.TMM.tsv'),sep='\t');exposure=n['lib.size'].to_numpy()*n['norm.factors'].to_numpy();z=v/exposure*1e7;eff=np.log2((z[:,3:].mean(1)+5)/(z[:,:3].mean(1)+5));g=pd.DataFrame({'gene':x.Gene[keep],'count_effect':eff[keep]});g['pr']=(-g.count_effect).rank(method='average')/(len(g)+1);g['nr']=g.count_effect.rank(method='average')/(len(g)+1);rows=[]
 for name,a in g.groupby('gene'):
  if len(a)<3 or name in ['LacZ','EGFP','luciferase']:continue
  k=np.arange(1,len(a)+1);p=min(1,float(stats.beta.cdf(np.sort(a.pr),k,len(a)-k+1).min()*len(a)));neg=min(1,float(stats.beta.cdf(np.sort(a.nr),k,len(a)-k+1).min()*len(a)));rows.append({'gene':name,'COUNT_RANK_SCORE':-np.log10(max(p,1e-300))+np.log10(max(neg,1e-300)),'median_count_log2FC':a.count_effect.median(),'role':'independent count-rank sensitivity; no QL variance or null P; no threshold training'})
 a=pd.DataFrame(rows).set_index('gene');a.to_csv(T/(cid+'.COUNT_RANK_SENSITIVITY.tsv'),sep='\t');primary=pd.read_csv(T/(cid+'.OBSERVED_GENE_SCORES.tsv'),sep='\t').set_index('gene');j=a.join(primary.PRIMARY_GENE_SCORE,how='inner');out.append({'contrast':cid,'n':len(j),'spearman':stats.spearmanr(j.COUNT_RANK_SCORE,j.PRIMARY_GENE_SCORE).statistic,'pearson':stats.pearsonr(j.COUNT_RANK_SCORE,j.PRIMARY_GENE_SCORE).statistic,'sign_concordance':(j.COUNT_RANK_SCORE*j.PRIMARY_GENE_SCORE>0).mean(),'top100_overlap':len(set(j.COUNT_RANK_SCORE.nlargest(100).index)&set(j.PRIMARY_GENE_SCORE.nlargest(100).index)),'bottom100_overlap':len(set(j.COUNT_RANK_SCORE.nsmallest(100).index)&set(j.PRIMARY_GENE_SCORE.nsmallest(100).index))})
pd.DataFrame(out).to_csv(T/'INDEPENDENT_COUNT_RANK_AGREEMENT.tsv',sep='\t',index=False);print('Independent count-based rank sensitivity completed')
