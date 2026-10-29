"""Reconstruct locked HC/tier labels using rerun scores and frozen ancillary evidence.
Ancillary baseline/temporal/author and module/external evidence are explicit inputs,
not regenerated from restricted third-party annotation by this lightweight runner.
"""
import pandas as pd,numpy as np
def classify(inference,robustness,ancillary,evidence,output):
 x=pd.read_csv(inference,sep='\t').set_index('gene'); a=pd.read_csv(ancillary,sep='\t').set_index('gene'); r=pd.read_csv(robustness,sep='\t').set_index('gene'); e=pd.read_csv(evidence,sep='\t').set_index('gene')
 if set(x.index)!=set(a.index):raise ValueError('Eligible gene universe changed; no silent adaptation')
 flags=['temporal_compatible','author_direction_compatible','rank_direction_compatible','primary_score_effect_sign_compatible','baseline_specificity_compatible','count_rank_direction_compatible']
 keep=(x.empirical_FDR<=.1)&(abs(x.median_log2FC)>=.5)&a[flags].all(axis=1)&r.sample_robust&r.guide_robust
 x['HIGH_CONFIDENCE']=keep;x['ancillary_provenance']='FROZEN_CASE_STUDY_BASELINE_TEMPORAL_AUTHOR_INPUTS';x['tier']=None;x['residual']=False
 selected=x.index[keep]
 if not set(selected)<=set(e.index):raise ValueError('New candidate outside frozen evidence; downstream rerun required')
 z=e.loc[selected]
 x.loc[selected,'residual']=~z.known_core & ~z.other_published_hit & ~z.dominant_known_member
 b=z.allHC_module_stable_coherent|z.stable_coherent_residual_module|z.Frangieh_mechanistic_support
 aa=(z.residual_level=='RESIDUAL_LEVEL_3')&z.independent_experimental_directional_support&z.mechanistic_coherence&z.negative_audit_sufficient&~z.direct_prior_immune_role
 x.loc[selected,'tier']=np.where(aa,'TIER_A',np.where(b,'TIER_B','TIER_C'))
 x.to_csv(output,sep='\t');return x
