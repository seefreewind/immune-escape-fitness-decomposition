"""Audit frozen evidence roles without recomputing or pooling external significance."""
import csv,json
from pathlib import Path
def audit(package,output):
 with (package/'data/CANDIDATE_EVIDENCE_MATRIX.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 result=dict(HC=len(rows),Frangieh_exact=sum(r['Frangieh_exact_overlap']=='True' for r in rows),Frangieh_phenotype=sum(r['Frangieh_mechanistic_support']=='True' for r in rows),Joung_role='CRISPRa directional benchmark, not direct KO replication',Yeh_role='targeted descriptive knockout support',Akana_role='CRISPRa context',Frangieh_role='mechanistic phenotype anchor',Tier_A=sum(r['tier']=='TIER_A' for r in rows))
 Path(output).write_text(json.dumps(result,indent=2));return result
