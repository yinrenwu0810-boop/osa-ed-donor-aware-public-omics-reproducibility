import csv, json, math, sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(r"C:\Users\22394\Documents\Codex\2026-07-12\acad"); VP=ROOT/'revision_v2'/'07_virtual_perturbation'; OUT=VP/'10_exploratory_hypoxia_sensitivity'/'attempt_20260910_02'
OUTS=[OUT/'05_target_vs_matched_controls.tsv',OUT/'05_leading_edge_stability.tsv',OUT/'05_leading_edge_overlap_with_ED.tsv']
def rd(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def wr(p,rows):
 if p.exists():raise RuntimeError('refusing overwrite '+p.name)
 with open(p,'x',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['status'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def le(x):return set(filter(None,x.split(';')))
def jac(a,b):return len(a&b)/len(a|b) if a|b else float('nan')
def main():
 if Path.cwd().resolve()!=ROOT.resolve() or any(x.exists() for x in OUTS):raise RuntimeError('wrong root or output exists')
 scope=json.load(open(OUT/'01_scope_freeze.json',encoding='utf-8'))
 if scope['status']!='FROZEN_EXPLORATORY_SPECIFICATION':raise RuntimeError('scope not frozen')
 seed=rd(OUT/'03_target_donor_seed_hypoxia.tsv'); matched=rd(VP/'07_enrichment'/'v1'/'hypoxia_matched_control_enrichment.tsv')
 if len(seed)!=45 or len(matched)!=33:raise RuntimeError('frozen source cardinality mismatch')
 controls=[]
 for t in sorted({r['target_gene'] for r in matched}):
  x=[r for r in matched if r['target_gene']==t]; target=[r for r in x if r['role']=='target'][0]; c=[r for r in x if r['role']=='control']
  tn=float(target['NES']); cn=[float(r['NES']) for r in c]; primary=[r for r in seed if r['target_gene']==t and r['seed']=='2026082701']
  controls.append({'target_gene':t,'condition':target['condition'],'matched_controls':len(c),'comparison_profile':'same_primary_seed_three_donor_consensus','target_NES':tn,'target_FDR_within_Hallmark':target['padj'],'signed_NES_empirical_percentile':sum(v<=tn for v in cn)/len(c),'absolute_NES_empirical_percentile':sum(abs(v)<=abs(tn) for v in cn)/len(c),'conservative_absolute_upper_tail':(sum(abs(v)>=abs(tn) for v in cn)+1)/(len(c)+1),'target_primary_seed_positive_NES_donor_proportion':sum(float(r['NES'])>0 for r in primary)/len(primary),'control_donor_positive_NES_proportion':'NOT_AVAILABLE_PRIMARY_SEED_DONOR_LEVEL_ENRICHMENT_NOT_RECOMPUTED','target_seed_NES_range':max(float(r['NES']) for r in seed if r['target_gene']==t)-min(float(r['NES']) for r in seed if r['target_gene']==t),'control_seed_NES_stability':'NOT_AVAILABLE_ONE_FROZEN_CONTROL_SEED','target_leading_edge_size':len(le(target['leadingEdge'])),'control_leading_edge_size_median':sorted([len(le(r['leadingEdge'])) for r in c])[len(c)//2],'evidence_label':'EXPLORATORY_SAME_SEED_CALIBRATION'})
 stab=[]
 for (t,c,d),rs in defaultdict(list,{}).items():pass
 by=defaultdict(list)
 for r in seed:by[(r['target_gene'],r['condition'],r['donor'])].append(r)
 for k,rs in sorted(by.items()):
  for i in range(len(rs)):
   for j in range(i+1,len(rs)):
    a,b=le(rs[i]['leading_edge']),le(rs[j]['leading_edge']);stab.append({'comparison_type':'within_target_donor_seed','target_gene':k[0],'condition':k[1],'comparison_id_a':k[2]+'|'+rs[i]['seed'],'comparison_id_b':k[2]+'|'+rs[j]['seed'],'leading_edge_jaccard':jac(a,b),'leading_edge_size_a':len(a),'leading_edge_size_b':len(b),'status':'PASS'})
 formal=rd(VP/'07_enrichment'/'v1'/'hypoxia_donor_enrichment.tsv'); primary=rd(VP/'07_enrichment'/'v1'/'hypoxia_primary_summary.tsv')
 for t in sorted({r['target_gene'] for r in formal}):
  rs=[r for r in formal if r['target_gene']==t]
  for i in range(3):
   for j in range(i+1,3):
    a,b=le(rs[i]['leadingEdge']),le(rs[j]['leadingEdge']);stab.append({'comparison_type':'between_donor_five_seed_consensus','target_gene':t,'condition':rs[i]['condition'],'comparison_id_a':rs[i]['donor'],'comparison_id_b':rs[j]['donor'],'leading_edge_jaccard':jac(a,b),'leading_edge_size_a':len(a),'leading_edge_size_b':len(b),'status':'PASS'})
 for i in range(3):
  for j in range(i+1,3):
   a,b=le(primary[i]['leadingEdge']),le(primary[j]['leadingEdge']);stab.append({'comparison_type':'between_target_formal_consensus','target_gene':primary[i]['target_gene']+'|'+primary[j]['target_gene'],'condition':'NOT_APPLICABLE','comparison_id_a':primary[i]['target_gene'],'comparison_id_b':primary[j]['target_gene'],'leading_edge_jaccard':jac(a,b),'leading_edge_size_a':len(a),'leading_edge_size_b':len(b),'status':'PASS'})
 ed=defaultdict(set)
 for r in rd(VP/'07_enrichment'/'v1'/'ed_signature_sets.tsv'):ed[r['signature']].add(r['gene'])
 overlap=[]
 for r in primary:
  a=le(r['leadingEdge'])
  for name,genes in ed.items():overlap.append({'target_gene':r['target_gene'],'reference_set':name,'reference_status':'FROZEN_VP_G06','overlap_size':len(a&genes),'overlap_genes':';'.join(sorted(a&genes)),'not_available_reason':''})
  overlap.append({'target_gene':r['target_gene'],'reference_set':'STRICT_OSA_CPAP_HYPOXIA_LEADING_EDGE','reference_status':'NOT_AVAILABLE','overlap_size':'','overlap_genes':'','not_available_reason':'No frozen traceable gene-level OSA/CPAP bridge set exists.'})
 wr(OUTS[0],controls);wr(OUTS[1],stab);wr(OUTS[2],overlap)
 print(json.dumps({'gate':'VP-ES01-05','status':'PASS_COMPUTED_PENDING_NEXT_GATE','control_rows':len(controls),'leading_edge_rows':len(stab),'overlap_rows':len(overlap)}))
if __name__=='__main__':main()
