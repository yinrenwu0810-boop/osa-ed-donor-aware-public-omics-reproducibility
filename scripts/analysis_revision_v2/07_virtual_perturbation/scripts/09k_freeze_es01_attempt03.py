import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
root=Path(r'C:\Users\22394\Documents\Codex\2026-07-12\acad'); vp=root/'revision_v2'/'07_virtual_perturbation'; out=vp/'10_exploratory_hypoxia_sensitivity'/'attempt_20260910_03'
pre=out/'00_preflight.json'; target=out/'01_scope_freeze.json'
if Path.cwd().resolve()!=root.resolve() or target.exists():raise SystemExit('wrong root or output exists')
p=json.loads(pre.read_text(encoding='utf-8'))
if p['status']!='PASS_INPUTS_READ_ONLY' or p['errors']:raise SystemExit('preflight not pass')
h=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
scope={'gate':'VP-ES01-01','attempt':'attempt_20260910_03','status':'FROZEN_CORRECTED_ES06_SPECIFICATION','created_at_utc':datetime.now(timezone.utc).isoformat(),'preflight_sha256':h(pre),'analysis_label':'EXPLORATORY_SENSITIVITY_NOT_CONFIRMATORY','retained_attempts':['attempt_20260910_01','attempt_20260910_02'],'formal_result_protection':'VP-G06 HALLMARK_HYPOXIA formal results remain unchanged.','fixed_specifications':['S0 distance plus five-seed median standardized rank plus three-donor consensus','S1 distance plus five-seed mean standardized rank plus three-donor mean rank','S2 Z descending plus S0 aggregation','S3 donor-first descriptive NES aggregation'],'required_es06_corrections':['HALLMARK_HYPOXIA mapping QC includes source size, universe overlap, missing symbols, leading-edge network degree and rank distributions.','S0 is independently compared numerically with VP-G06 formal baseline.','Specification concordance includes NES difference, direction concordance, and leading-edge Jaccard.'],'fdr_rule':'Same VP-G06 tested Hallmark family; no single-set or reduced-family FDR.','biological_unit':'donor; seeds are technical network-stability units only.','next_gate':'VP-ES01-06 corrected implementation only'}
target.write_text(json.dumps(scope,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
