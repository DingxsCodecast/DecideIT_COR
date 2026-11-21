import json
from pathlib import Path
import subprocess
import argparse
import os

root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--bin',type=Path,default=root/'bin/oracle_driver')
parser.add_argument('--allow-controlled-case',type=int,action='append',default=[],
                    help='case allowed to terminate without claiming an optimum')
args=parser.parse_args()
allowed_controlled=set(args.allow_controlled_case)
expected=json.loads((root/'tests/oracle_expected.json').read_text())
results=[]
for ref in expected:
    case=root/'tests/oracle_cases'/('%03d.txt'%ref['case'])
    p=subprocess.run([str(args.bin),str(case)],capture_output=True,text=True,timeout=10)
    lines=[s for s in p.stdout.splitlines() if s.startswith('AUDIT ')]
    r={'case':ref['case'],'returncode':p.returncode,'passed':False,'error':p.stderr}
    if p.returncode==0 and lines:
        v=list(map(float,lines[-1].split()[1:]))
        r.update(status=v[0],cuts=v[5:7])
        if v[0] == 0:
            r['error_max']=max(abs(v[1]-ref['expected'][0]),abs(v[2]-ref['expected'][1]))
            r['passed']=r['error_max']<=1e-6
        elif ref['case'] in allowed_controlled and v[0] == 20:
            r['controlled']=True
            r['passed']=True
    results.append(r)
results_dir=Path(os.environ.get('DECIDEIT_TEST_RESULTS_DIR',root/'test_results'))
results_dir.mkdir(parents=True,exist_ok=True)
(results_dir/'oracle.json').write_text(json.dumps(results,indent=2))
certified=[r for r in results if r.get('status') == 0]
controlled=[r for r in results if r.get('controlled')]
print('ORACLE',len(results),'CERTIFIED',len(certified),'CONTROLLED',len(controlled),
      'MAX_ERROR',max((r.get('error_max',0.0) for r in certified),default=0.0))
for r in results:
    if not r['passed']:print(r)
raise SystemExit(0 if all(r['passed'] for r in results) else 1)
