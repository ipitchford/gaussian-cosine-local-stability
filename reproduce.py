"""Fresh-copy numerical replay; no optimizer trajectory regeneration."""
from pathlib import Path
import sys,os,shutil,subprocess,hashlib,json,time
import numpy,scipy,sympy
root=Path(__file__).resolve().parent
if len(sys.argv)!=2:raise SystemExit('Usage: python reproduce.py NEW_OUTPUT_DIRECTORY')
out=Path(sys.argv[1]).resolve();out.mkdir()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={str(p.relative_to(root)):sha(p) for p in (root/'code').iterdir() if p.is_file()}
report={'scope':'Finite diagnostics, not analytic proof or global certification','python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,'sympy':sympy.__version__,'sourceHashes':sources,'runs':[]}
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
for mode,flags in [('normal',[]),('optimized',['-O'])]:
    copy=out/mode;shutil.copytree(root/'code',copy)
    for script in ['check_local_algebra.py','check_results.py','test_guards.py']:
        cmd=[sys.executable,'-B',*flags,str(copy/script)];start=time.monotonic()
        p=subprocess.run(cmd,env=env,capture_output=True,timeout=120)
        (copy/(script+'.stdout')).write_bytes(p.stdout);(copy/(script+'.stderr')).write_bytes(p.stderr)
        report['runs'].append({'mode':mode,'command':cmd,'exit':p.returncode,'elapsed':time.monotonic()-start,'stdoutSha256':hashlib.sha256(p.stdout).hexdigest(),'stderrSha256':hashlib.sha256(p.stderr).hexdigest()})
report['sourceUnchanged']=all(sha(root/n)==h for n,h in sources.items())
report['status']='passed' if report['sourceUnchanged'] and all(x['exit']==0 for x in report['runs']) else 'failed'
(out/'receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(report['status'])
raise SystemExit(0 if report['status']=='passed' else 1)
