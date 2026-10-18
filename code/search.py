import numpy as np, json, time, sys
from scipy.optimize import minimize
from objective import objective
from pathlib import Path
from validation import finite,require,NumericalDomainError
ROOT=Path(__file__).resolve().parent

def frame(x,d):
    require(d in (3,4),'dimension must be 3 or 4')
    x=np.asarray(x,dtype=float)
    require(x.shape==((5,) if d==3 else (6,)),'Wrong angular parameter shape')
    finite(x,'parameters')
    if d==3:
        return np.array([[1,0,0],[np.cos(x[0]),np.sin(x[0]),0],[np.cos(x[1]),np.sin(x[1])*np.cos(x[2]),np.sin(x[1])*np.sin(x[2])],[np.cos(x[3]),np.sin(x[3])*np.cos(x[4]),np.sin(x[3])*np.sin(x[4])]])
    return np.array([[1,0,0,0],[np.cos(x[0]),np.sin(x[0]),0,0],[np.cos(x[1]),np.sin(x[1])*np.cos(x[2]),np.sin(x[1])*np.sin(x[2]),0],[np.cos(x[3]),np.sin(x[3])*np.cos(x[4]),np.sin(x[3])*np.sin(x[4])*np.cos(x[5]),np.sin(x[3])*np.sin(x[4])*np.sin(x[5])]])
if __name__=='__main__':
    d=int(sys.argv[1]); start=int(sys.argv[2]); count=int(sys.argv[3])
    rng=np.random.default_rng(21402026+1000*d+start); results=[]
    for k in range(start,start+count):
        x0=rng.uniform(.15,np.pi-.15,5 if d==3 else 6)
        def f(x):
            try:val=objective(frame(x,d),False)
            except NumericalDomainError:val=float('nan')
            return val if np.isfinite(val) else 2.
        r=minimize(f,x0,method='BFGS', options={'gtol':2e-6,'maxiter':150})
        U=frame(r.x,d)
        rec={'d':d,'start':k,'value':float(r.fun),'success':bool(r.success),'parameters':r.x.tolist(),'gram':(U@U.T).tolist(),'singular_values':np.linalg.svd(U,compute_uv=False).tolist(),'nfev':int(r.nfev)}
        results.append(rec)
        (ROOT/f'search_{d}_{start}.json').write_text(json.dumps(results,indent=2))
        print(d,k,r.fun,r.success,'sv',rec['singular_values'],flush=True)
