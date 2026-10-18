"""Reproduce algebraic, numerical and search audits. Not a global proof checker."""
from __future__ import annotations
import json,glob,platform,time
from pathlib import Path
import numpy as np
import scipy
import sympy as sp
from scipy.stats import qmc,norm
from objective import objective as _objective
from intrinsic_objective import intrinsic_objective as _intrinsic_objective
from validation import finite,require
from search import frame

ROOT=Path(__file__).resolve().parent

def objective(U):return float(finite(_objective(U),'elementary result'))
def intrinsic_objective(U):return float(finite(_intrinsic_objective(U),'intrinsic result'))

def load_batches(d):
    records=[]
    for start,count in [(0,10),(10,90)]:
        p=ROOT/f'search_{d}_{start}.json'
        require(p.is_file(),'Missing required batch '+p.name)
        batch=json.loads(p.read_text())
        require(isinstance(batch,list) and len(batch)==count,'Incorrect batch count')
        require([r['start'] for r in batch]==list(range(start,start+count)),'Incorrect batch start coverage')
        for r in batch:
            require(r['d']==d and isinstance(r['success'],bool),'Invalid record dimension/status')
            for k in ['value','parameters','gram','singular_values']:finite(r[k],k)
            require(np.shape(r['gram'])==(4,4) and np.shape(r['singular_values'])==(d,),'Wrong stored shape')
            U=frame(np.asarray(r['parameters']),d)
            require(np.max(np.abs(U@U.T-np.asarray(r['gram'])))<1e-10,'Stored Gram mismatch')
            require(np.max(np.abs(np.linalg.svd(U,compute_uv=False)-np.asarray(r['singular_values'])))<1e-10,'Stored singular-value mismatch')
        records+=batch
    return records

def main():
    report={'status':'PARTIAL: global inequality not proved','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__}
    r=sp.sqrt(2);c=r/2;t=r-1
    V=sp.Matrix([[1,0],[r/2,r/2],[0,1],[-r/2,r/2]])
    A=sp.Matrix([[0,1,0,-1],[1,0,1,0],[0,1,0,1],[-1,0,1,0]])
    P=sp.eye(4)-V*V.T/2
    def iszero(M):return all(sp.simplify(x)==0 for x in M)
    checks={'frame_columns':iszero(V.T*V-2*sp.eye(2)),
      'projection_idempotent':iszero(P*P-P),'projection_symmetric':iszero(P-P.T),
      'signed_cycle_square':iszero(A*A-2*sp.eye(4)),
      'transverse_hessian_identity':iszero(t*(c*sp.eye(4)-A/2)-(2-r)*P),
      'rotation_kernel':iszero(P*V)}
    require(all(checks.values()),'Symbolic identity failure');report['exact_symbolic_checks']=checks
    f0=1-2*np.sqrt(2)/np.pi
    report['closed_constants']={
      'cosine':str(sp.N(1-2*r/sp.pi,40)),
      'planar_triple_plus_independent':str(sp.N(1-sp.Rational(3,2)/sp.pi*(4+sp.sqrt(3)-sp.sqrt(15)),40)),
      'four_independent':str(sp.N(1-(12-16/sp.sqrt(3))/sp.pi,40)),
      'transverse_quadratic_coefficient':str(sp.N((2-r)/sp.pi,40))}
    rng=np.random.default_rng(9926423);cross=[]
    for d in (3,4):
      errors=[]
      for k in range(100):
        U=rng.normal(size=(4,d));x=objective(U);y=intrinsic_objective(U)
        errors.append(float(finite(abs(x-y),'cross discrepancy')))
      require(len(errors)==100 and max(errors)<1e-10,'Cross discrepancy failure')
      cross.append({'rank':d,'count':100,'maximum_absolute_discrepancy':float(max(errors))})
    report['independent_evaluators']=cross
    Vn=np.array(V,dtype=float);Pn=np.array(P,dtype=float);theta=np.arange(4)*np.pi/4
    hess=[]
    rng=np.random.default_rng(8102026)
    for d in (3,4):
      for k in range(4):
        a=rng.normal(size=4);a-=a.mean();B=Pn@rng.normal(size=(4,d-2))
        a/=np.linalg.norm(a);B/=np.linalg.norm(B)
        expected=np.sqrt(2)/(4*np.pi)*sum((np.roll(a,-1)-a)**2)+(2-np.sqrt(2))/np.pi*np.linalg.norm(Pn@B)**2
        rows=[]
        for e in [.02,.01,.005,.002]:
          def path(t):
            R=np.sqrt(1-t*t*np.sum(B*B,axis=1));th=theta+t*a
            return np.column_stack([R*np.cos(th),R*np.sin(th),t*B])
          measured=(objective(path(e))+objective(path(-e))-2*f0)/(2*e*e)
          rows.append({'epsilon':e,'coefficient':float(measured),'absolute_error':float(abs(measured-expected))})
        finite([r['coefficient'] for r in rows],'second variations')
        require(rows[-1]['absolute_error']<2e-5,'Second variation mismatch')
        hess.append({'ambient_dimension':d,'path':k,'predicted_coefficient':float(expected),'measurements':rows})
    report['second_variation_numerical_checks']=hess
    # Independent direct integration sanity check. Scrambled Sobol estimates,
    # explicitly not deterministic quadrature error certificates.
    z=norm.ppf(qmc.Sobol(4,scramble=True,seed=89340).random_base2(19))
    candidates={'cosine':Vn,'independent':np.eye(4),'tetrahedron':np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3),
      'triple_plus_one':np.array([[1,0,0],[-.5,np.sqrt(3)/2,0],[-.5,-np.sqrt(3)/2,0],[0,0,1]])}
    integration=[]
    for name,U in candidates.items():
      val=objective(U);est=np.min((z[:,:U.shape[1]]@U.T)**2,axis=1).mean()
      finite([val,est,est-val],'integration')
      integration.append({'case':name,'elementary_formula':float(val),'sobol_estimate':float(est),'difference':float(est-val),'samples':len(z)})
    report['direct_integration_sanity_checks']=integration
    searches=[]
    for d in (3,4):
      results=load_batches(d)
      if results:
        vals=[r['value'] for r in results]
        block=float(1-sp.Rational(3,2)/sp.pi*(4+sp.sqrt(3)-sp.sqrt(15)))
        counts={'cosine_value':sum(abs(x-f0)<1e-6 for x in vals),
          'triple_plus_one_value':sum(abs(x-block)<1e-6 for x in vals),
          'independence_value':sum(abs(x-(1-(12-16/np.sqrt(3))/np.pi))<1e-6 for x in vals)}
        counts['other']=len(vals)-sum(counts.values())
        fresh=[]
        for rec in results:
          U=frame(np.array(rec['parameters']),d)
          fresh.append(float(finite(abs(rec['value']-intrinsic_objective(U)),'endpoint discrepancy')))
        require(len(fresh)==100 and max(fresh)<1e-9,'Saved endpoint objective mismatch')
        searches.append({'ambient_dimension':d,'starts':len(results),'minimum_reported':float(min(vals)),
          'optimiser_success_count':sum(r['success'] for r in results),
          'unsuccessful_start_ids':[r['start'] for r in results if not r['success']],
          'termination_detail':'Historical records lack termination message and final gradient; success flags retained, no trajectories rerun',
          'value_class_counts_tolerance_1e_6':counts,
          'maximum_independent_reevaluation_discrepancy':float(max(fresh)),
          'counterexample_below_target_minus_1e_6':any(x<f0-1e-6 for x in vals)})
    report['multistart_search']=searches
    report['not_checked']=['global rank-three inequality','global rank-four inequality','uniform interval bounds for the remaining covariance domain','explicit radius for the local stability neighbourhood','external review or proof-assistant formalisation']
    out=ROOT/'audit_report.json';out.write_text(json.dumps(report,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
    print(json.dumps({'exact_checks':checks,'cross_checks':cross,'searches':searches,'status':report['status']},indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
    print('LOCAL ALGEBRA AND NUMERICAL AUDITS PASSED; GLOBAL INEQUALITY UNRESOLVED')

if __name__=='__main__':main()
