"""Independent floating-point evaluator through second intrinsic volumes.
Full-rank and rank-three frames; NOT an interval certificate.
"""
import numpy as np
from itertools import combinations,product
from math import pi
from scipy.spatial import ConvexHull
from validation import frame_input,finite

def triangle_twice_area(V):
    a,b=V[1]-V[0],V[2]-V[0]
    aa=float(a@a);bb=float(b@b);ab=float(a@b)
    return np.sqrt(max(0.,aa*bb-ab*ab))

def intrinsic_objective(U):
    U=frame_input(U,exact_four=True)
    _,s,vh=np.linalg.svd(U,full_matrices=False)
    d=int(sum(s>1e-10));U=U@vh[:d].T
    if d==1:return 1.
    if d==2:
        theta=np.sort(np.arctan2(U[:,1],U[:,0])%pi)
        return 1.-sum(np.sin(np.diff(np.r_[theta,theta[0]+pi])))/pi
    C=U@U.T
    base=2*sum(np.sqrt(max(0.,1-C[i,j]**2)) for i,j in combinations(range(4),2))
    total=0.
    hull_area=ConvexHull(np.vstack([U,-U])).area if d==3 else 0.
    for T in combinations(range(4),3):
        k=next(iter(set(range(4))-set(T)))
        for s,t in product((-1,1),repeat=2):
            signs=np.array([1,s,t]);V=U[list(T)]*signs[:,None]
            r=triangle_twice_area(V)
            if d==3:
                wt=1.
            else:
                a=np.zeros(4);a[list(T)]=signs;a[k]=1.
                b=a.copy();b[k]=-1.
                x=np.linalg.solve(U,a);y=np.linalg.solve(U,b)
                x=x/np.linalg.norm(x);y=y/np.linalg.norm(y)
                ang=np.arctan2(np.linalg.norm(x-y)*np.linalg.norm(x+y),2*float(x@y))
                wt=1.-ang/pi
            total+=r*wt
    return float(finite(1.-(base-.5*total+.5*hull_area)/pi,'intrinsic objective'))

if __name__=='__main__':
    from objective import objective
    rng=np.random.default_rng(9926423)
    errors={}
    for d in (3,4):
        errs=[]
        for k in range(100):
            U=rng.normal(size=(4,d))
            x=objective(U);y=intrinsic_objective(U)
            errs.append(abs(x-y))
        errors[d]=max(errs)
    print(errors)
