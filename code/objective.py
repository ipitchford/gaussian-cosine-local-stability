"""Deterministic elementary evaluation; floating point, NOT a proof certificate."""
from __future__ import annotations
import numpy as np
from numpy.linalg import norm
from scipy.linalg import null_space
from math import pi, sqrt, atan, hypot
from validation import frame_input,finite,NumericalDomainError

def interval_mass(lo: float, hi: float)->float:
    def h(x): return np.sign(x) if np.isinf(x) else x/hypot(1.,x)
    return .5*(h(hi)-h(lo))

def rank3_weight(a:float,b:np.ndarray,c:np.ndarray)->float:
    intervals=[]
    for bk,ck in zip(b,c):
        if abs(ck)<1e-13:
            if abs(bk)<a-1e-13: return 0.
            continue
        q=sorted([(-a-bk)/ck,(a-bk)/ck])
        intervals.append(q)
    intervals.sort()
    merged=[]
    for lo,hi in intervals:
        if merged and lo<=merged[-1][1]: merged[-1][1]=max(hi,merged[-1][1])
        else: merged.append([lo,hi])
    return max(0.,1.-sum(interval_mass(lo,hi) for lo,hi in merged))

def polygon_t2_probability(vertices:np.ndarray)->float:
    """Integral (pi*(1+|x|^2)^2)^-1 over a convex CCW polygon."""
    total=0.
    for p,q in zip(vertices,np.roll(vertices,-1,axis=0)):
        v=q-p; length=norm(v)
        if length==0:continue
        v=v/length
        h=p[0]*v[1]-p[1]*v[0] # outward normal (v_y,-v_x), dot p
        t0=float(p@v); aa=hypot(1.,h)
        total += h/aa*(atan((t0+length)/aa)-atan(t0/aa))
    return total/(2*pi)

def rank4_weight(a:float,b:np.ndarray,B:np.ndarray)->float:
    lens=norm(B,axis=1)
    if min(lens)<1e-13:
        # Singular near-degeneracy should be handled by reducing ambient rank.
        raise NumericalDomainError('Near-singular rank-four residual lengths')
    m=[]
    for bk,ck in zip(b,lens):
        m.append(interval_mass((-a-bk)/ck,(a-bk)/ck))
    rhs=np.array([[-a,-a],[a,-a],[a,a],[-a,a]])-b
    if abs(np.linalg.det(B))<1e-14: raise NumericalDomainError('Near-singular rank-four residual determinant')
    corners=np.linalg.solve(B,rhs.T).T
    if np.linalg.det(B)<0: corners=corners[::-1]
    both=polygon_t2_probability(corners)
    return max(0.,1.-sum(m)+both)

def objective(U:np.ndarray, reduce=True)->float:
    U=frame_input(U)
    if reduce:
        _,s,vh=np.linalg.svd(U,full_matrices=False)
        rank=int(np.sum(s>1e-10))
        U=U@vh[:rank].T
    d=U.shape[1]
    if d==1:return 1.
    if d==2:
        angles=np.sort(np.mod(np.arctan2(U[:,1],U[:,0]),pi))
        gaps=np.diff(np.r_[angles,angles[0]+pi])
        return 1.-np.sum(np.sin(gaps))/pi
    n=len(U); total=0.
    for i in range(n):
        for j in range(i):
            cij=float(U[i]@U[j]); fac=sqrt(max(0.,1.-cij*cij))
            if fac<1e-12: raise NumericalDomainError('Duplicate or nearly duplicate lines need preprocessing')
            others=[k for k in range(n) if k not in (i,j)]
            for sign in (-1,1):
                normal=(U[i]-sign*U[j]); normal/=norm(normal)
                w=(U[i]+sign*U[j]); w/=norm(w)
                a=sqrt((1.+sign*cij)/2.)
                Z=null_space(np.stack([normal,w]))
                b=U[others]@w
                B=U[others]@Z
                if d==3: wt=rank3_weight(a,b,B[:,0])
                elif d==4: wt=rank4_weight(a,b,B)
                else:raise ValueError('d<=4 required')
                total+=fac*wt
    return float(finite(1.-total/pi,'elementary objective'))

if __name__=='__main__':
    from scipy.stats import qmc, norm as snorm
    seed=89340
    rng=np.random.default_rng(seed)
    z=snorm.ppf(qmc.Sobol(4,scramble=True,seed=seed).random_base2(19))
    tests={'identity':np.eye(4),'tetra':np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]]),'cosine':np.stack([np.cos(np.arange(4)*pi/4),np.sin(np.arange(4)*pi/4)],axis=1)}
    for d in (3,4):
        for k in range(3):tests[f'random{d}_{k}']=rng.normal(size=(4,d))
    for label,U in tests.items():
        U=U/norm(U,axis=1)[:,None]
        exact=objective(U)
        mc=np.min((z[:,:U.shape[1]]@U.T)**2,axis=1).mean()
        print(label,exact,mc,mc-exact)
