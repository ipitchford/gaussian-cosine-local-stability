"""Standard-library exact check of the local quadratic-form identities.
This does NOT check the analytic differentiation or the global conjecture.
Numbers are pairs of rational coefficients for a+b*sqrt(2).
"""
from fractions import Fraction as F
import json
from pathlib import Path

def q(a=0,b=0):return (F(a),F(b))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[0],-x[1])
def sub(x,y):return add(x,neg(y))
def mul(x,y):return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def zeros(n,m):return [[q() for _ in range(m)] for _ in range(n)]
def eye(n):return [[q(int(i==j)) for j in range(n)] for i in range(n)]
def tr(A):return [list(x) for x in zip(*A)]
def scale(A,x):return [[mul(v,x) for v in row] for row in A]
def plus(A,B):return [[add(x,y) for x,y in zip(a,b)] for a,b in zip(A,B)]
def minus(A,B):return plus(A,scale(B,q(-1)))
def mm(A,B):
    C=zeros(len(A),len(B[0]))
    for i in range(len(A)):
      for j in range(len(B[0])):
        for k in range(len(B)):C[i][j]=add(C[i][j],mul(A[i][k],B[k][j]))
    return C

def checks():
    r=q(0,1);half=q(F(1,2));I=eye(4)
    V=[[q(1),q()],[q(0,F(1,2)),q(0,F(1,2))],[q(),q(1)],[q(0,F(-1,2)),q(0,F(1,2))]]
    ints=[[0,1,0,-1],[1,0,1,0],[0,1,0,1],[-1,0,1,0]]
    A=[[q(x) for x in row] for row in ints]
    P=minus(I,scale(mm(V,tr(V)),half))
    lhs=scale(minus(scale(I,q(0,F(1,2))),scale(A,half)),q(-1,1))
    rhs=scale(P,q(2,-1))
    results={'sqrt2_squared':mul(r,r)==q(2),
      'frame_columns':mm(tr(V),V)==scale(eye(2),q(2)),
      'projection_symmetric':tr(P)==P,
      'projection_idempotent':mm(P,P)==P,
      'signed_cycle_square':mm(A,A)==scale(I,q(2)),
      'rotation_kernel':mm(P,V)==zeros(4,2),
      'transverse_quadratic_form':lhs==rhs}
    if not all(results.values()):raise ValueError(results)
    bad=[row[:] for row in A];bad[0][1]=q(2)
    if mm(bad,bad)==scale(I,q(2)):raise ValueError('Deliberate mutation escaped check')
    results['deliberate_mutation_rejected']=True
    return results

if __name__=='__main__':
    report={'scope':'Local algebra only; NOT a proof of the global conjecture','checks':checks()}
    Path(__file__).with_name('local_algebra_report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
