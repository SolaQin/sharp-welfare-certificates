from fractions import Fraction as F
from itertools import combinations
import numpy as np
import mpmath as mp
mp.mp.dps=55
# Exact arithmetic for the shortened 8/9 proof.
p,r,q=F(13,20),F(8,25),F(15749,25000)
cs=[F(1),F(1,3),F(2,9),F(14,81),F(35,243)]
E=lambda z:sum(c*z**(3*n+3)/(n+1) for n,c in enumerate(cs))
C=lambda z:sum(c/(n+F(2,3))*(q/F(2)**n-z**(3*n+2)) for n,c in enumerate(cs))
assert q**3<F(1,4)
bound=p*E(r)+r*E(p)+p*r*(C(r)+C(p))
assert bound>F(1901,5120)+F(3,100000)
print('RATIONAL: five-term lower bound exceeds threshold by',float(bound-F(1901,5120)))
# Independent numerical differentiation of the closed integral formula.
def A(a,p,r):return 1-(p-a*p**(1/a))/(1-a)
def I(a,p,r):
    k=1-a;x=p**(1/a);y=r**(1/a)
    T=lambda z:(1-(1-z)**k)/k
    assert x+y<1
    return p*T(y)+r*T(x)+p*r*mp.betainc(k,k,y,1-x)
def analytic(a,p,r):
    k=1-a;x=p**(1/a);y=r**(1/a);xa=-x*mp.log(x)/a;ya=-y*mp.log(y)/a
    T=lambda z,m:mp.quad(lambda t:(1-t)**(-a)*(-mp.log(1-t))**m,[0,z])
    J=[mp.quad(lambda t:(t*(1-t))**(-a)*(-mp.log(t*(1-t)))**m,[y,mp.mpf('.5'),1-x]) for m in range(3)]
    hx=(1-x)**(-a);hy=(1-y)**(-a)
    da=[-p*T(1-x,1),-(1-x**k)/k,mp.mpf(0)]
    di=[p*T(y,1)+r*T(x,1)+p*r*J[1],T(y,0)+r*J[0],T(x,0)+p*J[0]]
    ha=[-p*T(1-x,2)+x*mp.log(x)**2/a,-T(1-x,1)+xa/p,0,x/(a*p*p),0,0]
    hi=[p*T(y,2)+r*T(x,2)+p*r*J[2]-p*y*mp.log(y)**2/a*hy-r*x*mp.log(x)**2/a*hx,T(y,1)+r*J[1]-r*xa/p*hx,T(x,1)+p*J[1]-p*ya/r*hy,-r*x/(a*p*p)*hx,J[0],-p*y/(a*r*r)*hy]
    return da,di,ha,hi
orders1=[(1,0,0),(0,1,0),(0,0,1)]
orders2=[(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
maxerr=mp.mpf(0)
for point in [('.3231442122','.6544790312','.3251658297'),('.2','.4','.5'),('.6','.3','.2'),('.4','.7','.2')]:
    vals=tuple(map(mp.mpf,point));da,di,ha,hi=analytic(*vals)
    for f,orders,expected in [(A,orders1,da),(I,orders1,di),(A,orders2,ha),(I,orders2,hi)]:
        for order,target in zip(orders,expected):
            error=abs(mp.diff(f,vals,order)-target);maxerr=max(maxerr,error)
            assert error<mp.mpf('1e-45'),(point,order,error)
print('DERIVATIVES: all first and second partials at four points agree; max error',mp.nstr(maxerr,5))
# Verify the compact boundary-payment formula and BIC inequalities.
for a,p,r in [(.3231442122,.6544790312,.3251658297),(1/3,.5,.5),(.2,.2,.7)]:
    k=1-a;x=p**(1/a);y=r**(1/a);smax=1-x;theta=1/(2-p**(k/a)-r**(k/a));assert theta<=1
    B=np.unique(np.r_[np.linspace(0,1,151),y,1]);S=np.unique(np.r_[np.linspace(0,1,151),smax,0])
    X=np.where(B<y,0,np.where(B<1,p*theta,p+(1-p)*theta))
    Y=np.where(S>smax,0,np.where(S>0,r*theta,r+(1-r)*theta))
    tb=np.where(B<y,0,np.where(B<1,p*theta*y,p*theta*y+p+(1-2*p)*theta))
    H=r*theta*smax;ts=np.where(S>smax,0,H)
    ub=B[:,None]*X[None,:]-tb[None,:];us=ts[None,:]-S[:,None]*Y[None,:]
    assert np.max(ub-np.diag(ub)[:,None])<1e-12
    assert np.max(us-np.diag(us)[:,None])<1e-12
    assert min(np.min(np.diag(ub)),np.min(np.diag(us)))>-1e-12
    assert abs((1-r)*p*theta*y+r*(p*theta*y+p+(1-2*p)*theta)-H)<1e-12
print('PAYMENTS: boundary lottery, zero expected surplus, BIC and interim IR checks pass')
# Compare ironing with direct vertex enumeration of small finite allocation LPs.
rng=np.random.default_rng(17092026)
def pooled(values,weights,ascending):
    out=[]
    for i,(v,w) in enumerate(zip(values,weights)):
        out.append([v*w,w,[i]])
        while len(out)>1:
            left,right=out[-2:];lv,rv=left[0]/left[1],right[0]/right[1]
            if (lv<=rv if ascending else lv>=rv):break
            out[-2:]=[[left[0]+right[0],left[1]+right[1],left[2]+right[2]]]
    result=np.empty(len(values))
    for v,w,inds in out:result[inds]=v/w
    return result
maxdiff=0.
for t in range(12):
    n,m=(3,2) if t%2 else (2,3);dim=n*m
    b=np.sort(rng.random(n))[::-1];s=np.sort(rng.random(m));w=rng.dirichlet(np.ones(n));z=rng.dirichlet(np.ones(m))
    q=np.cumsum(w)-w;Q=np.cumsum(z)-z
    phi=b-(np.r_[b[0],b[:-1]]-b)*q/w
    psi=s+(s-np.r_[s[0],s[:-1]])*Q/z
    rows=list(np.eye(dim))+list(-np.eye(dim));rhs=[1.]*dim+[0.]*dim
    for i in range(n-1):
        row=np.zeros((m,n));row[:,i]=-z;row[:,i+1]=z;rows.append(row.ravel());rhs.append(0.)
    for j in range(m-1):
        row=np.zeros((m,n));row[j,:]=-w;row[j+1,:]=w;rows.append(row.ravel());rhs.append(0.)
    mat=np.array(rows);rhs=np.array(rhs);vertices=[]
    for active in combinations(range(len(rhs)),dim):
        try:v=np.linalg.solve(mat[list(active)],rhs[list(active)])
        except np.linalg.LinAlgError:continue
        if np.max(mat@v-rhs)<1e-9:vertices.append(v)
    vertices=np.array(vertices)
    for lam in [.1,1.,10.]:
        a=lam/(1+lam);u=np.maximum(pooled((1-a)*b+a*phi,w,False),0);v=pooled((1-a)*s+a*psi,z,True)
        direct=np.max(vertices@((b[None,:]-s[:,None]+lam*(phi[None,:]-psi[:,None]))*z[:,None]*w[None,:]).ravel())
        iron=np.sum(np.maximum(u[None,:]-v[:,None],0)*z[:,None]*w[None,:])/(1-a)
        maxdiff=max(maxdiff,abs(direct-iron));assert abs(direct-iron)<1e-8
print('FINITE DUAL: 36 LP/ironing comparisons pass; max difference',maxdiff)
print('ALL INDEPENDENT MATHEMATICAL SANITY CHECKS PASS (not substitutes for proofs)')
