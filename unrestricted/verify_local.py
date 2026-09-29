import mpmath as mp
import numpy as np
import math, time, os, multiprocessing as mproc
from fractions import Fraction
from constants import CENTER_STR, RADIUS_STR
mp.mp.dps=35;mp.iv.dps=25
iv=mp.iv

def lo(v):return mp.mpf(v._mpi_[0])
def hi(v):return mp.mpf(v._mpi_[1])
def box(l,u):return iv.mpf([l,u])

def Tpoint(a,z,m):
    a=iv.mpf(a);z=iv.mpf(z);k=1-a
    w=-k*iv.ln(1-z)
    pol=iv.mpf(1)
    for j in range(1,m+1):pol+=w**j/math.factorial(j)
    return math.factorial(m)*(1-iv.exp(-w)*pol)/k**(m+1)

def Tbound(a,z,m):
    l=Tpoint(lo(a),lo(z),m);u=Tpoint(hi(a),hi(z),m)
    return box(lo(l),hi(u))

def Jpoint(a,l,u,m,N=256):
    a=iv.mpf(a);l=iv.mpf(l);u=iv.mpf(u);h=(u-l)/N
    def f(t):
        L=-iv.ln(t*(1-t))
        return iv.exp(a*L)*L**m
    mid=iv.mpf(0);trap=(f(l)+f(u))/2
    for j in range(N):
        mid+=f(l+(iv.mpf(j)+iv.mpf('.5'))*h)
        if j:trap+=f(l+j*h)
    # f is convex: midpoint sum is a lower bound, trapezoid sum an upper bound.
    return box(lo(mid*h),hi(trap*h))

def Jbound(a,y,u,m):
    l=Jpoint(lo(a),hi(y),lo(u),m)
    h=Jpoint(hi(a),lo(y),hi(u),m)
    return box(lo(l),hi(h))

def hessian(center,rad):
    a,p,r=[box(float(t)-rad,float(t)+rad) for t in center]
    R=box(mp.mpf('0.8882516'),mp.mpf('0.8882517'));C=1-R;k=1-a
    x=p**(1/a);y=r**(1/a)
    assert hi(x+y)<1
    lx=iv.ln(x);ly=iv.ln(y)
    xa=-x*lx/a;ya=-y*ly/a
    J=[Jbound(a,y,1-x,m) for m in range(3)]
    t1x,t1y,t1s=Tbound(a,x,1),Tbound(a,y,1),Tbound(a,1-x,1)
    t2x,t2y,t2s=Tbound(a,x,2),Tbound(a,y,2),Tbound(a,1-x,2)
    sx=(1-x)**(-a);sy=(1-y)**(-a)
    aaa=-p*t2s+x*lx**2/a
    aap=-t1s+xa/p
    gaa=p*t2y+p*r*J[2]+r*t2x-p*y*ly**2/a*sy-r*x*lx**2/a*sx
    gap=t1y+r*J[1]-r*xa/p*sx
    gar=t1x+p*J[1]-p*ya/r*sy
    H=[[iv.mpf(0) for j in range(3)]for i in range(3)]
    H[0][0]=C*aaa+2*p*r/k**3-R*gaa
    H[0][1]=H[1][0]=C*aap+r/k**2-R*gap
    H[0][2]=H[2][0]=p/k**2-R*gar
    H[1][1]=x/(a*p**2)*(C+R*r*sx)
    H[2][2]=R*p*y/(a*r**2)*sy
    H[1][2]=H[2][1]=1/k-R*J[0]
    return H

def check(center):
    H=hessian(center,.00150000001)
    for i in range(3):H[i][i]-=iv.mpf('.01')
    d1=H[0][0]
    if lo(d1)<=0:return False,[float(lo(d1))],center.tolist()
    d2=H[1][1]-H[0][1]**2/d1
    if lo(d2)<=0:return False,[float(lo(d1)),float(lo(d2))],center.tolist()
    d3=H[2][2]-H[0][2]**2/d1-(H[1][2]-H[0][1]*H[0][2]/d1)**2/d2
    return lo(d3)>0,[float(lo(z))for z in [d1,d2,d3]],center.tolist()

if __name__=='__main__':
    data=np.load(os.path.join(os.path.dirname(__file__),'cover.npz'))
    center=data['center'];offset=np.arange(-7,8,2)*.0015
    centers=[center+np.array([x,y,z]) for x in offset for y in offset for z in offset]
    # Verify the slightly enlarged subboxes overlap and cover the large box.
    rad=.00150000001
    for j in range(3):
        points=sorted([float(center[j]+x)for x in offset])
        # Compare the actual binary endpoints to an exact rational cube.
        lows=[Fraction(t-rad) for t in points]
        highs=[Fraction(t+rad) for t in points]
        c=Fraction(CENTER_STR[j]); h=Fraction(RADIUS_STR)
        assert lows[0]<=c-h and highs[-1]>=c+h
        assert all(highs[i]>=lows[i+1] for i in range(7))
    t=time.time();mn=[float('inf')]*3;bad=[];cnt=0
    with mproc.Pool(min(12,os.cpu_count() or 1))as pool:
        for ok,ds,c in pool.imap_unordered(check,centers):
            cnt+=1
            if not ok:bad.append((ds,c))
            for j,z in enumerate(ds):mn[j]=min(mn[j],z)
    print('LOCAL',cnt,'min_LDL_pivots',mn,'bad',len(bad),'seconds',time.time()-t)
    if bad:
        print(bad[:10]);raise SystemExit(1)
