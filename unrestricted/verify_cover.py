import mpmath as mp
import numpy as np
from functools import lru_cache
import time, sys, os, multiprocessing as mproc
mp.mp.dps=35
mp.iv.dps=25
iv=mp.iv
R=iv.mpf('0.8882517'); C=1-R
N=40

def lo(v): return mp.mpf(v._mpi_[0])
def hi(v): return mp.mpf(v._mpi_[1])
def vv(x): return iv.mpf(float(x))
def box(l,u): return iv.mpf([l,u])
def pos(z): return box(max(0,lo(z)),max(0,hi(z)))
def vmin(x,y): return box(min(lo(x),lo(y)),min(hi(x),hi(y)))
def vmax(x,y): return box(max(lo(x),lo(y)),max(hi(x),hi(y)))
def clip01(x): return box(max(0,lo(x)),min(1,hi(x)))

@lru_cache(maxsize=30000)
def coeff(af):
    a=vv(af); k=1-a
    cs=[1/k]
    for n in range(1,N+1):cs.append(cs[-1]*(a+n-1)*(k+n-1)/(n*(k+n)))
    # B(k,k)=2 B_{1/2}(k,k): no gamma-function oracle is used.
    z=iv.mpf('.5'); P=cs[-1]
    for c in reversed(cs[:-1]): P=P*z+c
    tail=cs[-1]*z**(N+1)/(1-z)
    B=2*z**k*(P+box(0,hi(tail)))
    return a,k,B,cs

def Bz(af,z):
    a,k,B,cs=coeff(float(af)); z=clip01(z)
    comp=lo(z)>mp.mpf('.5')
    if comp:z=1-z
    # z <= 1/2, up to tiny interval-rounding uncertainty.
    P=cs[-1]
    for c in reversed(cs[:-1]):P=P*z+c
    tail=cs[-1]*z**(N+1)/(1-z)
    result=z**k*(P+box(0,hi(tail)))
    return B-result if comp else result

def powx(af,p):return vv(p)**(1/vv(af))
def A(af,p):
    a=vv(af);k=1-a; pp=vv(p)
    return (k-pp+a*powx(af,p))/k

def T(af,z):
    k=1-vv(af);z=clip01(z)
    return (1-(1-z)**k)/k

def L(af,p,r):
    a=vv(af);k=1-a;pp=vv(p);rr=vv(r)
    x=powx(af,p);y=powx(af,r)
    if hi(x+y)<1:
        return T(af,y)/rr+T(af,x)/pp+Bz(af,1-x)-Bz(af,y)
    if lo(x+y)>1:
        return ((rr-a*y)/k-A(af,p))/(pp*rr)
    # Extremely rare equality case: use continuous formula with a safe hull.
    v1=T(af,y)/rr+T(af,x)/pp+Bz(af,1-x)-Bz(af,y)
    v2=((rr-a*y)/k-A(af,p))/(pp*rr)
    return box(min(lo(v1),lo(v2)),max(hi(v1),hi(v2)))

def Ip(af,r,smax):
    y=powx(af,r);smax=clip01(smax)
    return T(af,vmin(y,smax))+vv(r)*pos(Bz(af,smax)-Bz(af,y))
def Ir(af,p,ymin):
    smax=1-powx(af,p);ymin=clip01(ymin)
    start=vmax(ymin,smax)
    k=1-vv(af)
    return (1-start**k)/k+vv(p)*pos(Bz(af,smax)-Bz(af,ymin))

def check_val(b):
    al,pl,rl,au,pu,ru=map(float,b)
    v=C*A(au,pu)/(vv(pu)*vv(ru))+1/(1-vv(al))-R*L(au,pl,rl)
    return lo(v)>0,float(lo(v))

def check_grad(b,typ):
    al,pl,rl,au,pu,ru=map(float,b)
    if typ==3:
        q=vv(pl)/(1-vv(al))-R*Ir(au,pu,powx(al,rl))
        return lo(q)>0,float(lo(q))
    if typ==4:
        q=vv(pu)/(1-vv(au))-R*Ir(al,pl,powx(au,ru))
        return hi(q)<0,float(-hi(q))
    if typ==1:
        ls=(1-vv(pl)**((1-vv(al))/vv(al)))/(1-vv(au))
        q=-C*ls+vv(rl)/(1-vv(al))-R*Ip(au,ru,1-powx(al,pl))
        return lo(q)>0,float(lo(q))
    ls=(1-vv(pu)**((1-vv(au))/vv(au)))/(1-vv(al))
    q=-C*ls+vv(ru)/(1-vv(au))-R*Ip(al,rl,1-powx(au,pu))
    return hi(q)<0,float(-hi(q))

def chunk(task):
    typ,bs,classes=task
    bad=[];mn=float('inf')
    for j,b in enumerate(bs):
        try:
            ok,v=check_val(b) if typ=='val' else check_grad(b,int(classes[j]))
        except Exception as e:
            bad.append((j,str(e),b.tolist()));continue
        mn=min(mn,v)
        if not ok:bad.append((j,v,b.tolist()))
    return len(bs),mn,bad

if __name__=='__main__':
    path=sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(__file__),'cover.npz')
    ds=np.load(path)
    lim=int(sys.argv[2]) if len(sys.argv)>2 else 0
    tasks=[]
    for typ in ['val','grad']:
        bs=ds[typ]
        if lim:bs=bs[:lim]
        cls=ds['classes'][:len(bs)] if typ=='grad' else np.zeros(len(bs),dtype=int)
        for i in range(0,len(bs),500):tasks.append((typ,bs[i:i+500],cls[i:i+500]))
    t=time.time();cnt=0;mn=float('inf');bads=[]
    workers=int(sys.argv[3]) if len(sys.argv)>3 else 8
    with mproc.Pool(workers) as pool:
        for n,v,bad in pool.imap_unordered(chunk,tasks):
            cnt+=n;mn=min(mn,v);bads.extend(bad)
            if cnt%10000==0 or bad: print('verified',cnt,'min_margin',mn,'bad',len(bads),'seconds',time.time()-t,flush=True)
    print('FINAL',cnt,'min_margin',mn,'bad',len(bads),'seconds',time.time()-t,flush=True)
    if bads:
        import json
        with open(os.path.join(os.path.dirname(__file__),'failed_cover_checks.json'),'w') as f:json.dump(bads,f)
        sys.exit(1)
