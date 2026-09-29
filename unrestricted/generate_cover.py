"""Propose a finite binary cover; all analytic tests are rechecked by verify_cover.py."""
import numpy as np
from scipy.special import beta, betainc
import os
R=.8882517
CENTER=np.array([.323144212181945529204914052954389108756940743,
                 .654479031166666393676935731392503971225432437,
                 .325165829719811989207067293084997587516390811])
RAD=.012
ROOT=np.array([.1,.02,.00005,.999999,.9,1.])

def A(a,p):
    k=1-a;x=np.exp(np.log(p)/a)
    return (k-p+a*x)/k

def G(a,p,r):
    k=1-a;x=np.exp(np.log(p)/a);y=np.exp(np.log(r)/a)
    J=beta(k,k)*(betainc(k,k,1-x)-betainc(k,k,y))
    v=(p*(-np.expm1(k*np.log1p(-y)))+r*(-np.expm1(k*np.log1p(-x))))/k+p*r*J
    return np.where(x+y>=1,(r-a*y)/k-A(a,p),v)

def kernel(a,l,u):
    k=1-a
    return np.maximum(0,beta(k,k)*(betainc(k,k,np.clip(u,0,1))-betainc(k,k,np.clip(l,0,1))))

def Ip(a,r,s):
    k=1-a;y=np.exp(np.log(r)/a);end=np.minimum(y,s)
    return -np.expm1(k*np.log1p(-end))/k+r*kernel(a,y,s)

def Ir(a,p,y):
    k=1-a;s=-np.expm1(np.log(p)/a);start=np.maximum(y,s)
    return -np.expm1(k*np.log(start))/k+p*kernel(a,y,s)

def gradients(bs):
    al,pl,rl=bs[:,:3].T;au,pu,ru=bs[:,3:].T
    xmin=np.exp(np.log(pl)/al);xmax=np.exp(np.log(pu)/au)
    ymin=np.exp(np.log(rl)/al);ymax=np.exp(np.log(ru)/au)
    Lmin=-np.expm1(((1-au)/au)*np.log(pu))/(1-al)
    Lmax=-np.expm1(((1-al)/al)*np.log(pl))/(1-au)
    fpl=-(1-R)*Lmax+rl/(1-al)-R*Ip(au,ru,1-xmin)
    fpu=-(1-R)*Lmin+ru/(1-au)-R*Ip(al,rl,1-xmax)
    frl=pl/(1-al)-R*Ir(au,pu,ymin)
    fru=pu/(1-au)-R*Ir(al,pl,ymax)
    return fpl,fpu,frl,fru

def generate():
    bs=ROOT.reshape(1,6);ids=np.array([1],dtype=np.uint64)
    leaves={t:[]for t in ['val','grad','local']};leafids={t:[]for t in leaves};classes=[]
    tested=0
    for depth in range(63):
        loc=np.all(bs[:,:3]>=CENTER-RAD,axis=1)&np.all(bs[:,3:]<=CENTER+RAD,axis=1)
        leaves['local'].append(bs[loc]);leafids['local'].append(ids[loc])
        bs=bs[~loc];ids=ids[~loc]
        if not len(bs):break
        tested+=len(bs)
        au,pu,ru=bs[:,3:].T;al,pl,rl=bs[:,:3].T
        lower=(1-R)*A(au,pu)/(pu*ru)+1/(1-al)-R*G(au,pl,rl)/(pl*rl)
        good=lower>1e-8
        leaves['val'].append(bs[good]);leafids['val'].append(ids[good])
        bs=bs[~good];ids=ids[~good]
        fpl,fpu,frl,fru=gradients(bs)
        good=(fpl>1e-7)|(fpu< -1e-7)|(frl>1e-7)|(fru< -1e-7)
        typ=np.where(frl>1e-7,3,np.where(fru< -1e-7,4,np.where(fpl>1e-7,1,2)))
        leaves['grad'].append(bs[good]);leafids['grad'].append(ids[good]);classes.append(typ[good])
        bs=bs[~good];ids=ids[~good]
        if not len(bs):break
        w=bs[:,3:]-bs[:,:3];dim=np.argmax(w,axis=1);ii=np.arange(len(bs))
        mid=(bs[ii,dim]+bs[ii,dim+3])/2
        left=bs.copy();right=bs.copy();left[ii,dim+3]=mid;right[ii,dim]=mid
        bs=np.concatenate([left,right]);ids=np.concatenate([ids*2,ids*2+1])
    else:raise RuntimeError('Cover did not terminate.')
    out={t:np.concatenate(leaves[t])for t in leaves}
    out.update({t+'_ids':np.concatenate(leafids[t])for t in leaves})
    out.update(classes=np.concatenate(classes),center=CENTER,rad=np.array(RAD),root=ROOT)
    path=os.path.join(os.path.dirname(__file__),'cover.npz')
    np.savez_compressed(path,**out)
    print('Cover generated:',{t:len(out[t])for t in leaves},'tested nodes',tested,'last depth',depth)

if __name__=='__main__':
    with np.errstate(all='ignore'):generate()
