import mpmath as mp
mp.mp.dps=85;mp.iv.dps=75
iv=mp.iv

def lo(v):return mp.mpf(v._mpi_[0])
def hi(v):return mp.mpf(v._mpi_[1])
def box(l,u):return iv.mpf([l,u])
a=iv.mpf('0.323144212181945529204914052954389108756940743')
p=iv.mpf('0.654479031166666393676935731392503971225432437')
r=iv.mpf('0.325165829719811989207067293084997587516390811')
k=1-a;x=p**(1/a);y=r**(1/a)
Rlo=iv.mpf('0.88825169032566246988432630129334211037')
Rhi=iv.mpf('0.88825169032566246988432630129334211038')
Rcover=iv.mpf('0.8882517')

def primitive(z,N=700):
    # Return B_z(1-a,1-a) and its a-derivative at a fixed endpoint z.
    c=1/k; H=iv.mpf(0);zn=iv.mpf(1)
    L=H+1/k-iv.ln(z)
    s0=c;s1=c*L
    for n in range(1,N+1):
        H+=1/(a+n-1)
        c=c*(a+n-1)*(k+n-1)/(n*(k+n))
        zn*=z
        L=H+1/(n+k)-iv.ln(z)
        s0+=c*zn;s1+=c*zn*L
    tail0=c*zn*z/(1-z)
    tail1=c*zn*(L*z/(1-z)+z/((a+N)*(1-z)**2))
    return z**k*(s0+box(0,hi(tail0))),z**k*(s1+box(0,hi(tail1)))

def T(z,m):
    w=-k*iv.ln(1-z)
    if m==0:return (1-iv.exp(-w))/k
    return (1-iv.exp(-w)*(1+w))/k**2
j0u,j1u=primitive(1-x);j0l,j1l=primitive(y)
J=j0u-j0l;J1=j1u-j1l
A=(k-p+a*x)/k
G=p*T(y,0)+r*T(x,0)+p*r*J
Aa=-p*T(1-x,1);Ap=-(1-x**k)/k
Ga=p*T(y,1)+r*T(x,1)+p*r*J1
Gp=T(y,0)+r*J;Gr=T(x,0)+p*J

def FG(R):
    F=(1-R)*A+p*r/k-R*G
    grad=[(1-R)*Aa+p*r/k**2-R*Ga,
          (1-R)*Ap+r/k-R*Gp,
          p/k-R*Gr]
    return F,grad
FL,gL=FG(Rlo);FU,gU=FG(Rhi);FC,gC=FG(Rcover)
# Local strong convexity m=0.01 was established in verify_local.py.
m=iv.mpf('.01')
correction=sum(g*g for g in gL)/(2*m)
lower=FL-correction
assert lo(lower)>0
assert hi(FU)<0
# The boundary of the large cube is >=0.0119 from this rational center,
# while coordinatewise displacement sums are <=0.037.
ginf=max(max(abs(lo(g)),abs(hi(g)))for g in gC)
boundary=FC-iv.mpf('.037')*ginf+m*iv.mpf('.0119')**2/2
assert lo(boundary)>0
# Boundary-tie allocation is feasible throughout the extremizer's box.
Q=iv.mpf('1.983') # deliberately lower than (1-a)/a in the large box
assert hi(iv.mpf('.667')**Q+iv.mpf('.338')**Q)<1
assert hi(p**(k/a)+r**(k/a))<1
print('rho_lower =',Rlo)
print('rho_upper =',Rhi)
print('lower_certificate =',lower)
print('upper_certificate =',FU)
print('large_box_boundary_margin =',boundary)
print('candidate_theta =',1/(2-p**(k/a)-r**(k/a)))
print('point_and_boundary_checks: PASS')

# Locate the unique minimizer, using strong convexity and the certified rho interval.
maxg=max(max(abs(lo(g)),abs(hi(g))) for g in gL+gU)
location_bound=iv.sqrt(3)*iv.mpf(maxg)/m
assert hi(location_bound)<mp.mpf('1e-30')
print('distance_of_minimizer_from_rational_center <',hi(location_bound))
