#define main verify_cover_embedded_main
#include "verify_cover_mpfr.cpp"
#undef main

extern "C" {
int mpfr_exp(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sqrt(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
}

namespace {
I dec_interval(const char* l,const char* u){
  I z;
  if(mpfr_set_str(z.l.x,l,10,MPFR_RNDD)||mpfr_set_str(z.u.x,u,10,MPFR_RNDU))
    throw std::runtime_error("bad decimal interval");
  return z;
}
I interval_d(double l,double u){I z;mpfr_set_d(z.l.x,l,MPFR_RNDD);mpfr_set_d(z.u.x,u,MPFR_RNDU);return z;}
I expi(const I&a){I z;mpfr_exp(z.l.x,a.l.x,MPFR_RNDD);mpfr_exp(z.u.x,a.u.x,MPFR_RNDU);return z;}
I sqrti(const I&a){if(mpfr_cmp_ui(a.l.x,0)<0)throw std::runtime_error("sqrt domain");I z;mpfr_sqrt(z.l.x,a.l.x,MPFR_RNDD);mpfr_sqrt(z.u.x,a.u.x,MPFR_RNDU);return z;}
I abs_interval(const I&a){
  if(mpfr_cmp_ui(a.l.x,0)>=0)return a;
  if(mpfr_cmp_ui(a.u.x,0)<=0)return neg(a);
  I z;mpfr_set_ui(z.l.x,0,MPFR_RNDN);Real ll,uu;mpfr_mul(ll.x,a.l.x,a.l.x,MPFR_RNDU);mpfr_mul(uu.x,a.u.x,a.u.x,MPFR_RNDU);z.u=lt(ll,uu)?uu:ll;return z;
}
I scale_ui(const I&a,unsigned long n){return mul(a,I(n));}
I square_interval(const I&a){
  I z;
  if(mpfr_cmp_ui(a.l.x,0)>=0){mpfr_mul(z.l.x,a.l.x,a.l.x,MPFR_RNDD);mpfr_mul(z.u.x,a.u.x,a.u.x,MPFR_RNDU);return z;}
  if(mpfr_cmp_ui(a.u.x,0)<=0){mpfr_mul(z.l.x,a.u.x,a.u.x,MPFR_RNDD);mpfr_mul(z.u.x,a.l.x,a.l.x,MPFR_RNDU);return z;}
  mpfr_set_ui(z.l.x,0,MPFR_RNDN);Real ll,uu;mpfr_mul(ll.x,a.l.x,a.l.x,MPFR_RNDU);mpfr_mul(uu.x,a.u.x,a.u.x,MPFR_RNDU);z.u=lt(ll,uu)?uu:ll;return z;
}
I divide_ui(const I&a,unsigned long n){return divi(a,I(n));}
I pow_general(const I&base,const I&expo){return expi(mul(log_pos(base),expo));}
I factorial(unsigned m){unsigned long v=1;for(unsigned j=2;j<=m;j++)v*=j;return I(v);}

I Tpoint(const I&a,const I&z,unsigned m){
  I k=sub(one(),a);
  I w=neg(mul(k,log_pos(sub(one(),z))));
  I pol=one(),wp=one();unsigned long fact=1;
  for(unsigned j=1;j<=m;j++){wp=mul(wp,w);fact*=j;pol=add(pol,divide_ui(wp,fact));}
  I num=sub(one(),mul(expi(neg(w)),pol));
  return divi(mul(factorial(m),num),pow_ui_pos(k,m+1));
}
I Tbound(const I&a,const I&z,unsigned m){
  I al=from_bounds(a.l,a.l),au=from_bounds(a.u,a.u),zl=from_bounds(z.l,z.l),zu=from_bounds(z.u,z.u);
  I loz=Tpoint(al,zl,m),hiz=Tpoint(au,zu,m);
  return from_bounds(loz.l,hiz.u);
}
I J_integrand(const I&a,const I&t,unsigned m){
  I L=neg(log_pos(mul(t,sub(one(),t))));
  return mul(expi(mul(a,L)),pow_ui_pos(L,m));
}
I Jpoint(const I&a,const I&l,const I&u,unsigned m,unsigned N=256){
  I h=divide_ui(sub(u,l),N);
  I mid=zero();
  I trap=divide_ui(add(J_integrand(a,l,m),J_integrand(a,u,m)),2);
  for(unsigned j=0;j<N;j++){
    I midcoef((unsigned long)(2*j+1));
    I tm=add(l,divide_ui(mul(midcoef,h),2));
    mid=add(mid,J_integrand(a,tm,m));
    if(j){I tj=add(l,scale_ui(h,j));trap=add(trap,J_integrand(a,tj,m));}
  }
  I lower=mul(mid,h),upper=mul(trap,h);
  return from_bounds(lower.l,upper.u);
}
I Jbound(const I&a,const I&y,const I&u,unsigned m){
  I al=from_bounds(a.l,a.l),au=from_bounds(a.u,a.u);
  I yl=from_bounds(y.l,y.l),yu=from_bounds(y.u,y.u);
  I ul=from_bounds(u.l,u.l),uu=from_bounds(u.u,u.u);
  I low=Jpoint(al,yu,ul,m),high=Jpoint(au,yl,uu,m);
  return from_bounds(low.l,high.u);
}

std::array<std::array<I,3>,3> hessian(const std::array<double,3>&center,double rad){
  I a=interval_d(center[0]-rad,center[0]+rad);
  I p=interval_d(center[1]-rad,center[1]+rad);
  I r=interval_d(center[2]-rad,center[2]+rad);
  I R=dec_interval("0.8882516","0.8882517"),C=sub(one(),R),k=sub(one(),a);
  I x=pow01(p,reciprocal_pos(a)),y=pow01(r,reciprocal_pos(a));
  if(mpfr_cmp_ui(add(x,y).u.x,1)>=0)throw std::runtime_error("local nonoverlap failed");
  I lx=log_pos(x),ly=log_pos(y);
  I xa=divi(neg(mul(x,lx)),a),ya=divi(neg(mul(y,ly)),a);
  std::array<I,3> J={Jbound(a,y,sub(one(),x),0),Jbound(a,y,sub(one(),x),1),Jbound(a,y,sub(one(),x),2)};
  I t1x=Tbound(a,x,1),t1y=Tbound(a,y,1),t1s=Tbound(a,sub(one(),x),1);
  I t2x=Tbound(a,x,2),t2y=Tbound(a,y,2),t2s=Tbound(a,sub(one(),x),2);
  I sx=pow_general(sub(one(),x),neg(a)),sy=pow_general(sub(one(),y),neg(a));
  I aaa=add(neg(mul(p,t2s)),divi(mul(x,square_interval(lx)),a));
  I aap=add(neg(t1s),divi(xa,p));
  I gaa=add(add(add(mul(p,t2y),mul(mul(p,r),J[2])),mul(r,t2x)),
            add(neg(divi(mul(mul(mul(p,y),square_interval(ly)),sy),a)),
                neg(divi(mul(mul(mul(r,x),square_interval(lx)),sx),a))));
  I gap=add(add(t1y,mul(r,J[1])),neg(mul(divi(mul(r,xa),p),sx)));
  I gar=add(add(t1x,mul(p,J[1])),neg(mul(divi(mul(p,ya),r),sy)));
  std::array<std::array<I,3>,3> H;
  H[0][0]=add(add(mul(C,aaa),divi(scale_ui(mul(p,r),2),pow_ui_pos(k,3))),neg(mul(R,gaa)));
  H[0][1]=H[1][0]=add(add(mul(C,aap),divi(r,pow_ui_pos(k,2))),neg(mul(R,gap)));
  H[0][2]=H[2][0]=add(divi(p,pow_ui_pos(k,2)),neg(mul(R,gar)));
  H[1][1]=mul(divi(x,mul(a,pow_ui_pos(p,2))),add(C,mul(mul(R,r),sx)));
  H[2][2]=mul(divi(mul(mul(R,p),y),mul(a,pow_ui_pos(r,2))),sy);
  H[1][2]=H[2][1]=sub(reciprocal_pos(k),mul(R,J[0]));
  return H;
}

void verify_local(){
  const std::array<double,3> c={0.323144212181945529204914052954389108756940743,
                                0.654479031166666393676935731392503971225432437,
                                0.325165829719811989207067293084997587516390811};
  const double offsets[8]={-0.0105,-0.0075,-0.0045,-0.0015,0.0015,0.0045,0.0075,0.0105};
  const double rad=.00150000001;
  std::array<double,3> mins={INFINITY,INFINITY,INFINITY};size_t bad=0,cnt=0;
  for(double da:offsets)for(double dp:offsets)for(double dr:offsets){
    std::array<double,3> z={c[0]+da,c[1]+dp,c[2]+dr};auto H=hessian(z,rad);
    for(int i=0;i<3;i++)H[i][i]=sub(H[i][i],I("0.01"));
    I d1=H[0][0];
    I d2=sub(H[1][1],divi(square_interval(H[0][1]),d1));
    I cross=sub(H[1][2],divi(mul(H[0][1],H[0][2]),d1));
    I d3=sub(sub(H[2][2],divi(square_interval(H[0][2]),d1)),divi(square_interval(cross),d2));
    I ds[3]={d1,d2,d3};for(int i=0;i<3;i++)mins[i]=std::min(mins[i],lower_double(ds[i]));
    if(!positive(d1)||!positive(d2)||!positive(d3))bad++;cnt++;
  }
  std::cout<<"MPFR_LOCAL "<<cnt<<" pivots "<<std::setprecision(17)<<mins[0]<<" "<<mins[1]<<" "<<mins[2]<<" bad "<<bad<<"\n";
  if(bad)throw std::runtime_error("local check failed");
}

std::pair<I,I> primitive(const I&a,const I&z,unsigned N=700){
  I k=sub(one(),a),c=reciprocal_pos(k),H=zero(),zn=one();
  I L=add(add(H,reciprocal_pos(k)),neg(log_pos(z)));
  I s0=c,s1=mul(c,L);
  for(unsigned n=1;n<=N;n++){
    H=add(H,reciprocal_pos(add(a,I((unsigned long)(n-1)))));
    I nn((unsigned long)n),nm1((unsigned long)(n-1));
    c=divi(mul(mul(c,add(a,nm1)),add(k,nm1)),mul(nn,add(k,nn)));
    zn=mul(zn,z);
    L=add(add(H,reciprocal_pos(add(k,nn))),neg(log_pos(z)));
    s0=add(s0,mul(c,zn));s1=add(s1,mul(mul(c,zn),L));
  }
  I tail0=divi(mul(mul(c,zn),z),sub(one(),z));
  I term1=divi(mul(L,z),sub(one(),z));
  I term2=divi(z,mul(add(a,I((unsigned long)N)),pow_ui_pos(sub(one(),z),2)));
  I tail1=mul(mul(c,zn),add(term1,term2));
  I e0,e1;mpfr_set_ui(e0.l.x,0,MPFR_RNDN);e0.u=tail0.u;mpfr_set_ui(e1.l.x,0,MPFR_RNDN);e1.u=tail1.u;
  I zk=pow01(z,k);return {mul(zk,add(s0,e0)),mul(zk,add(s1,e1))};
}
I Tsmall(const I&a,const I&z,unsigned m){
  I k=sub(one(),a),w=neg(mul(k,log_pos(sub(one(),z))));
  if(m==0)return divi(sub(one(),expi(neg(w))),k);
  return divi(sub(one(),mul(expi(neg(w)),add(one(),w))),pow_ui_pos(k,2));
}
std::pair<I,std::array<I,3>> FG(const I&R,const I&a,const I&p,const I&r,
                               const I&A,const I&G,const I&Aa,const I&Ap,const I&Ga,const I&Gp,const I&Gr){
  I k=sub(one(),a),C=sub(one(),R);
  I F=add(add(mul(C,A),divi(mul(p,r),k)),neg(mul(R,G)));
  std::array<I,3> g={
    add(add(mul(C,Aa),divi(mul(p,r),pow_ui_pos(k,2))),neg(mul(R,Ga))),
    add(add(mul(C,Ap),divi(r,k)),neg(mul(R,Gp))),
    add(divi(p,k),neg(mul(R,Gr)))
  };
  return {F,g};
}

void verify_boundaries(){
  Checker ck;
  const double al=.1,pl=.02,rl=.00005,au=.999999,pu=.9;
  I R("0.8882517"),C=sub(one(),R),k=sub(one(),point(al));
  I K=mul(k,ck.coeff(al).B);
  std::array<std::pair<std::string,I>,6> tests={{{"small_a",sub(reciprocal_pos(K),R)},
    {"small_p",sub(C,mul(point(pl),sub(one(),log_pos(point(pl)))))},
    {"large_p",sub(point(pu),R)},
    {"small_r",sub(mul(C,add(sub(one(),point(pu)),mul(point(pu),log_pos(point(pu))))),mul(mul(R,point(rl)),sub(one(),log_pos(point(rl)))))},
    {"large_a",sub(divi(mul(point(pl),point(rl)),sub(one(),point(au))),R)},
    {"r_equals_one",C}}};
  for(auto &kv:tests){std::cout<<"MPFR_BOUNDARY "<<kv.first<<" margin "<<std::setprecision(17)<<lower_double(kv.second)<<"\n";if(!positive(kv.second))throw std::runtime_error("boundary failed");}
}

void verify_point(){
  I a("0.323144212181945529204914052954389108756940743");
  I p("0.654479031166666393676935731392503971225432437");
  I r("0.325165829719811989207067293084997587516390811");
  I k=sub(one(),a),x=pow01(p,reciprocal_pos(a)),y=pow01(r,reciprocal_pos(a));
  I Rlo("0.88825169032566246988432630129334211037");
  I Rhi("0.88825169032566246988432630129334211038");
  I Rc("0.8882517");
  auto ju=primitive(a,sub(one(),x)),jl=primitive(a,y);
  I J=sub(ju.first,jl.first),J1=sub(ju.second,jl.second);
  I A=divi(add(sub(k,p),mul(a,x)),k);
  I G=add(add(mul(p,Tsmall(a,y,0)),mul(r,Tsmall(a,x,0))),mul(mul(p,r),J));
  I Aa=neg(mul(p,Tsmall(a,sub(one(),x),1)));
  I Ap=neg(divi(sub(one(),pow01(x,k)),k));
  I Ga=add(add(mul(p,Tsmall(a,y,1)),mul(r,Tsmall(a,x,1))),mul(mul(p,r),J1));
  I Gp=add(Tsmall(a,y,0),mul(r,J)),Gr=add(Tsmall(a,x,0),mul(p,J));
  auto [FL,gL]=FG(Rlo,a,p,r,A,G,Aa,Ap,Ga,Gp,Gr);
  auto [FU,gU]=FG(Rhi,a,p,r,A,G,Aa,Ap,Ga,Gp,Gr);
  auto [FC,gC]=FG(Rc,a,p,r,A,G,Aa,Ap,Ga,Gp,Gr);
  I m("0.01"),sum=zero();for(auto &g:gL)sum=add(sum,square_interval(g));
  I lower=sub(FL,divi(sum,mul(I(2UL),m)));
  if(!positive(lower)||!negative(FU))throw std::runtime_error("point sign failed");
  I ginf=zero();for(auto &g:gC){I ag=abs_interval(g);if(lt(ginf.u,ag.u))ginf=ag;}
  I boundary=add(sub(FC,mul(I("0.037"),ginf)),divi(mul(m,I("0.00014161")),I(2UL)));
  if(!positive(boundary))throw std::runtime_error("boundary margin failed");
  I Q("1.983");if(mpfr_cmp_ui(add(pow01(I("0.667"),Q),pow01(I("0.338"),Q)).u.x,1)>=0)throw std::runtime_error("theta box failed");
  if(mpfr_cmp_ui(add(pow01(p,divi(k,a)),pow01(r,divi(k,a))).u.x,1)>=0)throw std::runtime_error("theta point failed");
  I maxg=zero();for(auto &g:gL){I ag=abs_interval(g);if(lt(maxg.u,ag.u))maxg=ag;}for(auto &g:gU){I ag=abs_interval(g);if(lt(maxg.u,ag.u))maxg=ag;}
  I loc=divi(mul(sqrti(I(3UL)),maxg),m);if(upper_double(loc)>=1e-30)throw std::runtime_error("location failed");
  I theta=reciprocal_pos(sub(I(2UL),add(pow01(p,divi(k,a)),pow01(r,divi(k,a)))));
  std::cout<<"MPFR_POINT lower "<<std::setprecision(17)<<lower_double(lower)
           <<" upper "<<upper_double(FU)<<" boundary "<<lower_double(boundary)
           <<" theta ["<<lower_double(theta)<<","<<upper_double(theta)<<"] location "<<upper_double(loc)<<" PASS\n";
}
}

int main(){
  try{verify_boundaries();verify_local();verify_point();std::cout<<"MPFR_AUXILIARY_ALL_PASS\n";return 0;}
  catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
