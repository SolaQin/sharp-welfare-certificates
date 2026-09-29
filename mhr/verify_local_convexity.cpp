// Validated parameter derivatives of the scalar residual.
// A successful point evaluation alone is not a global sharpness certificate.
#include "adjacent_float.hpp"
#define nextafter interval_nextafter
#define main archived_cover_main
#include "verify_global_cover.cpp"
#undef main
#undef nextafter
#include <iomanip>
#include "candidate_box.hpp"
namespace localcheck {
thread_local I rho,kappa;
thread_local int param_order=2,time_order=4;
thread_local int unnormalized=0;
thread_local bool physical_second=false;
constexpr int N=4, C=15, NT=5;
int ij(int i,int j){if(i>j)swap(i,j);return 5+i*4-i*(i-1)/2+j-i;}
I sq(I x){if(x.l>=0)return x*x;if(x.h<=0)return (-x)*(-x);return I(0,up(max(x.l*x.l,x.h*x.h)));}
I lg(I x){return I(logarithm_point(x.l).l,logarithm_point(x.h).h);}
I hull(I x,I y){return I(min(x.l,y.l),max(x.h,y.h));}
I intersect(I x,I y){I z(max(x.l,y.l),min(x.h,y.h));if(!isfinite(z.l)||!isfinite(z.h)||z.l>z.h)throw runtime_error("invalid intersection");return z;}
struct P{
 array<I,C> c{};
 P(double x=0){c[0]=I(x);} P(I x){c[0]=x;}
 static P var(I x,int i){P r(x);r.c[i+1]=1;return r;}
};
thread_local P center_result;
P operator+(P x,P y){for(int j=0;j<C;j++)x.c[j]=x.c[j]+y.c[j];return x;}
P operator-(P x,P y){for(int j=0;j<C;j++)x.c[j]=x.c[j]-y.c[j];return x;}
P operator-(P x){for(auto&v:x.c)v=-v;return x;}
P operator*(P x,P y){
 P z;z.c[0]=x.c[0]*y.c[0];
 for(int i=0;i<N;i++)z.c[i+1]=x.c[i+1]*y.c[0]+x.c[0]*y.c[i+1];
 if(param_order==2)for(int i=0;i<N;i++)for(int j=i;j<N;j++){
  int q=ij(i,j);z.c[q]=x.c[q]*y.c[0]+x.c[0]*y.c[q]+x.c[i+1]*y.c[j+1];
  if(i!=j)z.c[q]=z.c[q]+x.c[j+1]*y.c[i+1];
 }return z;
}
P operator*(P x,I y){for(auto&v:x.c)v=v*y;return x;}
P operator/(P x,I y){for(auto&v:x.c)v=v/y;return x;}
P compose(P x,I f,I d,I dd){
 P r(f);for(int i=0;i<N;i++)r.c[i+1]=d*x.c[i+1];
 if(param_order==2)for(int i=0;i<N;i++)for(int j=i;j<N;j++)
  r.c[ij(i,j)]=d*x.c[ij(i,j)]+dd*x.c[i+1]*x.c[j+1]/I(i==j?2:1);
 return r;
}
P inv(P x){I v=I(1)/x.c[0];return compose(x,v,-v*v,I(2)*v*v*v);}
P operator/(P x,P y){return x*inv(y);}
P expP(P x){I e=ex(x.c[0]);return compose(x,e,e,e);}
P logP(P x){I v=I(1)/x.c[0];return compose(x,lg(x.c[0]),v,-v*v);}
thread_local int ord=NT;
struct T{array<P,NT+1> c{};T(P x=P()){c[0]=x;}};
T operator+(T x,T y){for(int j=0;j<=ord;j++)x.c[j]=x.c[j]+y.c[j];return x;}
T operator-(T x,T y){for(int j=0;j<=ord;j++)x.c[j]=x.c[j]-y.c[j];return x;}
T operator-(T x){for(int j=0;j<=ord;j++)x.c[j]=-x.c[j];return x;}
T operator*(T x,T y){T z;for(int i=0;i<=ord;i++)for(int j=0;j<=i;j++)z.c[i]=z.c[i]+x.c[j]*y.c[i-j];return z;}
T inv(T x){T y(inv(x.c[0]));for(int i=1;i<=ord;i++){P q;for(int j=1;j<=i;j++)q=q+x.c[j]*y.c[i-j];y.c[i]=-q*y.c[0];}return y;}
T operator/(T x,T y){return x*inv(y);}
T expT(T x){T y(expP(x.c[0]));for(int i=1;i<=ord;i++){P q;for(int j=1;j<=i;j++)q=q+x.c[j]*I(j)*y.c[i-j];y.c[i]=q/I(i);}return y;}
T logT(T x){T y(logP(x.c[0])),v=inv(x);for(int i=1;i<=ord;i++){P q;for(int j=1;j<=i;j++)q=q+x.c[j]*I(j)*v.c[i-j];y.c[i]=q/I(i);}return y;}

template<class F>P integrate(F fn,int pieces){
 P sum;double h=1./pieces;
 for(int i=0;i<pieces;i++){
  double l=i*h,r=(i+1)*h,m=(l+r)/2;
  ord=0;P f0=fn(T(P(l))).c[0],fm=fn(T(P(m))).c[0],f1=fn(T(P(r))).c[0];
  T v(P(I(l,r)));v.c[1]=P(1);ord=4;P fourth=fn(v).c[4];
  sum=sum+(f0+fm*I(4)+f1)*I(h)/I(6)-fourth*pw(I(h),5)/I(120);
 }return sum;
}
P pre_equation(P a,P ell,P y,int pieces){
 P k=P(1)-a,r=a-ell,den=P(1)-expP(-(y+r));
 auto integrand=[&](T v){return expT(logT((T(P(1))-expT(-(T(y)*v+T(r))))/T(den))/T(k));};
 return y-ell-y*integrate(integrand,pieces);
}
double approximate_y(double a,double ell){
 double r=a-ell,k=1-a;
 auto fun=[&](double y){
  double den=1-exp(-y-r),sum=0;int n=128;
  for(int j=0;j<=n;j++){double v=double(j)/n;
   sum+=(j==0||j==n?1:j%2?4:2)*pow((1-exp(-y*v-r))/den,1/k);
  }return y-ell-y*sum/(3*n);
 };
 double lo=ell,hi=1;while(fun(hi)<0)hi*=2;
 for(int j=0;j<42;j++){double m=(lo+hi)/2;if(fun(m)>0)hi=m;else lo=m;}
 return (lo+hi)/2;
}
P crossing(P a,P ell,int pieces){
 I ai=a.c[0],ei=ell.c[0];
 auto enclose=[&](P av,P ev){
 I av0=av.c[0],ev0=ev.c[0];double yc=approximate_y((av0.l+av0.h)/2,(ev0.l+ev0.h)/2);
 double rad=3*(av0.h-av0.l)+4*(ev0.h-ev0.l)+1e-8;I yi;
 bool bracket=false;
 for(int j=0;j<15;j++){
  yi=I(max(ev0.l,yc-rad),yc+rad);
  I left=pre_equation(av,ev,P(yi.l),pieces).c[0],right=pre_equation(av,ev,P(yi.h),pieces).c[0];
  if(left.h<0&&right.l>0){bracket=true;break;}rad*=2;
 }
 if(!bracket)throw runtime_error("crossing bracket");
 for(int j=0;j<4;j++){
  P e=pre_equation(av,ev,P::var(yi,3),pieces);
  if(e.c[4].l<=0)break;
  double m=(yi.l+yi.h)/2;
  I em=pre_equation(av,ev,P(m),pieces).c[0];
  yi=intersect(yi,I(m)-em/e.c[4]);
 }return yi;
 };
 I yi;
 if(ai.h>ai.l||ei.h>ei.l){
  I yl=enclose(P(ai.l),P(ei.l)),yu=enclose(P(ai.h),P(ei.h));
  // On this range the crossing equation decreases in a and ell:
  // d_a log(M(v)/M(t)) = integral_v^t j(x)[k*j(x)-a]/k^2 dx >= 0.
  if((I(ai.h)+I(yu.h)-I(ei.l)).h<(-lg(I(ai.h))).l)yi=I(yl.l,yu.h);
  else yi=enclose(a,ell);
 }else yi=enclose(a,ell);
 // Use b=y-ell. Differentiation of the moving lower integration limit
 // supplies the ell derivatives directly, avoiding cancellation.
 I b=yi-ei,k=I(1)-ai,r=ai-ei,t=ai+b;
 I jr=I(1)/(ex(r)-I(1)),jt=I(1)/(ex(t)-I(1)),m=jt/k;
 I dlog=lg(I(1)-ex(-r))-lg(I(1)-ex(-t));
 I ratio=ex(dlog/k),gb=b*m;
 if(gb.l<=0)throw runtime_error("crossing derivative");
 P gaJet=pre_equation(P::var(ai,0),P(ei),P(yi),pieces);
 I ga=gaJet.c[1],gaa=I(2)*gaJet.c[ij(0,0)];
 I ma=m/k-ex(t)*jt*jt/k,mb=-ex(t)*jt*jt/k;
 I gab=-ga*m+b*ma,gbb=m-b*m*m+b*mb;
 I gale=-ratio*(dlog/(k*k)+(jr-jt)/k);
 I gble=ratio*m,glele=ratio*jr/k;
 I ba=-ga/gb,be=ratio/gb;
 P y(yi);y.c[1]=ba;y.c[2]=I(1)+be;
 y.c[ij(0,0)]=-(gaa+I(2)*gab*ba+gbb*ba*ba)/gb/I(2);
 y.c[ij(0,1)]=-(gale+gab*be+gble*ba+gbb*ba*be)/gb;
 y.c[ij(1,1)]=-(glele+I(2)*gble*be+gbb*be*be)/gb/I(2);
 return y;
}
T rhs(T t,T s,P a,P ell,P y0,P length){
 T y=T(y0)+T(length)*t;
 T den=T(P(1)-a)*(expT(y+T(a)-T(ell))*(T(P(1))-T(P(rho))*(T(P(1))-expT(T(ell)-s)))-T(P(1)));
 return T(length)*(y-s)/den;
}
T series(P s,I time,P a,P ell,P y0,P length,int n){
 T state(s),t{P(time)};t.c[1]=P(1);
 for(int j=0;j<n;j++){ord=j;T f=rhs(t,state,a,ell,y0,length);state.c[j+1]=f.c[j]/I(j+1);}
 return state;
}
bool inside(P x,P y){for(int j=0;j<C;j++){if(!isfinite(x.c[j].l)||!isfinite(x.c[j].h)||!isfinite(y.c[j].l)||!isfinite(y.c[j].h))throw runtime_error("nonfinite Picard bound");if(x.c[j].l<y.c[j].l||x.c[j].h>y.c[j].h)return false;}return true;}
P inflate(P x){
 for(auto&v:x.c){if(!isfinite(v.l)||!isfinite(v.h))throw runtime_error("nonfinite tube");double m=(v.l+v.h)/2,r=(v.h-v.l)*.65+1e-13;v=I(dn(m-r),up(m+r));}
 return x;
}
P flow_step(P s,double t,double dt,P a,P ell,P y0,P length){
  I times(t,t+dt);
  ord=0;P f=rhs(T(P(times)),T(s),a,ell,y0,length).c[0];
  P tube=inflate(s+f*I(0,dt));bool valid=false;
  for(int k=0;k<12;k++){
   ord=0;P f=rhs(T(P(times)),T(tube),a,ell,y0,length).c[0];
   P update=s+f*I(0,dt);
   if(inside(update,tube)){valid=true;break;}
   for(int q=0;q<C;q++)tube.c[q]=hull(tube.c[q],update.c[q]);tube=inflate(tube);
  }
  if(!valid)throw runtime_error("Picard enclosure");
  T coef=series(s,I(t),a,ell,y0,length,time_order),rem=series(tube,times,a,ell,y0,length,time_order+1);
  P next=rem.c[time_order+1];
  for(int k=time_order;k>=0;k--)next=coef.c[k]+next*I(dt);
  return next;
}
P tight(P box,P center,array<I,3> delta){
 if(param_order==1){
  I v=center.c[0];for(int i=0;i<3;i++)v=v+box.c[i+1]*delta[i];
  box.c[0]=intersect(box.c[0],v);return box;
 }
 for(int i=0;i<3;i++){
  I v=center.c[i+1];
  for(int j=0;j<3;j++)v=v+box.c[ij(i,j)]*I(i==j?2:1)*delta[j];
  box.c[i+1]=intersect(box.c[i+1],v);
 }
 I v=center.c[0];
 for(int i=0;i<3;i++)v=v+center.c[i+1]*delta[i];
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)
  v=v+box.c[ij(i,j)]*I(i==j?2:1)*delta[i]*delta[j]/I(2);
 box.c[0]=intersect(box.c[0],v);return box;
}
P advance(P s,double t,double dt,P a,P ell,P y0,P length){
  P phi=flow_step(P::var(s.c[0],3),t,dt,a,ell,y0,length);
  double mid=(s.c[0].l+s.c[0].h)/2;
  P center=flow_step(P(mid),t,dt,a,ell,y0,length),next;
  I J=phi.c[4];
  next.c[0]=intersect(phi.c[0],center.c[0]+J*(s.c[0]-I(mid)));
  for(int i=0;i<3;i++)next.c[i+1]=phi.c[i+1]+J*s.c[i+1];
  if(param_order==2)for(int i=0;i<3;i++)for(int q=i;q<3;q++){
   I h=phi.c[ij(i,q)]*I(i==q?2:1)+J*s.c[ij(i,q)]*I(i==q?2:1)
    +phi.c[ij(i,3)]*s.c[q+1]+phi.c[ij(q,3)]*s.c[i+1]
    +I(2)*phi.c[ij(3,3)]*s.c[i+1]*s.c[q+1];
   next.c[ij(i,q)]=h/I(i==q?2:1);
  }return next;
}
P post(P a,P ell,P y0,P end,int steps,P ac,P ec,P yc,P endc,array<I,3> delta,P& final_center){
 if(param_order==2&&!physical_second){
  P length=end-y0,lc=endc-yc;
  if(length.c[0].l<=0)throw runtime_error("head precedes lower endpoint");
  P s=ell,sc=ec;double dt=1./steps,t=0;
  while(t<1){
   dt=min(dt,1-t);
   try{
    P next=advance(s,t,dt,a,ell,y0,length),nc=advance(sc,t,dt,ac,ec,yc,lc);
    s=tight(next,nc,delta);sc=nc;t+=dt;dt=min(1./steps,dt*2);
   }catch(exception&e){if(dt<1e-7)throw;dt/=2;}
  }final_center=sc;return s;
 }
 double first=ceil(y0.c[0].h*steps)/steps,last=floor(end.c[0].l*steps)/steps;
 if(last<=first)throw runtime_error("head precedes lower endpoint");
 P s=ell,sc=ec;
 auto bridge=[&](P start,P length,P startc,P lengthc){
  double t=0,dt=.125;
  while(t<1){dt=min(dt,1-t);
   try{P next=advance(s,t,dt,a,ell,start,length),nc=advance(sc,t,dt,ac,ec,startc,lengthc);
    s=tight(next,nc,delta);sc=nc;t+=dt;dt=min(.125,dt*2);
   }catch(exception&e){if(dt<1e-7)throw;dt/=2;}
  }
 };
 bridge(y0,P(first)-y0,yc,P(first)-yc);
 double dt=1./steps,t=first;
 while(t<last){
  dt=min(dt,last-t);
  try{
   P next=advance(s,t,dt,a,ell,P(0),P(1));
   P nc=advance(sc,t,dt,ac,ec,P(0),P(1));
   s=tight(next,nc,delta);sc=nc;t+=dt;dt=min(1./steps,dt*2);
  }catch(exception&e){
   if(dt<1e-7)throw;dt/=2;
  }
 }
 bridge(P(last),end-P(last),P(last),endc-P(last));
 final_center=sc;return s;
}
I lambert(I x){
 auto endpoint=[](double x,bool upper){
  double w=log1p(x);for(int i=0;i<12;i++)w-=(w*exp(w)-x)/(exp(w)*(w+1));
  I err=I(w)*ex(I(w))-I(x); // derivative of w exp(w) is at least one.
  return upper?up(w+max(0.,-err.l)):max(0.,dn(w-max(0.,err.h)));
 };return I(endpoint(x.l,false),endpoint(x.h,true));
}
P endpoint(P ell,P U){
 double lo=0,hi=U.c[0].h;
 for(int j=0;j<50;j++){
  double m=(lo+hi)/2;I f=kappa*I(m)-rho*stop(I(m),ell.c[0].l,U.c[0].l);
  if(f.h<0)lo=m;else hi=m;
 }
 double upper=U.c[0].h,lower=0;
 for(int j=0;j<50;j++){
  double m=(lower+upper)/2;I f=kappa*I(m)-rho*stop(I(m),ell.c[0].h,U.c[0].h);
  if(f.l>0)upper=m;else lower=m;
 }
 I L(lo,upper);
 P lp=P::var(L,3),e=P(kappa)*lp-P(rho)*(expP(ell-lp)-expP(ell-U)),ans(L);I dy=e.c[4];
 for(int i=0;i<3;i++)ans.c[i+1]=-e.c[i+1]/dy;
 for(int i=0;i<3;i++)for(int j=i;j<3;j++){
  I h=e.c[ij(i,j)]*I(i==j?2:1)+e.c[ij(i,3)]*ans.c[j+1]
   +e.c[ij(j,3)]*ans.c[i+1]+I(2)*e.c[ij(3,3)]*ans.c[i+1]*ans.c[j+1];
  ans.c[ij(i,j)]=-h/dy/I(i==j?2:1);
 }return ans;
}
P evaluate(array<I,3> box,int pieces=16,int steps=64){
 P a=P::var(box[0],0),ell=P::var(box[1],1),U=P::var(box[2],2);
 P y0=crossing(a,ell,pieces),length=U-a-y0;
 array<I,3> delta,cb;
 for(int i=0;i<3;i++){double m=(box[i].l+box[i].h)/2;cb[i]=I(m);delta[i]=box[i]-I(m);}
 P ac=P::var(cb[0],0),ec=P::var(cb[1],1),uc=P::var(cb[2],2),yc=crossing(ac,ec,pieces);
 y0=tight(y0,yc,delta);length=U-a-y0;
 P zc,z=post(a,ell,y0,U-a,steps,ac,ec,yc,uc-ac,delta,zc),L=endpoint(ell,U),k=P(1)-a;
 if(getenv("MHR_DEBUG")){cerr<<"cross "<<y0.c[0].l<<" "<<y0.c[0].h<<" join "<<z.c[0].l<<" "<<z.c[0].h<<" L "<<L.c[0].l<<" "<<L.c[0].h<<endl;for(int i=1;i<4;i++)cerr<<"cross grad "<<i<<" "<<y0.c[i].l<<" "<<y0.c[i].h<<endl;for(int i=0;i<3;i++)for(int j=i;j<3;j++)cerr<<"H "<<i<<j<<" crossing "<<y0.c[ij(i,j)].l<<" "<<y0.c[ij(i,j)].h<<" join "<<z.c[ij(i,j)].l<<" "<<z.c[ij(i,j)].h<<endl;}
 if(z.c[0].l<=ell.c[0].h||L.c[0].l<=ell.c[0].h)throw runtime_error("tail crosses lower endpoint");
 auto complete=[&](P av,P ev,P uv,P zv,P lv){
  P kv=P(1)-av;
  auto integrand=[&](T v){
   T s=T(zv)*(T(P(1))-v)+T(lv)*v;
   return (T(P(1))-T(P(rho))*(T(P(1))-expT(T(ev)-s)))*expT(-T(av)*logT(T(uv)-s));
  };
  P tail=integrate(integrand,pieces);
  P norm=expP(ev-uv)*expP(kv*logP(uv-zv));
  if(unnormalized==1)return norm-kv*(lv-zv)*tail;
  if(unnormalized==2)return (norm-kv*(lv-zv)*tail)/expP(ev-uv);
  return P(1)-kv*(lv-zv)*tail/norm;
 };
 P result=complete(a,ell,U,z,L),center=complete(ac,ec,uc,zc,endpoint(ec,uc));
 center_result=center;
 return tight(result,center,delta);
}
bool positive_hessian(P f, double &margin){
 I aa=I(2)*f.c[ij(0,0)],bb=I(2)*f.c[ij(1,1)],cc=I(2)*f.c[ij(2,2)];
 I ab=f.c[ij(0,1)],ac=f.c[ij(0,2)],bc=f.c[ij(1,2)];
 if(aa.l<=0)return false;I d2=bb-sq(ab)/aa;if(d2.l<=0)return false;
 I d3=cc-sq(ac)/aa-sq(bc-ab*ac/aa)/d2;
 margin=d3.l;return d3.l>0;
}
}
#ifndef MHR_LOCAL_LIBRARY
int main(int argc,char**argv){
 using namespace localcheck;
 if(fegetround()!=FE_TONEAREST)throw runtime_error("rounding mode");
 localcheck::rho=I(91138936)/I(100000000);localcheck::kappa=I(1)-localcheck::rho;
 array<I,3> point{I(.28523064649549684),I(.04945697286361709),I(1.1746687514380845)};
 int pieces=argc>1?stoi(argv[1]):16,steps=argc>2?stoi(argv[2]):64;
 if(argc>3){
  double rad=stod(argv[3]);for(auto &x:point)x=I(x.l-rad,x.h+rad);
  localcheck::rho=I((I(91138936)/I(100000000)).l,(I(91138938)/I(100000000)).h);localcheck::kappa=I(1)-localcheck::rho;
 }
 try{
  auto start=chrono::steady_clock::now();P f=evaluate(point,pieces,steps);
  cout<<setprecision(17)<<"residual "<<f.c[0].l<<" "<<f.c[0].h<<endl;
  for(int j=0;j<3;j++)cout<<"gradient "<<j<<" "<<f.c[j+1].l<<" "<<f.c[j+1].h<<endl;
  for(int i=0;i<3;i++)for(int j=i;j<3;j++){I h=f.c[ij(i,j)]*I(i==j?2:1);cout<<"Hessian "<<i<<j<<" "<<h.l<<" "<<h.h<<endl;}
  double margin;cout<<"positive "<<positive_hessian(f,margin)<<" pivot "<<margin<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 }catch(exception&e){cout<<"INCONCLUSIVE "<<e.what()<<endl;return 2;}
}
#endif
