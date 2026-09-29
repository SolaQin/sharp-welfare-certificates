// Directed-rounding verifier. The input tree supplies no trusted values.
#include <mpfr.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <thread>
#include <atomic>
#include <mutex>
#include <cfenv>
#include <limits>
using namespace std;
double dn(double x){return nextafter(x,-INFINITY);}
double up(double x){return nextafter(x,INFINITY);}
mpfr_t mx,my,mz;
double expo(double x,mpfr_rnd_t r){mpfr_set_d(mx,x,MPFR_RNDN);mpfr_exp(mz,mx,r);return mpfr_get_d(mz,r);}
double power(double x,double y,mpfr_rnd_t r){mpfr_set_d(mx,x,MPFR_RNDN);mpfr_set_d(my,y,MPFR_RNDN);mpfr_pow(mz,mx,my,r);return mpfr_get_d(mz,r);}
struct I{double l,h;I(double x=0):l(x),h(x){}I(double a,double b):l(a),h(b){}};
I operator+(I x,I y){return {dn(x.l+y.l),up(x.h+y.h)};}
I operator-(I x,I y){return {dn(x.l-y.h),up(x.h-y.l)};}
I operator-(I x){return {-x.h,-x.l};}
I operator*(I x,I y){double v[]={x.l*y.l,x.l*y.h,x.h*y.l,x.h*y.h};return {dn(*min_element(v,v+4)),up(*max_element(v,v+4))};}
I operator/(I x,I y){if(y.l<=0&&y.h>=0)throw runtime_error("zero denominator");return x*I(dn(1/y.h),up(1/y.l));}
I logtwo;
I fast_exp_point(double x){
 // Elementary IEEE roundoff bounds; see Appendix C of the companion paper.
 constexpr double eps=0x1p-53,ln2=.6931471805599453;
 int n=int(round(x/.6931471805599453));
 double r=x-n*ln2,p=1;
 if(abs(r)>.35||abs(n)>1000000)throw runtime_error("exponential range reduction");
 for(int j=18;j>=1;j--)p=1+r*p/j;
 double err=eps*(100+2*abs(n));
 return {dn(ldexp(dn(p-err),n)),up(ldexp(up(p+err),n))};
}
I ex(I x){return {fast_exp_point(x.l).l,fast_exp_point(x.h).h};}
I logarithm_point(double x){
 if(x<=0)throw runtime_error("nonpositive logarithm");
 constexpr double eps=0x1p-53,ln2=.6931471805599453;
 int n;double y=frexp(x,&n);
 if(y<.7071067811865475){y*=2;n--;}
 double t=(y-1)/(y+1),t2=t*t,p=1./25;
 if(abs(t)>.172)throw runtime_error("logarithmic range reduction");
 for(int j=11;j>=0;j--)p=1./(2*j+1)+t2*p;
 double val=2*t*p+n*ln2,err=eps*(100+3*abs(n));
 return {dn(val-err),up(val+err)};
}
I pw(I x,double y){
 if(x.l<=0)throw runtime_error("nonpositive power base");
 return ex(I(logarithm_point(x.l).l,logarithm_point(x.h).h)*I(y));
}
I exprel(I x){
 double M=max(abs(x.l),abs(x.h));
 if(M<1e-4){I p=1,sum=1;double fac=1;for(int j=1;j<=10;j++){p=p*x;fac*=j+1;sum=sum+p/I(fac);}I rem=ex(I(M))*pw(I(max(M,1e-300)),11)/I(479001600.);return sum+I(-rem.h,rem.h);}
 return (ex(x)-I(1))/x;
}
struct Box{array<double,3> l,h;int depth;};
vector<double> us,ms,bs;
I rho,kappa;
I g(I w){return (I(1)+w)*ex(-w);}
I flow(double A,double B,double C,double s,double end,double u){
 I h=I(end)-I(s);h.l=max(0.,h.l);
 return ex(I(A)*h)*I(u)-I(B)*h*exprel(I(A)*h)-I(C)*ex(I(s))*h*ex(h)*exprel((I(A)-I(1))*h);
}
I stop(I s,double ell,double U){return s.h<=ell?I(ell)-s+I(1)-ex(I(ell)-I(U)):ex(I(ell)-s)-ex(I(ell)-I(U));}
double cutoff(double ell,double U){
 // Floating point only proposes a location. The final residual and the
 // analytic derivative bound f' >= kappa certify the returned upper bound.
 double r=(rho.l+rho.h)/2,l=0,h=U;
 for(int i=0;i<40;i++){
  double m=(l+h)/2,H=m<=ell?ell-m+1-exp(ell-U):exp(ell-m)-exp(ell-U);
  if((1-r)*m-r*H>0)h=m;else l=m;
 }
 double x=(l+h)/2;I f=kappa*I(x)-rho*stop(I(x),ell,U);
 return f.l>=0?x:min(U,(I(x)+I(-f.l)/kappa).h);
}
double approximate_flow(double A,double B,double C,double s,double end,double u){
 double h=end-s;
 auto er=[](double x){return abs(x)<1e-8?1+x/2+x*x/6:expm1(x)/x;};
 return exp(A*h)*u-B*h*er(A*h)-C*exp(s)*h*exp(h)*er((A-1)*h);
}
double before_event(double A,double B,double C,double s,double end,double u,double target,bool join){
 // An untrusted guess followed by a verified one-sided inequality.
 double l=s,h=end;
 for(int i=0;i<40;i++){
  double m=(l+h)/2,threshold=join?(1+target-m)*exp(-(target-m)):target;
  if(approximate_flow(A,B,C,s,m,u)>threshold)l=m;else h=m;
 }
 double offset=1e-10*(end-s)+1e-12;
 for(int i=0;i<18;i++){
  double candidate=max(s,l-offset);
  I threshold=join?g(I(target)-I(candidate)):I(target);
  if(flow(A,B,C,s,candidate,u).l>threshold.h)return candidate;
  if(candidate==s)return s;
  offset*=4;
 }
 return s;
}
double zlower(double alo,double ahi,double elo,double tcap,double limit){
 if(tcap<=0)return 0;
 double d=((I(1)-I(alo))*ex(I(alo)-I(elo))).h,s=0,u=1;
 for(size_t j=0;j<ms.size();){
  double A=(I(1)-I(ahi)*I(ms[j])).l;
  bool low=s<elo;
  double B=(I(ahi)*I(bs[j])+(low?I(0):I(d)*rho*ex(I(elo)))).h;
  double C=(I(d)*(low?I(1):kappa)).h;
  if((I(A)*I(u)-I(B)-I(C)*ex(I(s))).h>0)return -1;
  double end=min(tcap,limit);if(low)end=min(end,elo);
  I val=flow(A,B,C,s,end,u);double ue=val.l,next=us[j+1];
  bool knot=ue<next;
  if(knot){
   end=before_event(A,B,C,s,end,u,next,false);ue=next;
  }
  if(ue<=g(I(tcap)-I(end)).h){
   return before_event(A,B,C,s,end,u,tcap,true);
  }
  if(end>=limit)return limit;
  if(end<=s)return s;
  s=end;u=ue;if(knot)j++;
  if(end>=tcap)return end;
 }
 return 0;
}
bool verify(Box b){
 b.l[0]=max(b.l[0],b.l[1]);b.h[0]=min(b.h[0],b.h[2]);b.h[1]=min(b.h[1],b.h[0]);b.l[2]=max(b.l[2],b.l[0]);
 for(int i=0;i<3;i++)if(b.l[i]>b.h[i])return true;
 double al=b.l[0],ah=b.h[0],el=b.l[1],eh=b.h[1],Ul=b.l[2],Uh=b.h[2];
 double L=cutoff(eh,Uh),tcap=(I(Ul)-I(ah)).l;
 double z=zlower(al,ah,el,tcap,L);
 if(z<0)return false;
 if(z>=L)return true;
 if(L>=Ul)return false;
 auto fun=[&](I s){
  I v=I(Ul)-s;if(v.l<=0)return I(0,INFINITY);
  I Q=s.h<=eh?I(1):ex(I(eh)-s);Q.h=min(1.,Q.h);
  I factor;
  if(v.h<=1)factor=pw(v,-ah);
  else if(v.l>=1)factor=pw(v,-al);
  else{I p1=pw(v,-al),p2=pw(v,-ah);factor=I(min(p1.l,p2.l),max(p1.h,p2.h));}
  return (kappa+rho*Q)*factor;
 };
 auto integ=[&](double l,double h){
  if(h<=l)return I(0);
  int n=32;I delta=(I(h)-I(l))/I(n),ans=(fun(I(l))+fun(I(h)))/I(2);
  for(int j=1;j<n;j++)ans=ans+fun(I(l)+delta*I(j));
  return ans*delta;
 };
 I integral=integ(z,min(L,eh))+integ(max(z,eh),L);
 I base=I(Ul)-I(z),logbase(logarithm_point(base.l).l,logarithm_point(base.h).h),exponent;
 if(base.h<=1)exponent=I(1)-I(al);
 else if(base.l>=1)exponent=I(1)-I(ah);
 else exponent=I((I(1)-I(ah)).l,(I(1)-I(al)).h);
 I pos=ex(I(el)-I(Uh))*ex(logbase*exponent);
 return (pos-(I(1)-I(al))*integral).l>0;
}
int main(int argc,char**argv){
 static_assert(numeric_limits<double>::is_iec559&&numeric_limits<double>::digits==53,"requires IEEE binary64");
 if(fegetround()!=FE_TONEAREST)throw runtime_error("requires rounding to nearest");
 mpfr_init2(mx,53);mpfr_init2(my,53);mpfr_init2(mz,53);
 mpfr_const_log2(mx,MPFR_RNDD);logtwo.l=mpfr_get_d(mx,MPFR_RNDD);
 mpfr_const_log2(mx,MPFR_RNDU);logtwo.h=mpfr_get_d(mx,MPFR_RNDU);
 // The central constant's error is below eps/4 (checked at 200 bits).
 mpfr_t exact;mpfr_init2(exact,200);mpfr_const_log2(exact,MPFR_RNDU);
 mpfr_sub_d(exact,exact,.6931471805599453,MPFR_RNDU);
 if(mpfr_cmp_d(exact,0x1p-55)>=0)throw runtime_error("log(2) upper constant");
 mpfr_const_log2(exact,MPFR_RNDD);mpfr_sub_d(exact,exact,.6931471805599453,MPFR_RNDD);
 if(mpfr_sgn(exact)<=0)throw runtime_error("log(2) lower constant");mpfr_clear(exact);
 if(argc<2)return 1;ifstream file(argv[1],ios::binary);double raw;file.read(reinterpret_cast<char*>(&raw),sizeof(double));
 if(raw!=.90&&raw!=.91)throw runtime_error("unsupported target");rho=raw==.90?I(9)/I(10):I(91)/I(100);kappa=I(1)-rho;
 const double al=.0875,ah=.985;I upper=I(911390)/I(1000000);
 if((I(2)/(I(1)+ex(I(2)*I(al)))).l<=upper.h)throw runtime_error("small multiplier boundary");
 if((I(1)/(I(1)+(I(1)-I(ah))*(ex(I(2))-I(1)))).l<=upper.h)throw runtime_error("large multiplier boundary");
 if(((I(1)-upper)*I(2.5)-upper*ex(I(ah)-I(2.5))).l<=0)throw runtime_error("seller support boundary");
 if((I(1)-I(2.5)/I(29)).l<=upper.h)throw runtime_error("cap boundary");
 vector<double> ehi;
 double step=argc>2?stod(argv[2]):.1;int segments=int(30/step);
 for(int j=0;j<=segments;j++){
  double w=j*step,u=(1+w)*exp(-w);if(j==0)u=1;
  if(!isfinite(u)||u<=0||(j&&u>=us.back()))throw runtime_error("invalid polygon nodes");
  us.push_back(u);
  double l=0,h=40;
  for(int n=0;n<70;n++){double m=(l+h)/2;if(g(I(m)).l>u)l=m;else h=m;}
  ehi.push_back(j==0?1:ex(I(-l)).h);
 }
 for(int j=0;j<segments;j++){
  I slope=(I(ehi[j])-I(ehi[j+1]))/(I(us[j])-I(us[j+1]));
  I intercept=I(ehi[j])-slope*I(us[j]);ms.push_back(slope.h);bs.push_back(intercept.h);
 }
 vector<Box> stack{{{al,0,al},{ah,ah,29},0}};long long nodes=0,leaves=0,failed=0;
 auto start=chrono::steady_clock::now();
 int workers=argc>3?stoi(argv[3]):min(8u,max(1u,thread::hardware_concurrency()));
 vector<Box> batch;batch.reserve(8192);
 auto flush=[&](){
  atomic<size_t> next(0);atomic<long long> bad(0);mutex output;
  vector<thread> pool;
  for(int worker=0;worker<workers;worker++)pool.emplace_back([&](){
   for(;;){size_t j=next.fetch_add(1);if(j>=batch.size())break;
    if(!verify(batch[j])){long long ordinal=bad.fetch_add(1);
     if(failed+ordinal<6){lock_guard<mutex> lock(output);cout<<"FAILED ";for(int i=0;i<3;i++)cout<<batch[j].l[i]<<":"<<batch[j].h[i]<<" ";cout<<endl;}
    }
   }
  });
  for(auto& worker:pool)worker.join();failed+=bad.load();batch.clear();
 };
 while(!stack.empty()){
  Box b=stack.back();stack.pop_back();int label=file.get();
  // Permit a concurrently generated tree; an incomplete tree still fails.
  for(int retry=0;label==EOF&&retry<100;retry++){
   this_thread::sleep_for(chrono::milliseconds(100));file.clear();label=file.get();
  }
  if(label!=0&&label!=1)throw runtime_error("bad or incomplete tree");nodes++;
  if(label==1){
   leaves++;batch.push_back(b);if(batch.size()==8192)flush();
  }else{
   double scale[]={.2,.2,1};int axis=0;for(int i=1;i<3;i++)if((b.h[i]-b.l[i])/scale[i]>(b.h[axis]-b.l[axis])/scale[axis])axis=i;
   double m=(b.l[axis]+b.h[axis])/2;Box l=b,r=b;l.h[axis]=m;r.l[axis]=m;l.depth++;r.depth++;stack.push_back(r);stack.push_back(l);
  }
  if(nodes%200000==0)cout<<"visited "<<nodes<<" failed "<<failed<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 }
 flush();
 if(file.get()!=EOF)throw runtime_error("trailing tree data");
 cout<<(failed?"FAIL":"PASS")<<" target "<<raw<<" nodes "<<nodes<<" leaves "<<leaves<<" failed "<<failed<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 return failed?2:0;
}
