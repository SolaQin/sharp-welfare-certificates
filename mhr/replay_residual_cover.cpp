// Independent leaf replay: paths and method labels supply no trusted values.
#define MHR_LOCAL_LIBRARY
#include "verify_local_convexity.cpp"
#include "fast_atom_enclosure.hpp"
#include "fast_sensitivity.hpp"
#include "interval_newton.hpp"
struct Leaf{int label;string path;};
array<I,3> decode(const string&path){
 array<I,3>b{I(.19,.39),I(.002,.16),I(.85,1.65)};double scale[]={.2,.158,.8};
 for(char c:path){int d=0;for(int i=1;i<3;i++)if((b[i].h-b[i].l)/scale[i]>(b[d].h-b[d].l)/scale[d])d=i;
  double m=(b[d].l+b[d].h)/2;if(c=='0')b[d].h=m;else if(c=='1')b[d].l=m;else throw runtime_error("invalid path");}
 return b;
}
bool local_box(array<I,3>b){const double* lo=candidate_box::lower;const double* hi=candidate_box::upper;for(int i=0;i<3;i++)if(b[i].l<lo[i]||b[i].h>hi[i])return false;return true;}
bool jets(array<I,3>b){
 double w=max({b[0].h-b[0].l,b[1].h-b[1].l,(b[2].h-b[2].l)*.25});
 int pieces=w<1e-5?128:w<1e-3?16:8,steps=w<1e-5?256:w<1e-3?64:32;
 for(int mode:{0,1,2,3})try{
  int order=mode==3?2:mode==2?1:(w<.004?2:1);localcheck::physical_second=mode!=1;
  localcheck::param_order=order;localcheck::time_order=w<1e-4?4:3;localcheck::unnormalized=2;
  auto f=localcheck::evaluate(b,pieces,steps);bool pass=false;
  for(int i=0;i<4;i++){
   if(!isfinite(f.c[i].l)||!isfinite(f.c[i].h))throw runtime_error("nonfinite jet");
   if(f.c[i].l>0||(i&&f.c[i].h<0))pass=true;
  }
  if(pass)return true;
  if(order==2&&localcheck::stationary_certificate(f,b))return true;
 }catch(exception&e){}
 return false;
}
bool sensitivity(array<I,3>b){
 double w=max({b[0].h-b[0].l,b[1].h-b[1].l,(b[2].h-b[2].l)*.25});
 for(int n:{w>.002?128:256,2048})try{
  auto f=directscore::evaluate({b[0],b[1],b[2],localcheck::rho},n,64);bool pass=false;I bound=f.center;
  if(!isfinite(f.f.l)||!isfinite(f.f.h))throw runtime_error("nonfinite value");if(f.f.l>0)pass=true;
  for(int i=0;i<3;i++){
   I g=f.gradient[i];if(!isfinite(g.l)||!isfinite(g.h))throw runtime_error("nonfinite gradient");
   if(g.l>0||g.h<0)pass=true;bound=bound+g*(b[i]-I((b[i].l+b[i].h)/2));
  }if(pass||bound.l>0||(isfinite(f.cap.l)&&isfinite(f.cap.h)&&(f.cap.l>0||f.cap.h<0)))return true;
  if(w>=.00015)break;
 }catch(exception&e){}
 return false;
}
bool check(const Leaf&leaf){
 auto b=decode(leaf.path);
 if(local_box(b))return true;
 localcheck::rho=I((I(9113893680.)/I(10000000000.)).l,(I(9113893683.)/I(10000000000.)).h);
 localcheck::kappa=I(1)-localcheck::rho;
 if(leaf.label==2)return local_box(b);
 if(leaf.label==0){Box raw{{b[0].l,b[1].l,b[2].l},{b[0].h,b[1].h,b[2].h},0};return verify(raw);}
 if(leaf.label==4||leaf.label==5)try{
  auto v=atomcheck::enclose(b,8,32);I f=localcheck::kappa-v.atom*(I(1)-ex(b[1]-b[0]));
  if(isfinite(f.l)&&isfinite(f.h)&&(f.l>0||f.h<0))return true;
  if(isfinite(v.cap.l)&&isfinite(v.cap.h)&&(v.cap.l>0||v.cap.h<0))return true;
 }catch(exception&e){}
 if(leaf.label==6&&sensitivity(b))return true;
 if((leaf.label==1||leaf.label==7||leaf.label==8)&&sensitivity(b))return true;
 if(jets(b))return true;
 return leaf.label!=6&&sensitivity(b);
}
int main(int argc,char**argv){
 static_assert(numeric_limits<double>::is_iec559&&numeric_limits<double>::digits==53,"requires binary64");
 if(fegetround()!=FE_TONEAREST)throw runtime_error("rounding mode");if(argc<2)return 1;
 ::rho=I(911390)/I(1000000);::kappa=I(1)-::rho;
 vector<double>ehi;double spacing=.05;int segments=int(30/spacing);
 for(int j=0;j<=segments;j++){
  double w=j*spacing,u=(1+w)*exp(-w);if(j==0)u=1;
  if(!isfinite(u)||u<=0||(j&&u>=us.back()))throw runtime_error("polygon nodes");us.push_back(u);
  double l=0,h=40;for(int n=0;n<70;n++){double m=(l+h)/2;if(g(I(m)).l>u)l=m;else h=m;}
  ehi.push_back(j==0?1:ex(I(-l)).h);
 }
 for(int j=0;j<segments;j++){
  I slope=(I(ehi[j])-I(ehi[j+1]))/(I(us[j])-I(us[j+1]));
  I intercept=I(ehi[j])-slope*I(us[j]);ms.push_back(slope.h);bs.push_back(intercept.h);
 }
 vector<Leaf>leaves;ifstream input(argv[1]);Leaf l;
 while(input>>l.label>>l.path){if(l.label!=0&&l.label!=1&&l.label!=2&&l.label!=4&&l.label!=5&&l.label!=6&&l.label!=7&&l.label!=8)throw runtime_error("unknown method");leaves.push_back(l);}
 if(!input.eof()||leaves.empty())throw runtime_error("bad or empty certificate");
 atomic<size_t>next(0),done(0),failed(0);mutex output;auto start=chrono::steady_clock::now();
 int workers=argc>2?stoi(argv[2]):8;vector<thread>pool;
 for(int i=0;i<workers;i++)pool.emplace_back([&]{for(;;){size_t j=next.fetch_add(1);if(j>=leaves.size())return;
  bool ok=false;try{ok=check(leaves[j]);}catch(exception&e){}
  if(!ok){size_t f=failed.fetch_add(1);if(f<5){lock_guard<mutex>lock(output);cout<<"FAILED "<<leaves[j].label<<" "<<leaves[j].path<<endl;}}
  size_t n=done.fetch_add(1)+1;if(n%2000==0){lock_guard<mutex>lock(output);cout<<"checked "<<n<<" failed "<<failed.load()<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;}
 }});
 for(auto&t:pool)t.join();cout<<(failed.load()?"FAIL":"PASS")<<" checked "<<done.load()<<" failed "<<failed.load()<<" seconds "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 return failed.load()?2:0;
}
