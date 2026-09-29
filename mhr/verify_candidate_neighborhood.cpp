// Local checks only. A completed global exclusion certificate is also needed.
#define MHR_LOCAL_LIBRARY
#include "verify_local_convexity.cpp"
int main(){
 using localcheck::P;
 cout<<setprecision(17);
 localcheck::rho=I((I(9113893679.)/I(10000000000.)).l,(I(9113893684.)/I(10000000000.)).h);
 localcheck::kappa=I(1)-localcheck::rho;
 array<I,3> box;for(int i=0;i<3;i++)box[i]=I(candidate_box::lower[i],candidate_box::upper[i]);
 localcheck::physical_second=true;
 try{
  P f=localcheck::evaluate(box,128,256);double margin=0;
  bool ok=localcheck::positive_hessian(f,margin);
  cout<<"LOCAL HESSIAN "<<(ok?"PASS":"INCONCLUSIVE")<<" final pivot "<<margin<<endl;
  for(int i=0;i<3;i++)for(int j=i;j<3;j++){I h=f.c[localcheck::ij(i,j)]*I(i==j?2:1);cout<<"H "<<i<<j<<" "<<h.l<<" "<<h.h<<endl;}
  if(!ok)return 2;
  double scales[]={1,1,10},strong=INFINITY;
  for(int i=0;i<3;i++){
   I row=I(2)*f.c[localcheck::ij(i,i)]*I(scales[i]*scales[i]);
   for(int j=0;j<3;j++)if(i!=j){I h=f.c[localcheck::ij(i,j)]*I(scales[i]*scales[j]);row=row-I(max(abs(h.l),abs(h.h)));}
   strong=min(strong,row.l);
  }
  cout<<"SCALED STRONG CONVEXITY "<<strong<<endl;
  if(strong<=0)return 2;
  localcheck::physical_second=false;
  array<I,3> point{I(.28523064649549684),I(.04945697286361709),I(1.1746687514380845)};
  for(long long numerator:{9113893680LL,9113893683LL}){
   localcheck::rho=I(double(numerator))/I(10000000000.);
   localcheck::kappa=I(1)-localcheck::rho;
   P p=localcheck::evaluate(point,256,512);
   cout<<"POINT "<<numerator<<" residual "<<p.c[0].l<<" "<<p.c[0].h<<endl;
   if(numerator==9113893680LL){
    I norm=0;for(int i=0;i<3;i++)norm=norm+localcheck::sq(p.c[i+1]*I(scales[i]));
    I lower=p.c[0]-norm/(I(2)*I(strong));
    cout<<"LOCAL LOWER RESIDUAL "<<lower.l<<endl;if(lower.l<=0)return 2;
   }else if(p.c[0].h>=0)return 2;
  }
  for(long long numerator:{911389368124LL,911389368127LL}){
   localcheck::rho=I(double(numerator))/I(1000000000000.);
   localcheck::kappa=I(1)-localcheck::rho;
   P p=localcheck::evaluate(point,256,512);
   cout<<"REFINED POINT "<<numerator<<" / 1000000000000 residual "<<p.c[0].l<<" "<<p.c[0].h<<endl;
   if(numerator==911389368124LL){
    I norm=0;for(int i=0;i<3;i++)norm=norm+localcheck::sq(p.c[i+1]*I(scales[i]));
    I lower=p.c[0]-norm/(I(2)*I(strong));
    cout<<"REFINED LOCAL LOWER RESIDUAL "<<lower.l<<endl;if(lower.l<=0)return 2;
   }else if(p.c[0].h>=0)return 2;
  }
  cout<<"PASS LOCAL CHECKS ONLY"<<endl;
 }catch(exception&e){cout<<"INCONCLUSIVE "<<e.what()<<endl;return 2;}
}
