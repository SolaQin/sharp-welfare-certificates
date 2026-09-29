// Fixed preconditioning of the mean-value equation for a stationary point.
// Returns 7 for no stationary point, 8 if all are in the certified local box.
namespace localcheck {
int stationary_certificate(P f,array<I,3> box){
 double C[3][3]={{.268177,.07958937,.86465812},{.07958937,.72821789,.72790652},{.86465812,.72790652,21.69919562}};
 array<I,3> delta,center,beta{};array<array<I,3>,3>B{};
 for(int i=0;i<3;i++){center[i]=I((box[i].l+box[i].h)/2);delta[i]=box[i]-center[i];
  for(int j=0;j<3;j++){
   I g=center_result.c[j+1];if(!isfinite(g.l)||!isfinite(g.h))throw runtime_error("nonfinite center gradient");
   beta[i]=beta[i]-I(C[i][j])*g;
   B[i][j]=I(i==j?1:0);
   for(int k=0;k<3;k++){
    I H=f.c[ij(k,j)]*I(k==j?2:1);if(!isfinite(H.l)||!isfinite(H.h))throw runtime_error("nonfinite Hessian");
    B[i][j]=B[i][j]-I(C[i][k])*H;
   }
  }
 }
 const double* nl=candidate_box::lower;const double* nh=candidate_box::upper;
 for(int iter=0;iter<24;iter++){
  array<I,3> next;bool inside=true;
  for(int i=0;i<3;i++){
   I image=beta[i];for(int j=0;j<3;j++)image=image+B[i][j]*delta[j];
   if(!isfinite(image.l)||!isfinite(image.h))throw runtime_error("nonfinite Newton enclosure");
   if(image.h<delta[i].l||image.l>delta[i].h)return 7;
   next[i]=intersect(image,delta[i]);I located=center[i]+next[i];
   if(located.l<nl[i]||located.h>nh[i])inside=false;
  }
  if(inside)return 8;
  bool unchanged=true;for(int i=0;i<3;i++)if(next[i].l!=delta[i].l||next[i].h!=delta[i].h)unchanged=false;
  if(unchanged)break;delta=next;
 }
 return 0;
}
}
