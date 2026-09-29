// Scalar interval Taylor flow; an off-root floor derivative certificate for R.
namespace atomcheck {
constexpr int NN=5;
thread_local int vord=NN;
struct V {array<I,NN+1>c{};V(I x=I()){c[0]=x;}};
V operator+(V x,V y){for(int j=0;j<=vord;j++)x.c[j]=x.c[j]+y.c[j];return x;}
V operator-(V x,V y){for(int j=0;j<=vord;j++)x.c[j]=x.c[j]-y.c[j];return x;}
V operator-(V x){for(auto&v:x.c)v=-v;return x;}
V operator*(V x,V y){V z;for(int i=0;i<=vord;i++)for(int j=0;j<=i;j++)z.c[i]=z.c[i]+x.c[j]*y.c[i-j];return z;}
V inv(V x){V y(I(1)/x.c[0]);for(int i=1;i<=vord;i++){I q;for(int j=1;j<=i;j++)q=q+x.c[j]*y.c[i-j];y.c[i]=-q*y.c[0];}return y;}
V operator/(V x,V y){return x*inv(y);}
V ev(V x){V y(ex(x.c[0]));for(int i=1;i<=vord;i++){I q;for(int j=1;j<=i;j++)q=q+I(j)*x.c[j]*y.c[i-j];y.c[i]=q/I(i);}return y;}
V lv(V x){V y(localcheck::lg(x.c[0])),v=inv(x);for(int i=1;i<=vord;i++){I q;for(int j=1;j<=i;j++)q=q+I(j)*x.c[j]*v.c[i-j];y.c[i]=q/I(i);}return y;}
struct Pair{V s,j;};
Pair rhs(V t,V s,I a,I ell,I y0,I len){
 V y=V(y0)+V(len)*t;
 V den=V(I(1)-a)*(ev(y+V(a)-V(ell))*(V(I(1))-V(localcheck::rho)*(V(I(1))-ev(V(ell)-s)))-V(I(1)));
 return {V(len)*(y-s)/den,V(len*a)/den};
}
Pair series(I s,I t,I a,I ell,I y0,I len,int order){
 V state(s),time(t);time.c[1]=I(1);
 for(int j=0;j<order;j++){vord=j;state.c[j+1]=rhs(time,state,a,ell,y0,len).s.c[j]/I(j+1);}
 vord=4;
 return {state,rhs(time,state,a,ell,y0,len).j};
}
I inflate(I x){double m=(x.l+x.h)/2,r=(x.h-x.l)*.6+1e-12;return I(dn(m-r),up(m+r));}
Pair step(I s,double t,double dt,I a,I ell,I y0,I len){
 vord=0;I times(t,t+dt);I f=rhs(V(times),V(s),a,ell,y0,len).s.c[0];
 I tube=inflate(s+I(0,dt)*f);bool ok=false;
 for(int j=0;j<12;j++){
  I update=s+I(0,dt)*rhs(V(times),V(tube),a,ell,y0,len).s.c[0];
  if(isfinite(update.l)&&isfinite(update.h)&&update.l>=tube.l&&update.h<=tube.h){ok=true;break;}
  tube=inflate(localcheck::hull(update,tube));
 }
 if(!ok)throw runtime_error("atom flow tube");
 auto cf=series(s,I(t),a,ell,y0,len,4),rem=series(tube,times,a,ell,y0,len,5);
 I next=rem.s.c[5];for(int j=4;j>=0;j--)next=cf.s.c[j]+next*I(dt);
 I jinc=rem.j.c[4]/I(5);for(int j=3;j>=0;j--)jinc=cf.j.c[j]/I(j+1)+jinc*I(dt);jinc=jinc*I(dt);
 // Dependence on the scalar initial state is bounded by the variational flow.
 I y=y0+len*times,w=y-tube;
 I den=(I(1)-a)*(ex(y+a-ell)*(I(1)-localcheck::rho*(I(1)-ex(ell-tube)))-I(1));
 I fs=len*(-den+w*(I(1)-a)*localcheck::rho*ex(y+a-tube))/(den*den);
 I jac=ex(I(dt)*fs);
 double mid=(s.l+s.h)/2;
 auto cc=series(I(mid),I(t),a,ell,y0,len,4);
 I cn=rem.s.c[5];for(int j=4;j>=0;j--)cn=cc.s.c[j]+cn*I(dt);
 next=localcheck::intersect(next,cn+jac*(s-I(mid)));
 return {V(next),V(jinc)};
}
struct Summary{I atom,z,L,cap;};
array<I,3> tail_integrals(I a,I ell,I U,I z,I L,int pieces){
 array<I,3>sum{};double dh=1./pieces;
 auto fun=[&](V v){
  V s=V(z)*(V(I(1))-v)+V(L)*v,w=V(U)-s;
  V logw=lv(w),base=(V(I(1))-V(localcheck::rho)*(V(I(1))-ev(V(ell)-s)))*ev(-V(a+I(1))*logw);
  return array<V,3>{base,base/w,base*logw};
 };
 for(int i=0;i<pieces;i++){
  double l=i*dh,r=(i+1)*dh,m=(l+r)/2;
  vord=0;auto fl=fun(V(I(l))),fm=fun(V(I(m))),fr=fun(V(I(r)));
  V cell(I(l,r));cell.c[1]=I(1);vord=4;auto fc=fun(cell);
  for(int j=0;j<3;j++)sum[j]=sum[j]+(fl[j].c[0]+I(4)*fm[j].c[0]+fr[j].c[0])*I(dh)/I(6)-fc[j].c[4]*pw(I(dh),5)/I(120);
 }
 for(auto&v:sum)v=v*(L-z);return sum;
}
Summary enclose(array<I,3>box,int pieces=8,int steps=32,bool force_cap=false){
 using localcheck::P;I a=box[0],ell=box[1],U=box[2],k=I(1)-a;
 int old=localcheck::param_order;localcheck::param_order=1;
 I y0=localcheck::crossing(P(a),P(ell),pieces).c[0];
 I len=U-a-y0;if(len.l<=0)throw runtime_error("atom head below floor");
 I s=ell,J=0;double t=0,dt=1./steps;
 while(t<1){
  dt=min(dt,1-t);
  try{auto q=step(s,t,dt,a,ell,y0,len);s=q.s.c[0];J=J+q.j.c[0];t+=dt;dt=min(1./steps,dt*2);}
  catch(exception&e){if(dt<1e-5)throw;dt/=2;}
 }
 I L=localcheck::endpoint(P(ell),P(U)).c[0];
 if(s.l<=ell.h||s.h>=L.l||L.h>=U.l)throw runtime_error("atom phase");
 I logpre=(localcheck::lg(I(1)-ex(ell-a))-localcheck::lg(I(1)-ex(ell-a-y0)))*a/k;
 I ans=ex(a*localcheck::lg((U-L)/(U-s))+logpre-J);
 I floor=localcheck::kappa-ans*(I(1)-ex(ell-a));
 if(!force_cap&&isfinite(floor.l)&&isfinite(floor.h)&&(floor.l>0||floor.h<0)){
  localcheck::param_order=old;return {ans,s,L,I(-INFINITY,INFINITY)};
 }
 // F=R/q has an algebraic cap derivative W^{-a} H.
 localcheck::param_order=1;
 auto cap_function=[&](P av,P ev,P uv,P zv){
  P kv=P(1)-av,lv=localcheck::endpoint(ev,uv),w=uv-zv,dv=uv-lv;
  P c=P(localcheck::kappa)/localcheck::expP(ev-uv);
  return w+P(1)-P(2)*av-kv*P(localcheck::rho)*localcheck::expP(w)
   +c*(localcheck::expP(av*localcheck::logP(w/dv))*(dv+kv*(lv+P(1)))-w-kv);
 };
 array<I,4> xx{a,ell,U,s},center,delta;
 for(int i=0;i<4;i++){center[i]=I((xx[i].l+xx[i].h)/2);delta[i]=xx[i]-center[i];}
 P capjet=cap_function(P::var(a,0),P::var(ell,1),P::var(U,2),P::var(s,3));
 I value=cap_function(P(center[0]),P(center[1]),P(center[2]),P(center[3])).c[0];
 for(int i=0;i<4;i++)value=value+capjet.c[i+1]*delta[i];
 I cap=localcheck::intersect(capjet.c[0],value);
 localcheck::param_order=old;return {ans,s,L,cap};
}
I atom(array<I,3>box,int pieces=8,int steps=32){return enclose(box,pieces,steps).atom;}
}
