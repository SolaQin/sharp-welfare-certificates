// Fixed-score first sensitivities, adapted from centered_certificate/high_order_score.hpp.
// This version returns the potential R/q and its derivatives.
namespace directscore {
I hull(I x,I y){return I(min(x.l,y.l),max(x.h,y.h));}
I meet(I x,I y){I z(max(x.l,y.l),min(x.h,y.h));if(z.l>z.h||!isfinite(z.l)||!isfinite(z.h))throw runtime_error("empty enclosure");return z;}
I lg(I x){if(x.l<=0)throw runtime_error("logarithm");return I(logarithm_point(x.l).l,logarithm_point(x.h).h);}
I powi(I x,I y){return ex(lg(x)*y);}
double mag(I x){return max(abs(x.l),abs(x.h));}
struct Param{I a,e,U,r;};
struct Field{I f,fx,fs;array<I,3> fp;I third;};
Field branch(I x,I s,const Param&p,bool above){
 I k=I(1)-p.a,h=p.U-p.a,y=h*x,w=y-s;w.l=max(0.,w.l);
 I E=ex(y+p.a-p.e),D,Ds,De,Eb;
 if(!above){Eb=E;D=k*(E-I(1));Ds=I(0);De=-k*E;}
 else{I EQ=ex(p.a+w);Eb=(I(1)-p.r)*E+p.r*EQ;D=k*(Eb-I(1));Ds=-k*p.r*EQ;De=-k*(I(1)-p.r)*E;}
 if(D.l<=0)throw runtime_error("denominator");
 I f=w/D,fs=-(I(1)+f*Ds)/D,fy=(I(1)-f*(D+k))/D;
 I fa=-f*(I(1)-p.a*Eb)/D,fe=-f*De/D;
 I along=fy+fs*f,Dp=D+k+Ds*f,Dpp=D+k+Ds*(I(2)*f-f*f+along);
 I third=h*h*h*(-along/D-I(2)*(I(1)-f)*Dp/(D*D)-f*Dpp/D+I(2)*f*Dp*Dp/(D*D));
 return {h*f,h*h*fy,h*fs,{-f+h*(fa-x*fy),h*fe,f+h*x*fy},third};
}
Field field(I x,I s,const Param&p){
 if(s.h<=p.e.l)return branch(x,s,p,false);
 if(s.l>=p.e.h)return branch(x,s,p,true);
 auto a=branch(x,s,p,false),b=branch(x,s,p,true);
 a.f=hull(a.f,b.f);a.fx=hull(a.fx,b.fx);a.fs=hull(a.fs,b.fs);
 for(int i=0;i<3;i++)a.fp[i]=hull(a.fp[i],b.fp[i]);a.third=hull(a.third,b.third);return a;
}
struct Step{I next;Field coefficients;bool smooth;};
Step step(I s,double t,double dt,const Param&p){
 I X(t,t+dt),h=p.U-p.a,Y=h*X;
 Field first=field(X,I(s.l,min(Y.h,s.h+dt*max(1.,h.h))),p);
 double top=min(Y.h,(s+I(dt)*first.f).h);
 I tube(s.l,up(top+1e-12));bool ok=false;Field c;
 for(int iter=0;iter<20;iter++){
  c=field(X,tube,p);I reach=s+I(0,dt)*c.f;
  reach.h=min(reach.h,Y.h);
  if(reach.h<=tube.h){ok=true;break;}
  tube.h=up(reach.h+max(1e-12,(reach.h-tube.h)*.25));
 }
 if(!ok)throw runtime_error("state tube");
 if((I(1)+I(dt)*c.fs).l<=0)throw runtime_error("Euler map not monotone");
 bool smooth=tube.h<=p.e.l||tube.l>=p.e.h;
 I remainder;
 if(smooth){Field at=field(I(t),s,p);remainder=I(dt)*I(dt)/I(2)*(at.fx+at.fs*at.f)+I(dt)*I(dt)*I(dt)/I(6)*c.third;}
 else remainder=I(dt)*I(dt)/I(2)*(c.fx+c.fs*c.f);
 I low=I(s.l)+I(dt)*field(I(t),I(s.l),p).f+remainder;
 I high=I(s.h)+I(dt)*field(I(t),I(s.h),p).f+remainder;
 I next(max(s.l,low.l),min((h*I(t+dt)).h,high.h));
 if(next.l>next.h||!isfinite(next.h))throw runtime_error("state step");
 return {next,c,smooth};
}
I point_step(I s,double t,double dt,const Param&p){
 auto result=step(s,t,dt,p);
 if(result.smooth||dt<=0x1p-17)return result.next;
 I mid=point_step(s,t,dt/2,p);return point_step(mid,t+dt/2,dt/2,p);
}
I linear(I v,I alpha,I beta,double dt){
 I x=alpha*I(dt);
 // exprel is increasing, so endpoint evaluation avoids division across zero.
 I er(exprel(I(x.l)).l,exprel(I(x.h)).h);
 return ex(x)*v+I(dt)*er*beta;
}
struct Head{I z,zc;array<I,3> dz;};
Head head(Param p,int n){
 Param pc=p;array<I,3> delta;
 I* vals[]={&p.a,&p.e,&p.U};I* cvs[]={&pc.a,&pc.e,&pc.U};
 for(int j=0;j<3;j++){double m=(vals[j]->l+vals[j]->h)/2;*cvs[j]=I(m);delta[j]=*vals[j]-I(m);}
 I s(0),sc(0);array<I,3> ds{I(0),I(0),I(0)};double dt=1./n;
 for(int i=0;i<n;i++){
  double t=i*dt;auto b=step(s,t,dt,p);I cn=sc;
  for(int j=0;j<2;j++)cn=point_step(cn,t+j*dt/2,dt/2,pc);
  for(int j=0;j<3;j++)ds[j]=linear(ds[j],b.coefficients.fs,b.coefficients.fp[j],dt);
  I centered=cn;for(int j=0;j<3;j++)centered=centered+ds[j]*delta[j];
  s=meet(b.next,centered);sc=cn;
 }return {s,sc,ds};
}
I endpoint(Param p){
 auto at=[&](double ell,double U){
  double lo=0,hi=U,r=(p.r.l+p.r.h)/2;
  for(int i=0;i<50;i++){double m=(lo+hi)/2,H=m<=ell?ell-m+1-exp(ell-U):exp(ell-m)-exp(ell-U);if((1-r)*m-r*H<0)lo=m;else hi=m;}
  double m=(lo+hi)/2;I f=(I(1)-p.r)*I(m)-p.r*stop(I(m),ell,U);
  I correction=f/(I(1)-p.r);
  return I(max(0.,dn(m-max(0.,correction.h))),min(U,up(m+max(0.,-correction.l))));
 };
 return I(at(p.e.l,p.U.l).l,at(p.e.h,p.U.h).h);
}
struct Tail{I integral,ilog,i1,ib;};
Tail tail(Param p,I z,I L,int n){
 if(z.l<=p.e.h||L.l<=z.h||L.h>=p.U.l)throw runtime_error("tail geometry");
 I len=L-z,k=I(1)-p.r;
 auto terms=[&](I s){
  I v=p.U-s,Q=ex(p.e-s),b=k+p.r*Q,power=powi(v,-p.a);
  return array<I,4>{b*power,b*power*lg(v),b*power/v,p.r*Q*power};
 };
 I sim=terms(z)[0]+terms(L)[0];
 for(int j=1;j<n;j++)sim=sim+terms(z*(I(1)-I(j)/I(n))+L*I(j)/I(n))[0]*I(j%2?4:2);
 array<I,4> mid{};
 for(int j=0;j<n;j++){auto term=terms(z*(I(1)-(I(j)+I(.5))/I(n))+L*(I(j)+I(.5))/I(n));for(int q=1;q<4;q++)mid[q]=mid[q]+term[q];}
 I s(z.l,L.h),v=p.U-s,Q=ex(p.e-s),b=k+p.r*Q,A=p.a,lv=lg(v),power=powi(v,-A);
 I rising=I(1),poly=I(1);int binom[]={1,4,6,4,1};
 for(int j=1;j<=4;j++){rising=rising*(A+I(j-1));poly=poly+I((j%2?-1:1)*binom[j])*rising/powi(v,I(j));}
 I fourth=power*(k*rising/powi(v,I(4))+p.r*Q*poly);
 I factor2=powi(len,I(3))/(I(24)*I(n)*I(n));
 I d2log=power*(p.r*Q*lv-I(2)*p.r*Q*(A*lv-I(1))/v+b*(A*(A+I(1))*lv-I(2)*A-I(1))/(v*v));
 I d2one=power/v*(p.r*Q-I(2)*p.r*Q*(A+I(1))/v+b*(A+I(1))*(A+I(2))/(v*v));
 I d2b=p.r*Q*power*(I(1)-I(2)*A/v+A*(A+I(1))/(v*v));
 return {sim*len/(I(3)*I(n))-powi(len,I(5))*fourth/(I(180)*powi(I(n),I(4))),
  mid[1]*len/I(n)+factor2*d2log,mid[2]*len/I(n)+factor2*d2one,mid[3]*len/I(n)+factor2*d2b};
}
struct Result{I f,center;array<I,3> gradient;I cap;};
Result evaluate(Param p,int steps=256,int pieces=32){
 Head h=head(p,steps);I L=endpoint(p),z=h.z,k=I(1)-p.a,q=ex(p.e-p.U),W=p.U-z,d=p.U-L;
 Tail t=tail(p,z,L,pieces);I scale=q*powi(W,k),R=scale-k*t.integral,f=R/q;
 I QL=ex(p.e-L),bL=I(1)-p.r+p.r*QL,bz=I(1)-p.r+p.r*ex(p.e-z);
 I Le=p.r*(QL-q)/bL,LU=p.r*q/bL,C=k*(bz-q)*powi(W,-p.a),bd=bL*powi(d,-p.a);
 array<I,3> Rp{
  -scale*lg(W)+t.integral+k*t.ilog+C*h.dz[0],
  scale-k*(bd*Le+t.ib)+C*h.dz[1],
  -scale+k*q*powi(W,-p.a)-k*bd*LU+k*p.a*t.i1+C*h.dz[2]};
 array<I,3> slog{-lg(W)-k*h.dz[0]/W,I(1)-k*h.dz[1]/W,-I(1)+k*(I(1)-h.dz[2])/W};
 array<I,3> gradient{Rp[0]/q,(Rp[1]-R)/q,(Rp[2]+R)/q};
 Param pc=p;pc.a=I((p.a.l+p.a.h)/2);pc.e=I((p.e.l+p.e.h)/2);pc.U=I((p.U.l+p.U.h)/2);
 I Lc=endpoint(pc);Tail tc=tail(pc,h.zc,Lc,pieces);
 I fc=powi(pc.U-h.zc,I(1)-pc.a)-(I(1)-pc.a)*tc.integral/ex(pc.e-pc.U);
 // Pull the algebraic cap condition back through the enclosed head flow.
 using localcheck::P;
 int old_order=localcheck::param_order;localcheck::param_order=1;
 auto H=[&](P av,P ev,P uv,P zv,P lv){
  P kv=P(1)-av,Wv=uv-zv,dv=uv-lv,c=P(I(1)-p.r)/localcheck::expP(ev-uv);
  return Wv+P(1)-P(2)*av-kv*P(p.r)*localcheck::expP(Wv)
   +c*(localcheck::expP(av*localcheck::logP(Wv/dv))*(dv+kv*(lv+P(1)))-Wv-kv);
 };
 P zj(z),Lj(L);for(int j=0;j<3;j++)zj.c[j+1]=h.dz[j];Lj.c[2]=Le;Lj.c[3]=LU;
 P Hj=H(P::var(p.a,0),P::var(p.e,1),P::var(p.U,2),zj,Lj);
 I cap=H(P(pc.a),P(pc.e),P(pc.U),P(h.zc),P(Lc)).c[0];
 array<I,3> delta{p.a-pc.a,p.e-pc.e,p.U-pc.U};
 for(int j=0;j<3;j++)cap=cap+Hj.c[j+1]*delta[j];
 cap=meet(cap,Hj.c[0]);localcheck::param_order=old_order;
 return {f,fc,gradient,cap};
}
}
