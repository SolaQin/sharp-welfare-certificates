#include <gmp.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Minimal ABI declarations for MPFR 4.x.  The archived environment ships
// libmpfr.so.6 but not the development header.  On a normal development
// installation this block may be replaced by #include <mpfr.h>.
extern "C" {
typedef long int mpfr_prec_t;
typedef long int mpfr_exp_t;
typedef int mpfr_sign_t;
typedef enum {
  MPFR_RNDN = 0, MPFR_RNDZ = 1, MPFR_RNDU = 2, MPFR_RNDD = 3,
  MPFR_RNDA = 4, MPFR_RNDF = 5, MPFR_RNDNA = -1
} mpfr_rnd_t;
typedef struct {
  mpfr_prec_t _mpfr_prec;
  mpfr_sign_t _mpfr_sign;
  mpfr_exp_t _mpfr_exp;
  mp_limb_t *_mpfr_d;
} __mpfr_struct;
typedef __mpfr_struct mpfr_t[1];
typedef __mpfr_struct *mpfr_ptr;
typedef const __mpfr_struct *mpfr_srcptr;
void mpfr_init2(mpfr_ptr, mpfr_prec_t);
void mpfr_clear(mpfr_ptr);
int mpfr_set(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_set_d(mpfr_ptr, double, mpfr_rnd_t);
int mpfr_set_ui(mpfr_ptr, unsigned long, mpfr_rnd_t);
int mpfr_set_str(mpfr_ptr, const char *, int, mpfr_rnd_t);
int mpfr_add(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_sub(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_mul(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_div(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_ui_div(mpfr_ptr, unsigned long, mpfr_srcptr, mpfr_rnd_t);
int mpfr_pow(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_pow_ui(mpfr_ptr, mpfr_srcptr, unsigned long, mpfr_rnd_t);
int mpfr_log(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_neg(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int mpfr_cmp(mpfr_srcptr, mpfr_srcptr);
int mpfr_cmp_ui(mpfr_srcptr, unsigned long);
double mpfr_get_d(mpfr_srcptr, mpfr_rnd_t);
}

namespace {
constexpr mpfr_prec_t PREC = 160;
constexpr int SERIES_N = 40;

struct Real {
  mpfr_t x;
  Real() { mpfr_init2(x, PREC); mpfr_set_ui(x, 0, MPFR_RNDN); }
  explicit Real(double v) : Real() { mpfr_set_d(x, v, MPFR_RNDN); }
  Real(const Real &o) : Real() { mpfr_set(x, o.x, MPFR_RNDN); }
  Real(Real &&o) noexcept : Real() { mpfr_set(x, o.x, MPFR_RNDN); }
  Real &operator=(const Real &o) {
    if (this != &o) mpfr_set(x, o.x, MPFR_RNDN);
    return *this;
  }
  Real &operator=(Real &&o) noexcept {
    if (this != &o) mpfr_set(x, o.x, MPFR_RNDN);
    return *this;
  }
  ~Real() { mpfr_clear(x); }
};

struct I {
  Real l, u;
  I() = default;
  explicit I(double v) {
    mpfr_set_d(l.x, v, MPFR_RNDD);
    mpfr_set_d(u.x, v, MPFR_RNDU);
  }
  explicit I(unsigned long v) {
    mpfr_set_ui(l.x, v, MPFR_RNDD);
    mpfr_set_ui(u.x, v, MPFR_RNDU);
  }
  explicit I(const char *s) {
    if (mpfr_set_str(l.x, s, 10, MPFR_RNDD) ||
        mpfr_set_str(u.x, s, 10, MPFR_RNDU))
      throw std::runtime_error("invalid decimal constant");
  }
};

I from_bounds(const Real &l, const Real &u) { I z; z.l = l; z.u = u; return z; }
I point(double x) { return I(x); }
I zero() { return I(0UL); }
I one() { return I(1UL); }

bool lt(const Real &a, const Real &b) { return mpfr_cmp(a.x, b.x) < 0; }
bool positive(const I &a) { return mpfr_cmp_ui(a.l.x, 0) > 0; }
bool negative(const I &a) { return mpfr_cmp_ui(a.u.x, 0) < 0; }
double lower_double(const I &a) { return mpfr_get_d(a.l.x, MPFR_RNDD); }
double upper_double(const I &a) { return mpfr_get_d(a.u.x, MPFR_RNDU); }

I add(const I &a, const I &b) {
  I z;
  mpfr_add(z.l.x, a.l.x, b.l.x, MPFR_RNDD);
  mpfr_add(z.u.x, a.u.x, b.u.x, MPFR_RNDU);
  return z;
}
I sub(const I &a, const I &b) {
  I z;
  mpfr_sub(z.l.x, a.l.x, b.u.x, MPFR_RNDD);
  mpfr_sub(z.u.x, a.u.x, b.l.x, MPFR_RNDU);
  return z;
}
I neg(const I &a) {
  I z;
  mpfr_neg(z.l.x, a.u.x, MPFR_RNDD);
  mpfr_neg(z.u.x, a.l.x, MPFR_RNDU);
  return z;
}
Real min_real(const std::array<Real,4> &v) {
  Real r = v[0];
  for (int i=1;i<4;i++) if (lt(v[i],r)) r=v[i];
  return r;
}
Real max_real(const std::array<Real,4> &v) {
  Real r = v[0];
  for (int i=1;i<4;i++) if (lt(r,v[i])) r=v[i];
  return r;
}
I mul(const I &a, const I &b) {
  std::array<Real,4> dl,du;
  mpfr_mul(dl[0].x,a.l.x,b.l.x,MPFR_RNDD);
  mpfr_mul(dl[1].x,a.l.x,b.u.x,MPFR_RNDD);
  mpfr_mul(dl[2].x,a.u.x,b.l.x,MPFR_RNDD);
  mpfr_mul(dl[3].x,a.u.x,b.u.x,MPFR_RNDD);
  mpfr_mul(du[0].x,a.l.x,b.l.x,MPFR_RNDU);
  mpfr_mul(du[1].x,a.l.x,b.u.x,MPFR_RNDU);
  mpfr_mul(du[2].x,a.u.x,b.l.x,MPFR_RNDU);
  mpfr_mul(du[3].x,a.u.x,b.u.x,MPFR_RNDU);
  return from_bounds(min_real(dl),max_real(du));
}
I reciprocal_pos(const I &b) {
  if (mpfr_cmp_ui(b.l.x,0)<=0) throw std::runtime_error("nonpositive denominator");
  I z;
  mpfr_ui_div(z.l.x,1,b.u.x,MPFR_RNDD);
  mpfr_ui_div(z.u.x,1,b.l.x,MPFR_RNDU);
  return z;
}
I divi(const I &a, const I &b) { return mul(a,reciprocal_pos(b)); }
I pow01(const I &base, const I &expo) {
  if (mpfr_cmp_ui(base.l.x,0)<0 || mpfr_cmp_ui(base.u.x,1)>0 ||
      mpfr_cmp_ui(expo.l.x,0)<=0)
    throw std::runtime_error("pow01 domain");
  I z;
  mpfr_pow(z.l.x,base.l.x,expo.u.x,MPFR_RNDD);
  mpfr_pow(z.u.x,base.u.x,expo.l.x,MPFR_RNDU);
  return z;
}
I pow_ui_pos(const I &base, unsigned long n) {
  if (mpfr_cmp_ui(base.l.x,0)<0) throw std::runtime_error("pow_ui domain");
  I z;
  mpfr_pow_ui(z.l.x,base.l.x,n,MPFR_RNDD);
  mpfr_pow_ui(z.u.x,base.u.x,n,MPFR_RNDU);
  return z;
}
I log_pos(const I &a) {
  if (mpfr_cmp_ui(a.l.x,0)<=0) throw std::runtime_error("log domain");
  I z;
  mpfr_log(z.l.x,a.l.x,MPFR_RNDD);
  mpfr_log(z.u.x,a.u.x,MPFR_RNDU);
  return z;
}
I hull(const I &a, const I &b) {
  Real l=lt(a.l,b.l)?a.l:b.l;
  Real u=lt(a.u,b.u)?b.u:a.u;
  return from_bounds(l,u);
}
I pospart(const I &a) {
  I z;
  if (mpfr_cmp_ui(a.l.x,0)>0) z.l=a.l; else mpfr_set_ui(z.l.x,0,MPFR_RNDN);
  if (mpfr_cmp_ui(a.u.x,0)>0) z.u=a.u; else mpfr_set_ui(z.u.x,0,MPFR_RNDN);
  return z;
}
I min_interval(const I &a,const I &b) {
  Real l=lt(a.l,b.l)?a.l:b.l;
  Real u=lt(a.u,b.u)?a.u:b.u;
  return from_bounds(l,u);
}
I max_interval(const I &a,const I &b) {
  Real l=lt(a.l,b.l)?b.l:a.l;
  Real u=lt(a.u,b.u)?b.u:a.u;
  return from_bounds(l,u);
}
I clip01(const I &a) {
  I z;
  if (mpfr_cmp_ui(a.l.x,0)<0) mpfr_set_ui(z.l.x,0,MPFR_RNDN); else z.l=a.l;
  if (mpfr_cmp_ui(a.u.x,1)>0) mpfr_set_ui(z.u.x,1,MPFR_RNDN); else z.u=a.u;
  return z;
}

struct Coeff { I a,k,B; std::vector<I> c; };
uint64_t bits_of(double x) { uint64_t u; std::memcpy(&u,&x,sizeof u); return u; }

class Checker {
 public:
  Checker() : R("0.8882517"), C(sub(one(),R)) {}
  const Coeff &coeff(double af) {
    uint64_t key=bits_of(af);
    auto it=cache.find(key);
    if(it!=cache.end()) return *it->second;
    auto cc=std::make_unique<Coeff>();
    cc->a=point(af); cc->k=sub(one(),cc->a); cc->c.reserve(SERIES_N+1);
    cc->c.push_back(reciprocal_pos(cc->k));
    for(int n=1;n<=SERIES_N;n++) {
      I nn(static_cast<unsigned long>(n));
      I nm1(static_cast<unsigned long>(n-1));
      I num=mul(mul(add(cc->a,nm1),add(cc->k,nm1)),cc->c.back());
      I den=mul(nn,add(cc->k,nn));
      cc->c.push_back(divi(num,den));
    }
    I z("0.5"),P=cc->c.back();
    for(int n=SERIES_N-1;n>=0;n--) P=add(mul(P,z),cc->c[n]);
    I tail=divi(mul(cc->c.back(),pow_ui_pos(z,SERIES_N+1)),sub(one(),z));
    I err; mpfr_set_ui(err.l.x,0,MPFR_RNDN); err.u=tail.u;
    cc->B=mul(I(2UL),mul(pow01(z,cc->k),add(P,err)));
    auto [jt,ok]=cache.emplace(key,std::move(cc));
    return *jt->second;
  }
  I Bz(double af,I z) {
    const Coeff &cc=coeff(af); z=clip01(z);
    I half("0.5"); bool comp=lt(half.u,z.l);
    if(comp) z=sub(one(),z);
    I P=cc.c.back();
    for(int n=SERIES_N-1;n>=0;n--) P=add(mul(P,z),cc.c[n]);
    I tail=divi(mul(cc.c.back(),pow_ui_pos(z,SERIES_N+1)),sub(one(),z));
    I err; mpfr_set_ui(err.l.x,0,MPFR_RNDN); err.u=tail.u;
    I result=mul(pow01(z,cc.k),add(P,err));
    return comp?sub(cc.B,result):result;
  }
  I powx(double af,double p) { return pow01(point(p),reciprocal_pos(point(af))); }
  I A(double af,double p) {
    I a=point(af),k=sub(one(),a),pp=point(p);
    return divi(add(sub(k,pp),mul(a,powx(af,p))),k);
  }
  I T(double af,I z) {
    I k=sub(one(),point(af)); z=clip01(z);
    return divi(sub(one(),pow01(sub(one(),z),k)),k);
  }
  I L(double af,double p,double r) {
    I a=point(af),k=sub(one(),a),pp=point(p),rr=point(r);
    I x=powx(af,p),y=powx(af,r),sum=add(x,y);
    I v1=add(add(divi(T(af,y),rr),divi(T(af,x),pp)),
              sub(Bz(af,sub(one(),x)),Bz(af,y)));
    I v2=divi(sub(divi(sub(rr,mul(a,y)),k),A(af,p)),mul(pp,rr));
    if(mpfr_cmp_ui(sum.u.x,1)<0) return v1;
    if(mpfr_cmp_ui(sum.l.x,1)>0) return v2;
    return hull(v1,v2);
  }
  I Ip(double af,double r,I smax) {
    I y=powx(af,r);smax=clip01(smax);
    return add(T(af,min_interval(y,smax)),
               mul(point(r),pospart(sub(Bz(af,smax),Bz(af,y)))));
  }
  I Ir(double af,double p,I ymin) {
    I smax=sub(one(),powx(af,p));ymin=clip01(ymin);
    I start=max_interval(ymin,smax),k=sub(one(),point(af));
    return add(divi(sub(one(),pow01(start,k)),k),
               mul(point(p),pospart(sub(Bz(af,smax),Bz(af,ymin)))));
  }
  std::pair<bool,double> check_val(const std::array<double,6>&b) {
    auto [al,pl,rl,au,pu,ru]=b;
    I v=add(add(divi(mul(C,A(au,pu)),mul(point(pu),point(ru))),
                reciprocal_pos(sub(one(),point(al)))),
            neg(mul(R,L(au,pl,rl))));
    return {positive(v),lower_double(v)};
  }
  std::pair<bool,double> check_grad(const std::array<double,6>&b,int typ) {
    auto [al,pl,rl,au,pu,ru]=b; I q;
    if(typ==3) {
      q=sub(divi(point(pl),sub(one(),point(al))),mul(R,Ir(au,pu,powx(al,rl))));
      return {positive(q),lower_double(q)};
    }
    if(typ==4) {
      q=sub(divi(point(pu),sub(one(),point(au))),mul(R,Ir(al,pl,powx(au,ru))));
      return {negative(q),-upper_double(q)};
    }
    if(typ==1) {
      I exponent=divi(sub(one(),point(al)),point(al));
      I ls=divi(sub(one(),pow01(point(pl),exponent)),sub(one(),point(au)));
      q=add(add(neg(mul(C,ls)),divi(point(rl),sub(one(),point(al)))),
            neg(mul(R,Ip(au,ru,sub(one(),powx(al,pl))))));
      return {positive(q),lower_double(q)};
    }
    if(typ==2) {
      I exponent=divi(sub(one(),point(au)),point(au));
      I ls=divi(sub(one(),pow01(point(pu),exponent)),sub(one(),point(al)));
      q=add(add(neg(mul(C,ls)),divi(point(ru),sub(one(),point(au)))),
            neg(mul(R,Ip(al,rl,sub(one(),powx(au,pu))))));
      return {negative(q),-upper_double(q)};
    }
    throw std::runtime_error("invalid derivative class");
  }
 private:
  I R,C;
  std::unordered_map<uint64_t,std::unique_ptr<Coeff>> cache;
};

struct Data {
  std::vector<std::array<double,6>> val,grad;
  std::vector<unsigned char> cls;
};
Data read_data(const std::string &path) {
  std::ifstream f(path,std::ios::binary);
  if(!f) throw std::runtime_error("cannot open data");
  char magic[8];f.read(magic,8);
  if(std::memcmp(magic,"SBWCOV1\0",8)!=0) throw std::runtime_error("bad magic");
  uint64_t nv=0,ng=0;f.read(reinterpret_cast<char*>(&nv),8);f.read(reinterpret_cast<char*>(&ng),8);
  Data d;d.val.resize(nv);d.grad.resize(ng);d.cls.resize(ng);
  auto read_boxes=[&](auto &arr){
    for(auto &b:arr)for(double &x:b){uint64_t u;f.read(reinterpret_cast<char*>(&u),8);std::memcpy(&x,&u,8);}
  };
  read_boxes(d.val);read_boxes(d.grad);
  f.read(reinterpret_cast<char*>(d.cls.data()),static_cast<std::streamsize>(ng));
  if(!f) throw std::runtime_error("truncated data");
  return d;
}
}

int main(int argc,char**argv) {
  try {
    std::string path=argc>1?argv[1]:"cover_boxes.bin";
    size_t limit=argc>2?std::stoull(argv[2]):0;
    Data d=read_data(path);Checker ck;
    size_t nval=limit?std::min(limit,d.val.size()):d.val.size();
    size_t ngrad=limit?std::min(limit,d.grad.size()):d.grad.size();
    size_t bad=0,done=0;double margin=std::numeric_limits<double>::infinity();
    auto t0=std::chrono::steady_clock::now();
    for(size_t i=0;i<nval;i++){
      auto [ok,m]=ck.check_val(d.val[i]);margin=std::min(margin,m);
      if(!ok){bad++;std::cerr<<"bad val "<<i<<" margin "<<m<<"\n";}done++;
      if(done%10000==0)std::cout<<"verified "<<done<<" min_margin "<<std::setprecision(17)<<margin<<" bad "<<bad<<"\n";
    }
    for(size_t i=0;i<ngrad;i++){
      auto [ok,m]=ck.check_grad(d.grad[i],d.cls[i]);margin=std::min(margin,m);
      if(!ok){bad++;std::cerr<<"bad grad "<<i<<" class "<<int(d.cls[i])<<" margin "<<m<<"\n";}done++;
      if(done%10000==0)std::cout<<"verified "<<done<<" min_margin "<<std::setprecision(17)<<margin<<" bad "<<bad<<"\n";
    }
    double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::cout<<"MPFR_FINAL "<<done<<" min_margin "<<std::setprecision(17)<<margin<<" bad "<<bad<<" seconds "<<sec<<"\n";
    return bad?1:0;
  } catch(const std::exception&e){
    std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;
  }
}
