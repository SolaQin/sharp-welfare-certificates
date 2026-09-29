#pragma once
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
// IEEE binary64 adjacency, exactly std::nextafter for finite inputs and infinities.
// memcpy avoids aliasing assumptions; NaNs retain the library's behavior.
inline double interval_nextafter(double x,double toward){
 if(std::isnan(x)||std::isnan(toward))return std::nextafter(x,toward);
 if(x==toward)return toward;
 std::uint64_t bits;std::memcpy(&bits,&x,sizeof bits);
 if((bits&UINT64_C(0x7fffffffffffffff))==0){
  bits=toward>0?UINT64_C(1):UINT64_C(0x8000000000000001);
 }else if((x<toward)==(x>0))++bits;else --bits;
 std::memcpy(&x,&bits,sizeof x);return x;
}
