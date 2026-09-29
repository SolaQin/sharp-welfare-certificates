import mpmath as mp
import numpy as np
import os
from verify_cover import coeff
mp.mp.dps=60;mp.iv.dps=50;iv=mp.iv
D=np.load(os.path.join(os.path.dirname(__file__),'cover.npz'))
al,pl,rl,au,pu,ru=[iv.mpf(float(x))for x in D['root']]
R=iv.mpf('0.8882517');C=1-R
k=1-al;K=k*coeff(float(D['root'][0]))[2]
tests={
 'small_a':1/K-R,
 'small_p':C-pl*(1-iv.ln(pl)),
 'large_p':pu-R,
 'small_r':C*(1-pu+pu*iv.ln(pu))-R*rl*(1-iv.ln(rl)),
 'large_a':pl*rl/(1-au)-R,
 'r_equals_one':C,
}
for name,v in tests.items():
    assert mp.mpf(v._mpi_[0])>0
    print(name, v)
assert float(ru.a)==1.
print('ALL EXTERIOR-BOUNDARY CHECKS PASS')
