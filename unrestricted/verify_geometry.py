"""Independently verify that the archived leaves form the claimed binary partition."""
import numpy as np
import os
from fractions import Fraction
from constants import CENTER_STR, RADIUS_STR
D=np.load(os.path.join(os.path.dirname(__file__),'cover.npz'))
arrays={t:D[t] for t in ['val','grad','local']}
lookup={}
for t in ['val','grad','local']:
    for i,node in enumerate(D[t+'_ids']):
        node=int(node)
        assert node not in lookup
        lookup[node]=(t,i)
center=D['center'];rad=float(D['rad']);root=D['root']
exact_center=list(map(Fraction,CENTER_STR))
exact_rad=Fraction(RADIUS_STR)
assert len(D['classes'])==len(D['grad'])
assert set(map(int,D['classes'])) <= {1,2,3,4}
assert all(len(D[t])==len(D[t+'_ids']) for t in arrays)
assert all(Fraction(float(root[i])) < exact_center[i]-exact_rad
           and exact_center[i]+exact_rad < Fraction(float(root[i+3]))
           for i in range(3))
stack=[(1,root.copy())];visited=set();steps=0
while stack:
    node,b=stack.pop();steps+=1
    assert steps<=2*len(lookup)+1
    if node in lookup:
        typ,j=lookup[node]
        assert np.array_equal(b,arrays[typ][j])
        if typ=='local':
            assert all(Fraction(float(b[i])) >= exact_center[i]-exact_rad
                       and Fraction(float(b[i+3])) <= exact_center[i]+exact_rad
                       for i in range(3))
        visited.add(node)
        continue
    w=b[3:]-b[:3];dim=int(np.argmax(w));mid=(b[dim]+b[dim+3])/2
    assert b[dim]<mid<b[dim+3]
    l=b.copy();u=b.copy();l[dim+3]=mid;u[dim]=mid
    stack.append((node*2+1,u));stack.append((node*2,l))
assert visited==set(lookup)
assert steps==2*len(lookup)-1
print('GEOMETRY PASS:',len(lookup),'leaves;',steps,'tree nodes; exact binary partition.')
