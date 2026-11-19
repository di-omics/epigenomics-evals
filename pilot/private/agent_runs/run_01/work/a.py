import bisect
from collections import defaultdict
def load(p):
    d=defaultdict(list);n=0
    for l in open(p):
        l=l.rstrip('\n')
        if not l.strip(): continue
        f=l.split('\t');d[f[0]].append((int(f[1]),int(f[2])));n+=1
    return d,n
def ov(a,b):
    # a,b dicts; return count of a intervals overlapping any b, and set of b indices
    ca=0;hit=set()
    for c,iv in a.items():
        bl=sorted(b.get(c,[]))
        for s,e in iv:
            f=False
            for i,(bs,be) in enumerate(bl):
                if bs<e and s<be: f=True;hit.add((c,bs,be))
            ca+=f
    return ca,hit
import sys
A,na=load('inputs/peaks_A.bed');B,nb=load('inputs/peaks_B.bed');R,nr=load('inputs/reference_peaks.bed')
print(na,nb,nr)
# quadratic per chrom is ok? use simple sort-sweep via brute with bisect-free; sizes small enough
for nm,X,nx in(('A',A,na),('B',B,nb)):
    ca,_=ov(X,R); _,rh=ov(R,X); cr,_=ov(R,X)
    print(nm,ca,cr,ca/nx)
