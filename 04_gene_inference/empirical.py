"""Exact finite-reference empirical P and BH utilities (no parametric tail extrapolation)."""
import math
def empirical_p(observed,null):
 if not null or any(not math.isfinite(x) for x in [observed,*null]):raise ValueError('Finite nonempty reference required')
 return (1+sum(abs(x)>=abs(observed) for x in null))/(len(null)+1)
def bh(p):
 if any(not math.isfinite(x) or x<0 or x>1 for x in p):raise ValueError('P must be in [0,1]')
 order=sorted(range(len(p)),key=p.__getitem__); out=[0.]*len(p);prev=1.
 for j in range(len(p)-1,-1,-1):
  i=order[j];prev=min(prev,p[i]*len(p)/(j+1));out[i]=prev
 return out
