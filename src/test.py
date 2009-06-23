import numpy,math
from random import *

nbins = 1000
lmax = 20
dtau = 1
p_B = 10.0

def gfn(dtau):
  wght = numpy.zeros(lmax)
  for L in range(lmax):
    wght[L] = (2.0*L+1.0)*math.exp(-dtau*p_B*L*(L+1))
  wght = wght / (4.0*math.pi)
  
  P = numpy.zeros(lmax)
  gfn = numpy.zeros(nbins)
  P[0] = 1.0
  
  for i in range(nbins):
    x = 2.0*(float(i)/float(nbins-1)) - 1.0
    P[1] = x
    for L in range(1,lmax-1):
      P[L+1] = ((2.0*L+1.0)*x*P[L] - L*P[L-1]) / (L+1)
  
    gfn[i] = abs(sum(P*wght))
    print x,gfn[i]

for i in range(1,10):
  dtau = 1e-2*float(1.0/i)
  gfn(dtau)
  print
