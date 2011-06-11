#! Just a title
#!---------------
from numpy import *
import scipy.special as sp
import matplotlib.pyplot as P


G = 4.306e-6
rho0 =  1.144e7 
Sigma0 = 6.594e8
RD = 2.63 #in kpc
rc = 6.414 #in kpc

def gdisk(r):
    y = r / (2 * RD)
    return 4 * pi * G * Sigma0 * (RD / r) * (y ** 2) * (sp.i0(y) * sp.k0(y) - sp.i1(y) * sp.k1(y))

def ghalo(r):
    return 4 * pi * G * rho0 / r * ((rc ** 2) - ((rc ** 3) / r) * arctan(r/rc))

def vc(r):
    return sqrt(r * (gdisk(r) + ghalo(r)))

#$ This is \LaTeX : $c = 2\cdot(a+b)$
radius = range(1,30)
speed = [vc(i) for i in radius]

P.plot(radius,speed)
P.show()
