
from matplotlib.pyplot import plot, show
from numpy import *

epsilon = 0.4
delta = 0.1
M = 20

def V(x):
    return 0.5 * pow(x,2)

def rhosmall(xb,xa):
    C = 1#m/(2*pi*pow(hbar,2)*epsilon)
    return sqrt(C) * exp(-C*pow((xb - xa),2) - epsilon/2 * (V(xb) + V(xa)))

def rho(id,jd,N):
    if N==0:
        return rhosmall(id,jd)
    else:
        return delta * sum( [rho(id,k*delta,N-1) * \
                rho(k*delta,jd,N-1) for k in range(-M,M+1)])

ntotal = 2
Z = delta / 2 * sum( [rho((k+1)*delta,(k+1)*delta,ntotal) + rho(k*delta,k*delta,ntotal) for k in range(-M,M+1)])

def wavefunction(x):
    return rho(x,x,ntotal) / Z

xarray = arange(-2,2,delta)
wavearray = [wavefunction(i) for i in xarray]

plot(xarray,wavearray)
show()
