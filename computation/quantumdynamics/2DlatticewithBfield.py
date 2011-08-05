### 2D lattice with B field transport

#1D chain toys.
from pyqm import *
import os, sys
import time

### Parameters:

L = 101 #mesh size
sigma = 4.21 #gaussian standard deviation
X,Y  = create2Dmesh(1,L)
Kx,Ky  = create2DmeshK(1,L)
q = 3 #number of flux/unit cell
T = 20 #number of timesteps
t1 = 1.
t2 = 3.
kx = 1.
deltat = 1.

## Hamiltonian

def Hamiltonian(kx, ky): #returns a N by N matrix for a given kx, ky
    H1 = matrix(diag(ones(q-1, complex), 1))
    H1[q-1,0] = exp(kx*1j)
    H2 = matrix(diag([cos(2*pi*j/q + ky) for j in range(q)]))
    H = t1*H1+ t2*H2
    return trace(H+ H.getH())/q

print "hello world"

#def Hamiltonian(kx,ky,q):
    #return -cos(kx) - cos(ky)

def gaussian(x, sigma):
    return  exp(-x*x/(2.*sigma**2))/(2.*pi*sigma**2)**.5

#http://docs.scipy.org/doc/numpy/reference/routines.fft.html#background-information
#def FT(wavefn): return fftshift(abs(fft2(wavefn)))

### Evolution
def evolve1step(wavefn,H,deltat):
    '''H must be a matrix in K basis'''
    def timestep(wavefn):
        return wavefn * (1 - 1j * deltat * H)
    return renormalize(ifft2(timestep(fft2(wavefn))))



def evolve(initial, t):
    allvec = [initial]
    for i in range(t):
        allvec.append(evolve1step(allvec[i],Hamiltonian(Kx,Ky),deltat))
    return allvec

   
### Actions

R = sqrt(X**2 + Y**2)
Initial = renormalize(gaussian(R,sigma))

ion()
def plotsurface(Z):
    fig = figure()
    ax = Axes3D(fig)
    ax.set_zlim3d(0,.12)
    for i in Z:
        ax.plot_surface(X,Y,i)
        draw()


amplitudes = map(abs,evolve(Initial,T))
plotsurface(amplitudes)


#createvideo(spectrums,plotter)
#createvideo(amplitudes,plotsurface)



