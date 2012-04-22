from pyqm import *
import os, sys
import time

### Parameters:

nx = 500 #mesh size
sigma = .4 #gaussian standard deviation
dx = .01 #spatial resolution
X = create1Dmesh(dx,nx)
X = X.reshape(len(X),-1)
T = 10 #number of timesteps
dt = 1/(2/(dx**2)) #.08
#python toolbox nibot crank-nicolson  ilikear

## Hamiltonian
def potential(x): return renormalize(x**2)
Hamiltonian = -laplacematrix(dx,nx) #+ X*X

### Time Evolution
def propagate(H,dt):
    '''H must be a matrix in K basis'''
    #return (1 - 1j * dt * H)
    from scipy.linalg import expm 
    return matrix(expm(-1j*H))#* dt * H)
    #return inv((1 + .5j * dt * H)) * (1 - .5j * dt * H)

def evolve(initial, t):
    U = propagate(Hamiltonian,dt)
    allvec = [initial]
    for i in range(t):
        #allvec.append(renormalize(U*allvec[i]))
        allvec.append(U*allvec[i])
    return allvec
    #return [initial, U*initial, U*U*initial, U*U*U*initial]

   
### Actions

Initial = renormalize(gaussian(X,sigma))
#Initial = renormalize(sin(X))

ion()
amplitudes = map(abs,evolve(Initial,T))
animate1D(X,amplitudes)
#waitforbuttonpress()


#createvideo(amplitudes,plotsurface)



