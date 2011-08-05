### 2D lattice with B field transport

#1D chain toys.

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
import scipy.linalg as ln
import scipy.sparse as sp
import scipy.stats as st
import os, sys

##Pauli spin matrices

sx = matrix([[0, 1],[ 1, 0]])
sy = matrix([[0, -1j],[1j, 0]])
sz = matrix([[1, 0],[0, -1]])

smin = matrix([[0,1],[0,0]])
smax = matrix([[0,0],[1,0]])

### Parameters:

L = 101 #mesh size
sigma = .05 #gaussian standard deviation
X = arange(-(L-1)/2, (L+1)/2)
Xmax = arange(-3/2*(L-1), 3/2*(L+1))
Y = X
K = linspace(-pi, pi, num = L)
q = 3 #number of flux/unit cell
T = 20 #number of timesteps
t1 = 1.
t2 = 1.



## Hamiltonian

def Hamiltonian(kx, ky, N): #returns a N by N matrix for a given kx, ky
    H1 = matrix(diag(ones(N-1, complex), 1))
    H1[N-1,0] = exp(kx*1j)
    H2 = matrix(diag([cos(2*pi*j/N + ky) for j in range(N)]))
    H = t1*H1+ t2*H2
    return H+ H.getH()


### Initial state:

Initial = array(zeros(q*L**2)).reshape(L,L,q)
def gaussian(x, sigma):
    y = e**(-x**2/4.*sigma**2)/(2.*pi*sigma**2)**.5

    return y

for i in range(L):
    for j in range(L):
        for k in range(q):
            x = ln.norm(array([X[i]+q, Y[j]]))
            Initial[i,j,k] = gaussian(x, sigma)


### Fourier transforms:

def FT(wavefn): #return the normalized array of fourier transformed wavefunction
    FT = array(zeros(q*L**2), complex).reshape(L,L,q)
    for k in range(q):
        FT[:,:,k]= fft.fft2(wavefn[:,:,k])

    return FT

def IFT(wavefn):

    IFT = array(zeros(q*L**2), complex).reshape(L,L,q)
    for k in range(q):
        IFT[:,:,k]= fft.ifft2(wavefn[:,:,k])

    IFT = IFT/ln.norm(IFT)

    return IFT
    

### Evolution

def timestep(wavefn, H):
    nextwavefn = array(zeros(q*L**2), complex).reshape(L,L,q)

    for i in range(L):
        for j in range(L):
            vec = wavefn[i,j,:]
            vec = vec + dot(H(K[i],K[j],q),vec)/(5.*1j)
            nextwavefn[i,j,:] = vec
    
    nextwavefn= nextwavefn/ln.norm(nextwavefn)
    
    return nextwavefn

    
def evolve(initial, t):
    allvec = [initial]
    FTallvec = [FT(initial)]
    for i in range(t):
        FTnextwavefn = timestep(allvec[i], Hamiltonian)
        nextwavefn = IFT(FTnextwavefn)
        allvec.append(nextwavefn)
        FTallvec.append(FTnextwavefn)

    return FTallvec

### Actions

def plotsurface(function):
    x,y = meshgrid(X,Y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x,y,function[:,:,0])
    plt.show()

plotsurface(abs(FT(Initial)))




##wavefunctions = evolve(Initial, T)
##amplitudes = [abs(vec) for vec in wavefunctions]
##
##for i in range(len(amplitudes)):
##    plotsurface(amplitudes[i])
##    filename = str('%03d' % i) + '.png'
##    plt.savefig(filename, dpi=100)
##    print 'Wrote file', filename
##    plt.clf()

##command = ('mencoder',
##           'mf://*.png',
##           '-mf',
##           'type=png:w=800:h=600:fps=25',
##           '-ovc',
##           'lavc',
##           '-lavcopts',
##           'vcodec=mpeg4',
##           '-oac',
##           'copy',
##           '-o',
##           'output.avi')
##
##os.spawnvp(os.P_WAIT, 'mencoder', command)






