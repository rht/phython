###2D lattice B field brute force

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

L = 21 #mesh size
sigma = 1. #gaussian standard deviation
X = arange(-(L-1)/2, (L+1)/2)
Xmax = arange(-3/2*(L-1), 3/2*(L+1))
Y = X
K = linspace(-pi, pi, num = L)
q = 3 #number of flux/unit cell
T = 300 #number of timesteps
t1 = 3.
t2 = 3.
E_0 = -1
kx = .5

### Initial state:

Initial = array(zeros(q*L**2), complex).reshape(L,L,q)
def gaussian(x, sigma):
    y = e**(-x**2/(4.*sigma**2))/(2.*pi*sigma**2)**.5

    return y

for i in range(L):
    for j in range(L):
        for k in range(q):
            x = ln.norm(array([X[i]+q, Y[j]]))
            Initial[i,j,k] = gaussian(x, sigma)*exp(kx*(X[i]+k)*1j)

### Hamiltonian:

def Hamiltonian():
    H1 = E_0*diag(ones(q*L**2, complex))
    
    H2q = t1*diag(ones(q-1, complex), -1)
    H2xy= diag(ones(L*L, complex)).reshape(L**2,L**2)
    H2 = ln.kron(H2xy, H2q)

    
    H3q = diag([t2*exp(2*pi*k*1j/q) for k in range(q)])
    H3x = diag(ones(L, complex))
    H3y = diag(ones(L-1, complex), -1)
    H3xy = ln.kron(H3x, H3y)
    H3 = ln.kron(H3xy, H3q)

    H4q = zeros(q**2).reshape(q,q)
    H4q[0,q-1]= t1
    H4x = diag(ones(L-1, complex), -1)
    H4y = diag(ones(L, complex))
    H4xy = ln.kron(H4x, H4y)
    H4 = ln.kron(H4xy, H4q)

    H = H1+ H2 + H3+ H4

    H = H + array(matrix(H).getH())

    return H



### evolving the Hamiltonian

def timestep(vec, H):
    vec = vec + dot(H,vec)/(.05*1j)
    normedvec = vec/(ln.norm(vec))

    return normedvec



def evolve(initial, t):
    initial = initial.flatten()
    allvec = [initial]
    for i in range(t):
        allvec.append(timestep(allvec[i], Hamiltonian()))

    for i in range(t):
        allvec[i] = allvec[i].reshape(L, L, q)
 
    return allvec


### actions


wavefunctions = evolve(Initial, T)
amplitudes = [abs(vec) for vec in wavefunctions]


    

def plotsurface(function):
    x,y = meshgrid(X,Y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.plot_surface(x,y,function[:,:,0])



for i in range(T):
    plotsurface(amplitudes[i])
    filename = str('%03d' % i) + '.png'
    plt.savefig(filename, dpi=100)
    print 'Wrote file', filename
    plt.clf()

command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=25',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

os.spawnvp(os.P_WAIT, 'mencoder', command)




    





    
