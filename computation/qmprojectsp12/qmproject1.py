from pylab import *
from scipy.linalg.matfuncs import expm2

def propagator(H,t):
    return expm2(-1j*H*t) #expm2: compute the matrix exponential using eigenvalue decomposition

def qitpropagator(H,t):
    #only works for time-independent H
    #alternative implementation from pyqit
    #convince yourself why this code works
    d, v = eig(H)
    return dot(dot(v,diag(exp(-1j*t*d))), v.conj().transpose())

#initialize a two-state system at |0>
ket0 = array([1.,0])
Hamiltonian = matrix([[0, 1], [1, 0]])
timestep = .1
U = propagator(Hamiltonian,timestep) #time evolution operator
states = [ket0]

#time coordinates
time = arange(0,10,timestep)

for i in range(len(time)-1):
    states.append(array(dot(U,states[-1])))

#list of the states
groundstate = map(real,zip(*states)[0])
excitedstate = map(imag,zip(*states)[1])

plot(time,groundstate,'blue')
plot(time,excitedstate,'red')
show()
