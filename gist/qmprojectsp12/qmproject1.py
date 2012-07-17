from pylab import *
from scipy.linalg.matfuncs import expm2

def propagator(H,t):
    return expm2(-1j*H*t) #expm2: compute the matrix exponential using eigenvalue decomposition

def qitpropagator(H,t):
    #alternative implementation of propagator from pyqit
    d, v = eig(H)
    return dot(dot(v,diag(exp(-1j*t*d))), v.conj().transpose())

def expect(observable,ket):
    #return expectation value at state ket
    return dot(ket.conj().transpose(), dot(observable,ket))

def evolve(initial,H,tfinal,timestep):
    time = arange(0,tfinal,timestep)
    states = [ket0]
    U = propagator(H,timestep)
    #evolve the system
    for i in range(len(time)-1):
        states.append(array(dot(U,states[-1])))
    return time, states

def E_spectrum_plotter(H,param):
    eigvalues = array([eigvals(H(i)) for i in param])
    for i in eigvalues.T:
        plot(param,i)


# Week 1, two state system
#initialize a two-state system at |0>
ket0 = array([1.,0])
#Hamiltonian = matrix([[0, 1], [1, 0]])
Hamiltonian = lambda x: matrix([[4+x, 1], [1, 4-x]])

#evolve
#time, states = evolve(ket0,Hamiltonian,10,.1)
#groundstate = map(real,zip(*states)[0])
#excitedstate = map(imag,zip(*states)[1])
#plot(time,groundstate,'blue')
#plot(time,excitedstate,'red')


#spectrum
#feynman ch9
E_spectrum_plotter(Hamiltonian,arange(0,10,.5))





show()


# week 3
#perturbation theory stuff
def energy1(H0,H1,ket):
    return expect(H1,ket)

def ket1(H0,H1,Espectrum,ketsk0,ketn0,ketn0i):
    return sum([dot(ketsk0[i].conj().transpose(),dot(H1,ketn0)) / (Espectrum[ketn0i] - Espectrum[i]) for i in range(len(kets0)) if ketsk0[i] != ketn0])

