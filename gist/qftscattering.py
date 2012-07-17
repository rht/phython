from pylab import *
from scipy.stats import norm
from scipy.linalg.matfuncs import expm2

def phi(x):
    k = 5
    return norm.pdf(x-4) * exp(1j * k * x) + norm.pdf(x-40) * exp(-1j * k * x)

def evolve(initial,H,tfinal,timestep):
    def propagator(H,t):
        return expm2(-1j*H*t) #expm2: compute the matrix exponential using eigenvalue decomposition
    time = arange(0,tfinal,timestep)
    states = [initial]
    U = propagator(H,timestep)
    #evolve the system
    for i in range(len(time)-1):
        states.append(array(dot(U,states[-1])))
    return time, states


def laplacematrix(dx,nx):
    return matrix((diag(ones(nx-1),1) + diag(ones(nx-1),-1) + diag(-2*ones(nx))) / (dx*dx))


X = arange(0,50,.1)

H = -laplacematrix(.1, len(X))

time, states = evolve(phi(X), H, 5, 1)

for i in states:
    figure()
    plot(X,i)

show()
