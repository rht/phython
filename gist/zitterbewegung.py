from pylab import *
from scipy.linalg.matfuncs import expm2
##Pauli spin matrices

sx = matrix([[0, 1],[ 1, 0]])
sy = matrix([[0, -1j],[1j, 0]])
sz = matrix([[1, 0],[0, -1]])

##Dirac matrices
def mgamma(x):
    sigma = select([x==0,x==1,x==2,x==3], [eye(2,2),sx,sy,sz])
    return kron(sx,sigma)


#Dirac hamiltonian
#http://en.wikipedia.org/wiki/Dirac_equation
mass = 1.
p = array([1,0,0])
Ep = sqrt(sum(p*p) + mass*mass)
u_p0 = array([ 1, 0, 0, p[0]/(Ep+mass) ])
def expect(observable,ket):
    #return expectation value at state ket
    return dot(ket.conj() , dot(observable,ket).T)

def H_dirac(p):
    #return mgamma(0).I *(sum([ p[i]*mgamma(i+1) for i in range(3)]) + mass)
    return (p[0]*mgamma(1) + p[1]*mgamma(2) + p[2]*mgamma(3) + mass)


H = H_dirac(p)
HI = inv(H)
alpha1 = mgamma(0) * mgamma(1)

def evolve_x1(x0,t):
    return x0 + expect(p[0] * HI * t + .5j * HI * (alpha1 - HI) * (expm2(-2j*H*t) - 1), u_p0)
plot(*array([(t,real(evolve_x1(0,t))) for t in linspace(0,30) ]).T)
show()

