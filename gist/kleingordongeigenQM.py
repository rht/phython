#reference https://kluedo.ub.uni-kl.de/files/1171/no_series_145.pdf
#http://en.wikipedia.org/wiki/Square_root_of_a_matrix
#Dirac oscillator http://iopscience.iop.org/0305-4470/22/17/002/pdf/0305-4470_22_17_002.pdf
from pylab import *
from scipy.linalg.matfuncs import sqrtm

N = 400
s = 1 # 1/sqrt(m\omega_0)

n = arange(1,N)
m = sqrt(n)

x = matrix(s /sqrt(2) * (diag(m,-1) + diag(m,1)))
p = matrix(1j/(s*sqrt(2)) * (diag(m,-1) - diag(m,1) ))
#H = (p*p + x*x) /2.
#H = (p*p + x**4) /2.

#Klein-Gordon Hamiltonian
mass = 1000
H_KG = sqrtm( p*p + x*x + mass)  - mass
H = H_KG
def eigval_KG_analytic(n): return sqrt(mass**2 + 2 * (n + .5) * mass) - mass
def eigval_classical(n): return (n+.5)

value, vector = eig(H)
value = sort(value)[1:]
plot(value)
plot(*array([(i,eigval_KG_analytic(i)) for i in range(len(value))]).T)
plot(*array([(i,eigval_classical(i)) for i in range(len(value))]).T)
legend(['numerical','analytical','nonrelativistic'])
show()
