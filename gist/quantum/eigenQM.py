#reference https://kluedo.ub.uni-kl.de/files/1171/no_series_145.pdf
#"We observe a fast convergence towards the exact eigenvalues E_n for N > E_n +5"
from pylab import *

N = 4
s = 1 # 1/sqrt(m\omega_0)

n = arange(1,N)
m = sqrt(n)

x = matrix(s /sqrt(2) * (diag(m,-1) + diag(m,1)))
p = matrix(1j/(s*sqrt(2)) * (diag(m,-1) - diag(m,1) ))
#HO
#H = (p*p + x*x) /2.
#heliumtable
H = p*p/2. + where(x < 0, - 1000, 3) * x
#H = (p*p + x**4) /2.
#print H
value, vector = eig(H)
print real(min(value))
