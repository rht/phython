from pylab import *

K = arange(-100,100,.01)
def psi(x):
    return sum([exp(1j * x * k) for k in K])

def phi(x):
    return sum([exp(1j * x * k) / (k * k + 10) if k >.2 else 0 for k in K])


X = arange(-10,10,.1)
psis = [psi(x) for x in X]
phis = [phi(x) for x in X]
plot(X, psis / norm(psis))
plot(X, phis / norm(phis))
savefig("qftfree.png")
#show()
