from pylab import *
from itertools import product
#nested loops with product
#http://docs.python.org/library/itertools.html#itertools.product
from functools import reduce

weyl = matrix([[0, 1],[ -1, 0]])
sx = matrix([[0, 1],[ 1, 0]])
sy = matrix([[0, -1j],[1j, 0]])
sz = matrix([[1, 0],[0, -1]])
def pauli(x):
    return [sx, sy, sz][x - 1]


#gamma0 = kron(sz, eye(2))
#gamma1 = kron(sy, eye(2))

# Weyl basis
index = [0, 1, 2, 3] 
gamma0 = kron(sx, eye(2))
gamma5 = -kron(sz, eye(2))
mgamma = lambda x: gamma0 if x == 0 else kron(weyl, pauli(x))

tgamma = [mgamma(i) for i in index]
tgamma5 = [gamma5 for i in index]
def binaryproduct(a, b):
    #return array([a[i,] * b for i in index])
    return [reduce(dot, i) for i in product(a,b)]

def traceofgammas(*args):
    return reshape([trace(reduce(dot, i)) for i in product(*args)],
            len(args) * [4])

#possibly useful
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.einsum.html

#print(shape(tproducts(tgamma,tgamma,tgamma)))
#print(tproducts(tgamma,tgamma))
print(traceofgammas(tgamma,tgamma,tgamma))
