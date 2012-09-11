from pylab import *

#utilities
#1. timing



#physics
weyl = matrix([[0, 1],[ -1, 0]])
sx = matrix([[0, 1],[ 1, 0]])
sy = matrix([[0, -1j],[1j, 0]])
sz = matrix([[1, 0],[0, -1]])
pauli = lambda x: [sx, sy, sz][x - 1]

# Weyl basis
index = [0, 1, 2, 3] 
gamma0 = kron(sx, eye(2))
gamma5 = -kron(sz, eye(2))
gamma = lambda x: gamma0 if x == 0 else kron(weyl, pauli(x))


def fourdot(p, k):
    # four vector dot product
    return p[0] * k[0] - dot(p[1:], k[1:])


def fourvector(threevector, m):
    # generate on-shell fourvector only
    return array([sqrt(m ** 2 + norm(threevector) ** 2)] + list(threevector))


def fermion_u(p, spin, m):
    # Peskin 46
    #http://en.wikipedia.org/wiki/Dirac_spinor
    E, p1, p2, p3 = p
    if spin == 1:
        return sqrt(E + m) * array([1, 0, p3 / (E + m),
            (p1 + 1j * p2) / (E + m)])
    else:
        return sqrt(E + m) * array([0,
            1, p1 - 1j * p2 / (E + m),
            -p3 / (E + m)])

def slash(k):
    from numpy import add
    return matrix(add.reduce([gamma(i) * k[i] for i in
        range(4)]))


def fermion_ubar(p, spin, m):
    return fermion_u(p, spin, m).T


def antifermion_v(p, spin, m):
    # Peskin 46
    #http://en.wikipedia.org/wiki/Dirac_spinor
    E, p1, p2, p3 = p
    if spin == 1:
        return sqrt(E + m) * array([p1 - 1j * p2 / (E + m), -p3 / (E + m), 0, 1])
    else:
        return sqrt(E + m) * array([p3 / (E + m), (p1 + 1j * p2) / (E + m), 1,
            0])


def antifermion_vbar(p, spin, m):
    return antifermion_v(p, spin, m).T


