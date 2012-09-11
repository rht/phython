from pylab import *
import time
ion()
#failed attempts
#1. values, vectors = eigsh(laplacematrix(dx, nx), which='SA')
#2. http://sfepy.org/doc-devel/_modules/sfepy/linalg/eigen.html

#remarks
#1. the numerov/shooting method is only suitable for even-parity potential http://www.hep.vanderbilt.edu/~maguirc/physics245Fall06/p245_lect26.pdf
 #-> not necessarily true, because we can start from far negative x axis

# solving 1D schrodinger equation numerically
#http://en.wikipedia.org/wiki/Numerov's_method
#http://iopscience.iop.org/0022-3700/15/6/009/pdf/0022-3700_15_6_009.pdf
#http://physics.bu.edu/~py502/lectures4/schrod.pdf
#http://k2.chem.uh.edu/quantum/Supplement/Numerov/
# for time dependent http://en.wikipedia.org/wiki/Pseudo-spectral_method

#improve the way we pick up E

tic = time.time()
# variables
k = .01
#m = 2000
m = 1
xmin = -3
xmax = 3.
dx = 0.01
omega = sqrt(k/m)
nsteps = int((xmax-xmin)/dx)
eps = omega / 100.  # tolerance
print "Energy step", eps
X = arange(xmin,xmax,dx)


def V(x):
    #return exp(-x*x)
    return .5*k*x*x
    #return where(x < .1, -1000, 3) * x
    #return .5 * k * x * x * x * x


def G(x, energy, V):
    #almost inv green function
    #(-d^2 +(V-E))psi = 0
    return 2 * m * (energy - V(x))


def generatepsi(E, V):
    delta_nodes = 0
    psi = [.0, .0001]
    nodes = 0

    for n in range(1,nsteps-1):
        # the meat of numerov method
        numerator = (2 * psi[n] - psi[n-1] -
                dx*dx/12.*(10*G(X[n],E, V)*psi[n]+G(X[n-1], E, V)*psi[n-1]))
        denominator = 1 + dx * dx / 12. * G(X[n+1], E, V)
        psip = numerator / denominator
        psi.append(psip)
        if psip * psi[-2] < 0:
            nodes += 1
    return psi, nodes


def find_eigenfunction(Emin, Emax, V):
    #choosing energy using bisection method
    tolerance = 1e-4
    MAX_ITER = 50
    E = (Emin + Emax) / 2.
    psi, nodes = generatepsi(E, V)
    i = 0
    while (abs(psi[-1]) > tolerance) and (i <= MAX_ITER):
        E = (Emin + Emax) / 2.
        psi, nodes = generatepsi(E, V)
        if psi[-1] * generatepsi(Emax, V)[-1] < 0:
            Emin = E
        else:
            Emax = E
        i += 1
    print "Energy", E/omega, "after", i, "iterations"
    print "nodes", nodes
    # renormalize psi here
    psi /= norm(psi)
    return psi, E




#Energy = 1 * omega / 2.
energy_intervals = arange(10)
for i in range(len(energy_intervals) - 1):
    psi, E  = find_eigenfunction(energy_intervals[i],energy_intervals[i+1], V)
    figure()
    plot(X, psi)
print time.time() - tic

raw_input()
exit()


# further exploration
# http://www.mathworks.ir/downloads/Spectral%20Methods%20in%20MATLAB%5Bwww.mathworks.ir%5D.pdf


# Several methods to compute the eigenfunctions
#1. lanczos algorithm
#2. chebyshev spectral method
#3. 
