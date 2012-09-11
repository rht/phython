from pylab import *
ion()
#failed attempts
#1. values, vectors = eigsh(laplacematrix(dx, nx), which='SA')
#2. http://sfepy.org/doc-devel/_modules/sfepy/linalg/eigen.html

# solving 1D schrodinger equation numerically
#http://en.wikipedia.org/wiki/Numerov's_method
#http://iopscience.iop.org/0022-3700/15/6/009/pdf/0022-3700_15_6_009.pdf
#http://physics.bu.edu/~py502/lectures4/schrod.pdf
#http://k2.chem.uh.edu/quantum/Supplement/Numerov/
# for time dependent http://en.wikipedia.org/wiki/Pseudo-spectral_method

#improve the way we pick up E

# variables
k = .01
m = 2000
xmin = -3
xmax = 3.
dx = 0.01
omega = sqrt(k/m)
nsteps = int((xmax-xmin)/dx)
#eps = 0.003
eps = omega / 100.  # tolerance
print "Energy step", eps
X = arange(xmin,xmax,dx)


def V(x):
    #return exp(-x*x)
    return 1/2.*k*x*x


def G(x, energy, V):
    #almost inv green function
    #(-d^2 +(V-E))psi = 0
    return 2 * m * (energy - V(x))


def generatepsi(E, V):
    delta_nodes = 0
    psi = [.0, .0001]

    for n in range(1,nsteps-1):
        # the meat of numerov method
        numerator = (2 * psi[n] - psi[n-1] -
                dx*dx/12.*(10*G(X[n],E, V)*psi[n]+G(X[n-1], E, V)*psi[n-1]))
        denominator = 1 + dx * dx / 12. * G(X[n+1], E, V)
        psip = numerator / denominator
        psi.append(psip)
        # check sign flip
        if psip*psi[n] < 0:
            delta_nodes += 1
    return psi, delta_nodes


def find_eigenfunction(Emin, Emax, V, quantumnumber):
    #old algorithm
    Energy = Emin
    nodes = 0
    count = 1
    psis = []
    energies = []
    while not ((nodes==quantumnumber+1) or (Energy > Emax)):
        count += 1
        psi, delta_nodes = generatepsi(Energy, V)
        nodes += delta_nodes
        psi /= norm(psi)
        psis.append(psi)
        Energy += eps
        energies.append(Energy)
    print count
    return psis, energies


def find_eigenfunction(Emin, Emax, V, quantumnumber):
    #choosing energy using bisection method
    tolerance = 1e-4
    while True:
        E = (Emin + Emax) / 2.
    energies = []
    while not ((nodes==quantumnumber+1) or (Energy > Emax)):
        count += 1
        psi, delta_nodes = generatepsi(Energy, V)
        nodes += delta_nodes
        psi /= norm(psi)
        psis.append(psi)
        Energy += eps
        energies.append(Energy)
    print count
    return psis, energies




for n in range(1,3,2):
    quantumnumber = n
    Energy = quantumnumber * omega / 2.
    print Energy
    #psi, E  = find_eigenfunction(Energy, V, quantumnumber)
    psis, energies  = find_eigenfunction(Energy-50*eps, Energy+50*eps, V, quantumnumber)
    for i in range(len(psis)):
        figure()
        plot(X, psis[i])
        title("%f" %(energies[i] / omega))
    #figure()
    #plot(X, psi)

raw_input()
exit()


# further exploration
# http://www.mathworks.ir/downloads/Spectral%20Methods%20in%20MATLAB%5Bwww.mathworks.ir%5D.pdf


# Several methods to compute the eigenfunctions
#1. lanczos algorithm
#2. chebyshev spectral method
#3. 
