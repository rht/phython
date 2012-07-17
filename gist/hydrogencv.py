#http://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1#Thermo-Gas

#selected properties on hydrogen
#http://www.boulder.nist.gov/div838/Hydrogen/PDFs/McCarty.Hord.Weber.1981.monograph168.pdf
#selected topics on hydrogen fuel
#http://www.boulder.nist.gov/div838/Hydrogen/PDFs/Hord(ed.)1975.NBSIR%2075-803.pdf
#polarizability of the hydrogen molecule
#http://jcp.aip.org/resource/1/jcpsa6/v46/i4/p1426_s1

#low temperature specific heat (ortho vs para)
#http://ocw.tudelft.nl/fileadmin/ocw/courses/AdvancedStatisticalMechanics/res00030/!486973746f7279206f6620746865206465736372697074696f6e206f662074686520687964726f67656e20676173.pdf
#ionization of hydrogen
#http://csep10.phys.utk.edu/astr162/lect/stars/spectra.html

#less useful links
#http://www.wag.caltech.edu/home/jsu/Thesis/node31.html

from pylab import *
from scipy.constants import *
from scipy.misc import derivative as deriv
from scipy.integrate import quad
#from PyQuante.Molecule import Molecule
#from PyQuante.hartree_fock import rhf

#brute force anharmonic specific heat calculation with partition function
def anharmoniccvib(T):
    k2 = 0.782
    k3 = -34.1
    k4 = 34.88
    def logZ(T):
        foo = lambda x: exp(( -k2*x*x -k3*x*x*x - k4*(x**4) ) /  T )
        return log(quad(foo,-Inf,Inf)[0])
    energy = lambda x: x**2 * deriv(logZ,x,dx=1e-6)
    #return R*( 1/2. + deriv(energy,T,dx=1e-6 )) # 1/2 from kinetic energy
    #from perturbation theory
    return R * (1 - (k*T/1.6e-19) * 9. / 4 * k4 / (k2**2))

#def bondpotential(x):
    #h2 = Molecule('h2',[(1,(0,0,0)),(1,(x,0,0))])
    #return rhf(h2)[0]

def NISTcp(T):
    param = []
    if 298 < T < 1001:
        param += [33.066178, -11.363417, 11.432816, -2.772874, -0.158558]
    elif 1000 < T < 2501:
        param += [ 18.563083, 12.257357, -2.859786, 0.268238, 1.977990  ]
    elif 2500 < T < 6000:
        param += [ 43.413560, -4.293079, 1.272428, -0.096876, -20.533862 ]
    t = T/1000.
    return param[0] + param[1]*t + param[2]*t*t + param[3]*t*t*t + param[4]/(t*t)

def atomfraction(T):
    l = 1.47 * 5e-11
    sigma = l * l
    coeff = .5 * 1e5 /(k*298) * l * sigma
    if T > 950:
        b = coeff * exp( 52452.3639093575 /T) #4.52eV/k_B
        return (sqrt(4*b + 1) -1)/(2*b)
    else: return 0

cHO = lambda x: R * (x**2) * exp(x)/((exp(x)-1)**2)

def theorycpatom(T):
    # angular frequency from polarizability 6.9e16 Hz
    ratiovibatom = 527037/T #hbar * 6.9e16 / k
    return 5.*R/2 #+ cHO(ratiovibatom)

def crot(ratiorot):
    A = sum([(2*l+1)* ((l*(l+1))**2) * exp(-l*(l+1)*ratiorot) for l in range(30) ])
    B = sum([(2*l+1)* (l*(l+1)) * exp(-l*(l+1)*ratiorot) for l in range(30) ])
    Z = sum([(2*l+1)*  exp(-l*(l+1)*ratiorot) for l in range(30) ])
    return (A/Z- (B/Z)**2) * ratiorot * ratiorot * R

def theorycpmolecule(T):
    ratiorot = 85./T
    ratiorotelectron = 1000./T
    ratiovib = 6215. /T
    ratiodipole = 26000. /T
    cvib = cHO(ratiovib)
    cdipole = cHO(ratiodipole)
    return R + 1.5*R + crot(ratiorot) + 2* cdipole + cvib  #+ anharmoniccp(T)  #+cvib

def theorycp(T):
    f = atomfraction(T)
    return (1-f)*theorycpmolecule(T) + f*theorycpatom(T)

xaxis = range(300,6000,20)
p1= plot(xaxis, [NISTcp(i) for i in xaxis])
#p2= plot(xaxis, [theorycp(i) for i in xaxis])
p2 = plot(xaxis,[theorycpmolecule(i) for i in xaxis])
xlabel("Temperature (K)")
ylabel("Specific heat (J/(kg mol))")
title("Hydrogen molecule specific heat with Ctrans,Crot,Cvib,Cdip")
#legend([p1, p2], ["exp", "theory"], loc = "lower right")
show()
