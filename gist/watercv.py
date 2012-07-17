#how to get data
#http://webbook.nist.gov/cgi/fluid.cgi?P=0.1&TLow=273.16&THigh=1275.0&TInc=1&Digits=5&ID=C7732185&Action=Load&Type=IsoBar&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm&RefState=DEF

#direct link to data
#http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C7732185&Type=IsoBar&Digits=5&P=0.1&THigh=1275.0&TLow=273.16&TInc=1&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm

#http://en.wikipedia.org/wiki/Water_(data_page)

#SEE THIS
#http://www.lsbu.ac.uk/water/vibrat.html

from pylab import *
from scipy.constants import *
from scipy.misc import derivative as deriv
from scipy.integrate import quad
#from PyQuante.Molecule import Molecule
#from PyQuante.hartree_fock import rhf

def NISTcp(T):
    #http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=2#Thermo-Condensed
    param = [-203.6060, 1523.290, -3196.413, 2474.455, 3.855326]
    t = T/1000.
    cp = param[0] + param[1]*t + param[2]*t*t + param[3]*t*t*t + param[4]/(t*t)
    return cp /1000 * 55.55


waterdata = loadtxt("waterdata.txt", unpack=1, usecols=arange(10), skiprows=2)
T = waterdata[0]
Cv = waterdata[7]
Cp = waterdata[8]

#plot(T,Cp)
plot(*zip(*[(i,NISTcp(i)) for i in arange(300, 500, 0.1)]))
show()
exit()









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

#xaxis = range(300,6000,20)
#p1= plot(xaxis, [NISTcp(i) for i in xaxis])
##p2= plot(xaxis, [theorycp(i) for i in xaxis])
#p2 = plot(xaxis,[theorycpmolecule(i) for i in xaxis])
#xlabel("Temperature (K)")
#ylabel("Specific heat (J/(kg mol))")
#title("Hydrogen molecule specific heat with Ctrans,Crot,Cvib,Cdip")
##legend([p1, p2], ["exp", "theory"], loc = "lower right")
#show()
