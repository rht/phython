#!/usr/bin/env python
""" Hartree-Fock program for H2 molecule. """
from scipy import *
import symeig
from pylab import *
import copy

def F0(x):
    "Special integral 1/sqrt(x)*Integrate[exp^{-y^2}, {y,0,sqrt(x)}]"
    if abs(x)<1e-12: return 1.0
    sx = sqrt(x)
    return sqrt(pi)*special.erf(sx)/(2*sx)

def CmpOverlapHam0(base, Rpos, R0):
    """ Calculates Hamiltonian K(a,b) = 3*a*b/(a+b)*(Pi/(a+b))^{3/2}
                               Vext(a,b) = -4*Pi/(a+b)  """
    Olap = zeros((len(base),len(base)), dtype=float)
    Ham = zeros((len(base),len(base)), dtype=float)
    for i,ai in enumerate(base):
        for j,bj in enumerate(base):
            R_pq = 0.5*R0*((2*Rpos[i]-1)*ai+(2*Rpos[j]-1)*bj)/(ai+bj) # center of mass
            a_ij = ai*bj/(ai+bj)

            Olap[i,j] = (pi/(ai+bj))**(3/2.)
            if Rpos[i]!=Rpos[j]:
                Olap[i,j] *= exp(-R0**2*a_ij)
            
            Ham[i,j] = 3*Olap[i,j]*a_ij # first part of kinetic term
            if Rpos[i]!=Rpos[j]:
                Ham[i,j] *= (1-2/3.*a_ij*R0**2)

            F1 = F0((ai+bj)*(R_pq+R0/2.)**2) # F0((a+b)*(R1-R_pq)^2)
            F2 = F0((ai+bj)*(R_pq-R0/2.)**2) # F0((a+b)*(R2-R_pq)^2)
            
            Ham[i,j] += -2*Olap[i,j]*sqrt((ai+bj)/pi)*(F1+F2)
    return (Olap, Ham)
    
def CmpU(base, Rpos, Olap, R0):
    " Coulomb matrix elements U_{abcd} for H2 "
    Uc = zeros((len(base),len(base),len(base),len(base)), dtype=float)
    for i,ai in enumerate(base):
        for j,bj in enumerate(base):
            for k,ck in enumerate(base):
                R_ik = 0.5*R0*((2*Rpos[i]-1)*ai+(2*Rpos[k]-1)*ck)/(ai+ck)
                prf_Olap_ik = 2*Olap[i,k]/sqrt(pi)
                b_ik = ai + ck
                for l,dl in enumerate(base):
                    
                    R_jl = 0.5*R0*((2*Rpos[j]-1)*bj+(2*Rpos[l]-1)*dl)/(bj+dl)
                    b_jl = bj + dl
                    a_ijkl = b_ik*b_jl/(b_ik+b_jl)
                    
                    Uc[i,j,k,l] = prf_Olap_ik*Olap[j,l]*sqrt(a_ijkl)*F0(a_ijkl*(R_ik-R_jl)**2)

    return Uc
    

def CmpHeffHartree(Ham0, Uc, rho):
    Heff = copy.deepcopy(Ham0)
    for i,ai in enumerate(base):      
        for j, bj in enumerate(base):
            dsum=0
            for k,ck in enumerate(base):      
                for l, dl in enumerate(base): 
                    dsum += (2*Uc[k,i,l,j]-Uc[k,i,j,l])*rho[k,l]                    
            Heff[i,j] += dsum
            
    return Heff

def Eigensystem(Heff, Olap):
    """ General eigenvalue problem solved"""
    w,Z = symeig.symeig(Heff, Olap, type=1) # symmetric generalized eigenvalue problem
    Energy = w[:2]
    return (Energy, Z[:,:2])


def CmpDensity(Z, Nell):
    """ Computes the density matrix """
    nrho = zeros((len(Z),len(Z)),dtype=float)
    for a in range(Nell):
        for i in range(len(Z)):
            for j in range(len(Z)):
                nrho[i,j] += 0.5*Z[i,a/2]*Z[j,a/2]
    return nrho

def TotEnergy(E0, rho, Ham0, Nell):
    Et=0
    for a in range(Nell):  Et +=  0.5*E0[a/2]
    return Et + sum(rho.transpose()*Ham0)

def NucleiRepulsion(R0):
    return 1*1./R0


def Solution(R0, rho, small=1e-6):
    
    (Olap, Ham0) = CmpOverlapHam0(base, Rpos, R0)
    Uc = CmpU(base, Rpos, Olap, R0)

    Et = 0
    for i in range(nmax):
        Heff = CmpHeffHartree(Ham0, Uc, rho)
        (E0, Z) = Eigensystem(Heff, Olap)
        nrho = CmpDensity(Z, Nell)
        nEt = TotEnergy(E0, rho, Ham0, Nell) + NucleiRepulsion(R0)
        rho = nrho
        if abs(nEt-Et)<small: break
        Et = nEt
    return (Et, rho)

def Root(R0, rho, small=1e-6):
    return Solution(R0, rho, small)[0]


def RadialDensity(r, z, rho, base, Rpos, R0):
    dsum=0
    for i,ai in enumerate(base):
        R_ii = 0.5*R0*(2*Rpos[i]-1)
        for j,bj in enumerate(base):
            R_jj = 0.5*R0*(2*Rpos[j]-1)
            dsum += rho[i][j]*exp(-ai*(r**2+(z-R_ii)**2)-bj*(r**2+(z-R_jj)**2))
    return dsum


if __name__ == '__main__':

    base = [13.00773, 1.962079, 0.444529, 0.1219492]*2
    Rpos = [0,0,0,0,1,1,1,1]
    Nell = 2
    nmax = 100
    
    R = 1.388

    rho = zeros((len(base),len(base)), dtype=float)
    
    (Et, rho) = Solution(R, rho, small=1e-6)


    print 'Computing E(R_HH) to plot energy profile'
    r = arange(0.5, 3., 0.1)
    Er = [Solution(x, rho, small=1e-6)[0] for x in r]
    plot(r, Er, 'o-')
    xlabel('R[a.u.]')
    ylabel('Energy[Hartree]')
    show()
    
    print 'Minimizing E(R_HH) to find precise minimum'
    Rmin = optimize.fmin_powell(Root, 1., args=(rho,))    
    print 'Rmin=', Rmin
    (Et, rho) = Solution(Rmin, rho, small=1e-6)


    print 'Ploting radial density distribution'
    wr = arange(-1.5,1.5,0.05)
    wz  =arange(-1.5,1.5,0.05)
    dens=zeros((len(wr),len(wz)), dtype=float)
    for ir,r in enumerate(wr):
        for iz,z in enumerate(wz):
            dens[ir,iz] = RadialDensity(r, z, rho, base, Rpos, Rmin)

    X, Y = meshgrid(wr, wz)
    contour(X, Y, dens)
    xlabel('r')
    ylabel('z')
    show()

    
    
