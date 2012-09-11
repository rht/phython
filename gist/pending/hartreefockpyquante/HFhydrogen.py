from PyQuante.Molecule import Molecule
from PyQuante import SCF

#Intializing the molecules
H2 = Molecule('H2', [ (1, (0,0,0)), (1,(1.4,0,0)) ] )
OHmin = Molecule( 'OHmin', [ (8, (0., 0., -0.08687037)) , (1, (0., 0., 0.86464814)) ], units='Angstrom', charge=-1  )


#HF method example
#Different basis can be used, sto-3g
#H2solver = SCF(H2, method="HF")
H2solver = SCF(H2, method="HF")
H2solver.iterate()
print "HF energy = ", H2solver.energy

#DFT method example
lda = SCF(H2,method="DFT")
lda.iterate()
blyp = SCF(H2,method="DFT",functional="BLYP")
blyp.iterate()
print "DFT results: LDA = ", lda.energy, " BLYP = ", blyp.energy

#Open shell Hartree Fock
# method='UHF' or method='ROHF'
