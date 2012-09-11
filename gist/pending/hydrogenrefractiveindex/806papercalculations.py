from scipy.constants import *

a0 = (value("Bohr radius") * 100 ) # in cm

#static hydrogen molecule polarizability

alphaexp = 0.8042e-24 # in cm3

print alphaexp / (a0 ** 3)
