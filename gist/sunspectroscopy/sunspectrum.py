from pylab import *
# http://rredc.nrel.gov/solar/spectra/am1.5/
# http://www.krioma.net/pyspec.php

wavelength, etr, globaltilt, direct = loadtxt("ASTMG173.csv", skiprows=2,
        unpack=1, delimiter=',')

print(len(etr))
plot(wavelength[30:], etr[30:])
xlabel("wavelength(nm)")
#savefig("sunspectrum.png")
show()
