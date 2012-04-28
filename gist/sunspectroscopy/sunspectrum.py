from pylab import *
# http://rredc.nrel.gov/solar/spectra/am1.5/

wavelength, etr, globaltilt, direct = loadtxt("ASTMG173.csv", skiprows=2,
        unpack=1, delimiter=',')

plot(wavelength, etr)
savefig("sunspectrum.png")
#show()
