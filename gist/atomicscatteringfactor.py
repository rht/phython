from pylab import *
import urllib

#http://henke.lbl.gov/optical_constants/asf.html
#i.e. atomic form factor

a = urllib.request.urlopen("http://henke.lbl.gov/optical_constants/sf/hg.nff")

E, f1, f2 = loadtxt(a, skiprows=1)[100:].T
plot(E,f1)
plot(E,f2)
legend(["f1","f2"])
xlabel("E(eV)")
show()
