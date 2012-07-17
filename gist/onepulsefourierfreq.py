from pylab import *
from scipy import fftpack

#Just want to see how the fourier frequencies are distributed inside a single square pulse
N=1000
pulse = zeros(N)
pulse[80:91] = 1.
freqs = fftpack.rfftfreq(N)

plot(pulse)
figure()
plot(freqs,fftpack.rfft(pulse))
show()
#savefig("pulsefft.png")
