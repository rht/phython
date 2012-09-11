from pylab import *
from time import time
class checktime:
    def __init__(self):
        self.tic = time()
    def time(self):
        newtic = time()
        print newtic - self.tic
        self.tic = newtic
check = checktime()

#http://www.physics.ohio-state.edu/~ntg/780/readings/hjorth-jensen_notes2011_14.pdf

def energy_local(x, alpha):
    return alpha **2 + x**2 * (1 - alpha**4)

def prob_density(x, alpha):
    return sqrt(alpha) / pi**.25 * exp(-x*x*alpha*alpha/2.)

#exact_value of energy
energy_exact = lambda x: x * x / 2. + 1. / 2 / x / x

mcs = 500  # number of Monte Carlo samplings
#here it is as if we're assuming a gaussian wavefunction, since x is normally
#distributed
def find_energy_old(alpha):
    energy_total = 0
    energy_total2 = 0
    for i in range(mcs):
        x = random()
        E_l = energy_local(x, alpha)
        energy_total += E_l
        energy_total2 += E_l * E_l
    return energy_total / mcs

def find_energy(alpha):
    x = random(mcs) / sqrt(2) / alpha
    energy_locals = energy_local(x, alpha)
    energy = sum(energy_locals) / mcs
    energy2 = sum(energy_locals ** 2) / mcs
    return energy

#metropolis
def find_energy_metropolis_old(alpha):
    x = random(mcs) / sqrt(2) / alpha
    x = [random()]
    step = .1
    for i in range(mcs-1):
        lastx = x[-1]
        newx = lastx + random() * .1
        ratio = prob_density(newx, alpha) / prob_density(lastx, alpha)
        if ratio > 1:
            x.append(newx)
        elif random() < ratio:
            x.append(newx)
    x = array(x)
    energy_locals = energy_local(x, alpha)
    energy = sum(energy_locals) / mcs
    energy2 = sum(energy_locals ** 2) / mcs
    return energy



ion()
X = arange(0.2,2.4,.1)
plot(X, map(find_energy, X))
plot(X, map(find_energy_metropolis_old, X))
plot(X, map(energy_exact, X))

check.time()

raw_input()
exit()
