from pylab import *
from scipy.integrate import odeint
from numpy.fft import rfft
import os
from pyqm import createvideo, centralderivative


ion()
ioff()
# initialization
N = 32
L = 1.
dx = L / N
Ncycle = 10
#initial phase space coordinates
Xs = [sin(pi * linspace(0,L,num=N))]
dt = sqrt(.125)
Time = linspace(0, Ncycle * dt, num=Ncycle)
Cycles = Time * 2 * sin(pi / N / 2)
alpha = 2.25
beta = 8


def evolve(X, t=0):
    D = centralderivative(dx, N)
    v = -alpha * dot(D, X) - 1. / 24 * dot(D * D * D, X)
    return v


# the meat of the code is in this very single line
Xs = odeint(evolve, Xs[0], Time)



figures = []
xsmin, xsmax = Xs.min(), Xs.max()

# it turns out python can't accomodate a list bigger than 3000
for i in range(Ncycle):
    fig = figure()
    ylim(xsmin, xsmax)
    plot(Xs[i])
    print "video %d" %i
    figures.append(fig)
createvideo(figures,prefix=0)

