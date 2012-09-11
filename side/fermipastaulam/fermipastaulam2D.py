from pylab import *
from scipy.integrate import odeint
from numpy.fft import rfft
import os
# "The only computer experiments worth doing are those that yield a surprise"
# "Metropolis And von Neumann Install Awful Computer"
# Our problem turned out to have been felicitously chosen. The results were
# entirely different qualitatively from what even Fermi, with his great knowledge of wave motions, had expected.

#Computational physics was born! http://www-dft.ts.infn.it/~resta/tcse/lec1.pdf
# For more than two decades, most of genuine computational physics adressed classical statistical mechanics (or classical dynamical systems).
# Modern quantum mechanics of materials (alias the birth of computational electronic structure)


ion()
# initialization
N = 32
L = 1.
Ncycle = 30000
#initial phase space coordinates
XVs = [concatenate((sin(pi * linspace(0,L,num=N)), zeros(N)))]
dt = sqrt(.125)
Time = linspace(0, Ncycle * dt, num=Ncycle)
Cycles = Time * 2 * sin(pi / N / 2)
alpha = 2.25
beta = 8


def fouriertransform(X):
    """fourier transform for sinusoidal initial condition"""
    indices = arange(N)
    return array([sqrt(2. / N) * sum(dot(X, sin(indices * k * pi / N))) for k
        in indices])


def evolve(XV, t=0):
    a = []
    for i in range(N):
        if (i == 0) or (i == N-1) :
            a.append(0)
        else:
            #nonlinear terms
            #nonlinear = beta * ((XV[i+1] - XV[i])**3 - (XV[i] - XV[i-1])**3)
            nonlinear = alpha * ((XV[i+1] - XV[i])**2 - (XV[i] - XV[i-1])**2)
            #nonlinear = 0

            xdotdot = (XV[i+1] + XV[i-1] - 2 * XV[i]) + nonlinear
            a.append(xdotdot)
    return concatenate((XV[N:], a))


def energyk(akdot, ak, k):
    omegak = 2 * sin(pi * k * 1./ N / 2)
    return .5 * akdot**2 + .5 * (omegak * ak)**2


def energyx(X, V, i):
    if i < len(X)-1:
        return .5 * V[i]**2 + .5 * ((X[i+1] -X[i])**2 + (X[i] - X[i-1])**2)
    else:
        return 0


# the meat of the code is in this very single line
XVs = odeint(evolve, XVs[0], Time)


Xs, Vs = split(array(XVs).T, 2)
Xs, Vs = Xs.T, Vs.T



#plotting

#1. plotting the movement of the waves
#todo: create movie
for i in range(Ncycle):
    if not i % 10:
        plot(Xs[i])


#2. energy_x
figure()
plot(Cycles, energyx(array(Xs).T, array(Vs).T, 2))
title("Energy of the 2nd chain vs cycles")


#3. energy of the normal modes
figure()
Xks = array([fouriertransform(i) for i in Xs]).T
Vks = array([fouriertransform(i) for i in Vs]).T

# plotting the normal modes from 1 to 5
for i in range(1,11):
    plot(Cycles,energyk(Vks[i], Xks[i], i))
ylabel("Energy")
xlabel("cycles")


raw_input()
exit()
