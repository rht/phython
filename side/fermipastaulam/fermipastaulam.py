from pylab import *
from scipy.integrate import odeint
from numpy.fft import rfft
from pyqm import createvideo, createvideofromdirectory
import os
import time
# "The only computer experiments worth doing are those that yield a surprise"
# "Metropolis And von Neumann Install Awful Computer"
# Our problem turned out to have been felicitously chosen. The results were
# entirely different qualitatively from what even Fermi, with his great knowledge of wave motions, had expected.

#Computational physics was born! http://www-dft.ts.infn.it/~resta/tcse/lec1.pdf
# For more than two decades, most of genuine computational physics adressed classical statistical mechanics (or classical dynamical systems).
# Modern quantum mechanics of materials (alias the birth of computational electronic structure)
#http://stackoverflow.com/questions/6611678/need-fast-c-qt-qwt-scatter-plot


tic = time.time()

# initialization
N = 32
L = 1.
Ncycle = 30000
mode = 1
BC = "fixed"
#initial phase space coordinates
XVs = [concatenate((sin(mode * pi * linspace(0,L,num=N)), zeros(N)))]
dt = sqrt(.125)
Time = linspace(0, Ncycle * dt, num=Ncycle)
Cycles = Time * 2 * sin(pi / N / 2)
alpha = .25#1.25
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
            nonlinear = alpha * ((XV[(i+1) % N] - XV[i])**2 -
                                 (XV[i] - XV[(i-1) %N])**2)
            #nonlinear = 0

            xdotdot = (XV[(i+1) % N] + XV[(i-1) % N] - 2 * XV[i]) + nonlinear
            a.append(xdotdot)
    return concatenate((XV[N:], a))


def evolvePBC(XV, t=0):
    '''periodic boundary condition'''
    a = []
    for i in range(N):
        #nonlinear terms
        #nonlinear = beta * ((XV[i+1] - XV[i])**3 - (XV[i] - XV[i-1])**3)
        nonlinear = alpha * ((XV[(i+1) % N] - XV[i])**2 -
                             (XV[i] - XV[(i-1) % N ])**2)
        #nonlinear = 0

        xdotdot = (XV[(i+1) % N] + XV[(i-1) % N] - 2 * XV[i]) + nonlinear
        a.append(xdotdot)
    return concatenate((XV[N:], a))


def evolvevectorizedPBC(XV, t=0):
    def laplacematrix(dx,nx):
            return matrix((diag(ones(nx-1),1) + diag(ones(nx-1),-1) +
                diag(-2*ones(nx))) / (dx*dx))
    # from
    # http://iopscience.iop.org/0295-5075/64/5/606/pdf/0295-5075_64_5_606.pdf
    X = XV[:N]
    a = ravel(dot(laplacematrix(L / N, N), (X + alpha * X * X)))
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
#question 6
#XVs = odeint(evolvePBC, XVs[0], Time)


Xs, Vs = split(array(XVs).T, 2)
Xs, Vs = array(Xs.T), array(Vs.T)




dirname = "mode%dcycle%dalpha%.2f%s" %(mode, Ncycle, alpha, BC)
try:
    os.mkdir(dirname)
except:
    pass

#plotting

#1. plotting the movement of the waves
#todo: create movie
#for i in range(Ncycle):
    #if not i % 10:
        #plot(Xs[i])
contourf(Xs)
savefig(dirname + "/xsvst.png")


#2. energy_x
figure()
plot(Cycles, energyx(array(Xs).T, array(Vs).T, 2))
title("Energy of the 2nd chain vs cycles")
savefig(dirname + "/energyxvst.png")


#3. energy of the normal modes
figure()
Xks = array([fouriertransform(i) for i in Xs]).T
Vks = array([fouriertransform(i) for i in Vs]).T
# plotting the normal modes from 1 to 5
for i in range(1,11):
    plot(Cycles,energyk(Vks[i], Xks[i], i))
ylabel("Energy"); xlabel("cycles")
savefig(dirname + "/normalmodes.png")




# 6. plotting a travelling soliton
print "creating video"
figures = []
xsmin, xsmax = Xs.min(), Xs.max()

# it turns out python can't accomodate a list bigger than 3000
if Ncycle < 1000:
    for i in range(Ncycle):
        fig = figure()
        ylim(xsmin, xsmax)
        plot(Xs[i])
        print "video %d" %i
        figures.append(fig)
    createvideo(figures,prefix=0, outputdir=dirname)
else:
    #we need to parallelize the loop, or else it will be too slow
    #http://stackoverflow.com/questions/6652124/naive-and-easiest-way-to-decompose-independent-loop-into-parallel-threads-proces
    cycles = split(arange(Ncycle), Ncycle / 1000)
    import time
    import tempfile
    directory = tempfile.mkdtemp()
    from multiprocessing import Pool

    def f(c):
        for i in c:
            figure()
            ylim(xsmin, xsmax)
            plot(Xs[i])
            print "video %d" %i
            #pref = str(c[0]) + time.strftime("%b%d%Y")
            pref = ''
            filename = directory + '/%s%03d.png'%(pref, i)
            savefig(filename)
            clf()

    #pool = Pool(processes=4)              # start 4 worker processes
    #pool.map(f, cycles)
    #pool.close()
    #map(f, cycles)
    for i in cycles:
        f(i)
    createvideofromdirectory(directory, outputdir=dirname)
    #http://stackoverflow.com/questions/581851/in-python-how-do-i-make-a-temp-file-that-persists-until-the-next-run
#    import shutil
    #import subprocess
    #try:
        #subprocess.check_call(['/bin/echo', 'Directory:', directory])
    #finally:
            #shutil.rmtree(directory)


total_time = time.time() - tic
print "total time spent", total_time / 60., "min"

#storing datas
#import shelve
#data = shelve.open(dirname + "/data.dat")
#data['total_time'] = total_time
#data['Xs'] = Xs
#data.close()
