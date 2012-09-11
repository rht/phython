#why multidimensional monte carlo integral is more efficient
#http://farside.ph.utexas.edu/teaching/329/lectures/node109.html

#based on this pdf
#http://www.andrew.cmu.edu/user/cmorning/QM_MonteCarlo.pdf
#montepython    
#http://arxiv.org/pdf/physics/0609191v1


#euclidean action
from pylab import *
from scipy import *

epsilonomega = .25
kappa = 1/4. * epsilonomega **2
# xk = dk sqrt(epsilon hbar / m)
g = (1-kappa)/(1+kappa)

Nsweep = 10000 # number of iteration
dimension = 1000

def S(u):
    #return 1/2. * sum( (d[i+1] -d[i])**2 + kappa*(d[i+1]+d[i]) for i in range(len(d)-1)  )
    return sum([ u[j]**2 - g * (u[j]* u[j+1])  for j in range(len(u) -1)  ])

def deltaS(u,delta,j):
    return delta * (delta + 2*u[j] - g*(u[j-1] + u[j+1]))

#one sweep
def calculateev(observable,time):
    accumulator = 0
    U = random.uniform(-100.,100.,dimension)
    for i in range(Nsweep):
        for j in range(dimension-1):
            delta = random.uniform(-1.5,1.5)
            if random.random() < min(1,exp(-deltaS(U,delta,j))):
                accumulator += observable(time,U)
                U[j] += delta
    return accumulator

def pos(time,locations): return locations[time] * locations[0]
plot([calculateev(pos,i) for i in range(10)])
show()
