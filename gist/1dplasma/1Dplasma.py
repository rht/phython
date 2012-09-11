from __future__ import division
from math import sqrt
from numpy import roots
from random import random,choice
import matplotlib.pyplot as P


#contants
#c = 299792458.6
#h = 6.6 / (10 ** 34)
NUM=100
numarray = range(NUM)
numarraymin1 = range(NUM-1)
class Particles:
    def __init__(self):
        self.x= [(i+random())/NUM for i in numarray]
        self.c = [((i<(NUM/2)) and 1) or (-1)       for i in numarray]
        self.a = [0.01*self.c[i] * (sum(self.c[:i]) - sum(self.c[i+1:])) for i in numarray]
        while (self.potentialenergy() > 100):
            self.x = [(i+random())/NUM for i in numarray]
        self.v = [ (2.0*random() -1) for i in numarray]
        k = self.kineticenergy()

        #renormalization
        self.v = [ 0.2 * i / sqrt(2.0*k/NUM) for i in self.v]
        u = self.potentialenergy()
        print u, self.kineticenergy()
        self.karray = [k]
        self.timearray = [0]
        self.t = 0

    def computea(self):
        return [0.01*self.c[i] * (sum(self.c[:i]) - sum(self.c[i+1:])) for i in numarray]

    def potentialenergy(self):
        def potential(particle):
            return -1 * sum([ self.c[i] * abs(self.x[particle] - self.x[i]) for i in numarray])
        return 0.5  * 0.01* sum(self.c[i]*potential(i) for i in numarray)

    def kineticenergy(self):
         return 0.5 *  sum([(self.v[i])**2 for i in numarray])

    def quadratic(self,a,b,c):
        d = b*b-4*a*c
        if (d >= 0):
            solarray = roots([a,b,c])
            solmin = min(solarray)
            solmax = max(solarray)
            if solmin > 0: return solmin
            elif solmax > 0: return solmax
            else: return 300.0
        else: #discriminant < 0
            return 300.0

    def deltat(self):
        t1 = [self.quadratic(((self.a[i+1]-self.a[i])/2.0)
                             ,(self.v[i+1]-self.v[i]),
                             (self.x[i+1] - self.x[i])) for i in numarraymin1]
        t1min = min(t1)
        for i in numarraymin1:
            if t1[i] == min(t1):
                index = i
        #left
        tl = self.quadratic(self.a[0]/2.0,self.v[0],self.x[0])

        #right
        tr = self.quadratic(self.a[NUM-1]/2.0,self.v[NUM-1],(self.x[NUM-1] -1))
        
        tmin = min(t1min,tl,tr)
        #print t1min, "         ", tl,"                ", tr, index
        if tmin==t1min:return t1min,"one",index
        elif tmin==tl:return tl,"two",0
        elif tmin==tr:return tr,"three",NUM-1
        
    def update(self):
        #TRY TO FIT THE UNITS
        dt, state, index = self.deltat()
        self.a = self.computea()
        self.x = [(self.x[i] + self.v[i] * dt + (0.5 * self.a[i] * dt * dt)) for i in numarray]
        self.v = [(self.v[i] + self.a[i] * dt) for i in numarray]
        self.t += dt
        self.timearray.append(self.t)
        self.karray.append(self.kineticenergy())
        #case2 and 3
        if (state=="two") or (state=="three"):
            self.v[index] = -1.0 *self.v[index]
            if self.x[index] < 2e-19: self.x[index]=0
            print self.x[index],"HEREE",index, self.v[index]
        elif (state=="one"):
            #self.x[index],self.x[index+1] = self.x[index+1],self.x[index]
            self.v[index],self.v[index+1] = self.v[index+1],self.v[index]
            self.c[index],self.c[index+1] = self.c[index+1],self.c[index]
            #correction = 2*self.c[index]*self.c[index+1]
           # self.a[index],self.a[index+1] = self.a[index+1]-correction,self.a[index]+correction



def simulate():
    length =200
    sysnum = 1
    simlength = xrange(length)
    ps = [Particles() for i in xrange(sysnum)]

    #for i in simlength:
    while ps[0].t < 1:
        for j in ps:
            j.update()

    averagekarray = [0]*length
    averagetimearray = [0]*length
    ## for i in simlength:
    ##     for j in ps:
    ##         averagekarray[i]+= j.karray[i] / sysnum
    ##         averagetimearray[i]+=j.timearray[i] / sysnum
    print ps[0].c
    P.subplot(221)
    P.scatter(ps[0].x,ps[0].v)    
    
    P.subplot(223)
    P.plot(ps[0].timearray,[10]*len(ps[0].timearray))
    P.plot(ps[0].timearray,ps[0].karray)
    P.subplot(224)
    #P.plot(averagetimearray,averagekarray)
    #P.axis('equal')
    P.show()
simulate()
