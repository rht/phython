from math import sqrt, cos, sin, atan
from random import random
import matplotlib.pyplot as P

try:
    import psyco
    psyco.full()
except ImportError:
    print 'Psyco not installed, the program will just run slower'

#contants
#c = 299792458.6
#h = 6.6 / (10 ** 34)
NUM=100
numarray = range(NUM)
numarraymin1 = range(NUM-1)
class Particles:
    def __init__(self):
        self.x= [(i+random())/NUM for i in numarray]
        self.c = [signchooser(i) for i in numarray]
        self.a = [self.c[i] * (sum(self.c[:i]) - sum(self.c[i+1:])) for i in numarray]
        self.v = [ (2.0*random() -1) for i in numarray]
        k = self.kineticenergy()
        self.v = [ 0.2 * i / sqrt(k) for i in self.v]
        u = self.potentialenergy()
        print u, self.kineticenergy()
        self.karray = [k]
        self.timearray = [0]
        self.t = 0

    def potentialenergy(self):
        def potential(particle):
            return -1 * sum([ self.c[i] * abs(self.x[particle] - self.x[i]) for i in numarray])
        return 0.5  * 0.01* sum(self.c[i]*potential(i) for i in numarray)

    def kineticenergy(self):
         return 0.5 * 0.1* sum([(self.v[i])**2 for i in numarray])

    def quadratic(self,a,b,c):
        d = b*b-4*a*c
        if (d >= 0):
            if a==0:
                return -c*1.0/b
            else:
                if c==0:
                    solarray = [-b/a,-b/a]
                else:
                    solarray = [((sqrt(d)-b)/(2*a)), ((-sqrt(d)-b)/(2*a))]
                solmin = min(solarray)
                solmax = max(solarray)
                if solmin > 0:
                    return solmin
                elif solmax > 0:
                    return solmax
                else:
                    return 20.0
        else:
            return 10.0

    def deltat(self):
        t1 = [self.quadratic((self.a[i+1]-self.a[i])/2.0,(self.v[i+1]-self.v[i]),(self.x[i+1] - self.x[i])) for i in numarraymin1]
        t1min = min(t1)
        for i in numarraymin1:
            if t1[i] == t1min:
                index = i
        #left
        xmin = min(self.x)
        for i in numarray:
            if xmin == self.x[i]:
                index2 = i
        tl = self.quadratic(self.a[index2]/2.0,self.v[index2],self.x[index2])

        #right
        xmax = max(self.x)
        for i in numarray:
            if xmax == self.x[i]:
                index3 = i
        tr = self.quadratic(self.a[index3]/2.0,self.v[index3],(self.x[index3]-1))
        
        tmin = min(t1min,tl,tr)
        #print t1min, "         ", tl,"                ", tr, index
        if tmin==t1min:return t1min,"",index
        elif tmin==tl:return tl,"two",index2
        elif tmin==tr:return tr,"three",index3
        
    def update(self):
        dt, state, index = self.deltat()
        self.x = [(self.x[i] + self.v[i] * dt + 0.5 * self.a[i] * dt * dt) for i in numarray]
        self.v = [(self.v[i] + self.a[i] * dt) for i in numarray]
        self.t+=dt
        self.timearray.append(self.t)
        self.karray.append(self.kineticenergy())
        #case2 and 3
        if (state=="two") or state=="three":
            self.v[index] = -1.0 *self.v[index]
        else:
            self.x[index],self.x[index+1] = self.x[index+1],self.x[index]
            self.v[index],self.v[index+1] = self.v[index+1],self.v[index]
            self.c[index],self.c[index+1] = self.c[index+1],self.c[index]
            correction = 2*self.c[index]*self.c[index+1]
            self.a[index],self.a[index+1] = self.a[index+1]-correction,self.a[index]+correction



def signchooser(x):
    return ((x<(NUM/2)) and 1) or (-1)

def simulate():
    ps = Particles()
    #testpos = [ps.x[48]]
    testpos2 = [ps.x[49]]
    testpos3 = [ps.x[50]]
    #testpos4 = [ps.x[51]]
    for i in xrange(800000):
    #while ps.t < 1:
        ps.update()
        #print ps.a[49]
        #testpos.append(ps.x[48])
        testpos2.append(ps.x[49])
        testpos3.append(ps.x[50])
        #testpos4.append(ps.x[51])
    print ps.c
    P.subplot(221)
    #P.plot(ps.v,ps.x)
    P.scatter(ps.v,ps.x)
    P.subplot(222)
    #P.plot(ps.timearray,testpos)
    P.plot(ps.timearray,testpos2)
    P.plot(ps.timearray,testpos3)
    #P.plot(ps.timearray,testpos4)
    P.subplot(223)
    P.plot(ps.timearray,ps.karray)
    #P.axis('equal')
    P.show()
simulate()
