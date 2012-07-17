from pylab import *

#calculate beam's mass

#contants
dt = 0.001
#c = 299792458.6
#h = 6.6 / (10 ** 34)

class Particle:
    def __init__(self):
        self.v = array([random(), 10. + random()])
        self.x = array([10.0 + random(), random()])
        self.xideal = array([10.0,0])
        self.t = 0

    def radius(self):
        return norm(self.x)
    def radiusideal(self):
        return norm(self.xideal)
    def vradial(self):
        return norm(self.v) * cos(self.calc_theta(self.v) - self.theta())

    def calc_theta(self,foo):
        return arctan(foo[1] / foo[0])
    def theta(self):
        return self.calc_theta(self.x)
    def thetaideal(self):
        return self.calc_theta(self.xideal
    def vtangential(self):
        return norm(self.v) * sin(self.calc_theta(self.v) - self.theta())

    def force(self):
        def Fmagnet():
            return [-self.v[1],self.v[0]]
        def RF(state):
            if (state=='on') and (self.x[0] > 0) and (abs(self.theta()) < 0.05):
                #print(self.theta())
                return array([0,10.0*cos(self.t)])
            return array([0,0])
        def Fradiation():
            return array([-0.01 * self.v[0], -0.01 * self.v[1]])
        return array([sum(i) for i in zip(Fmagnet(),RF("on"),Fradiation())])

        
    def update(self):
        self.x += self.v * dt
        self.v += self.force() * dt
        self.xideal = array([10.0 - sin(self.t) , cos(self.t)])
        self.t += dt

    def deviation(self):
        return [self.radius() - self.radiusideal(),
                self.theta() - self.thetaideal()]

def average(values):
    return sum(values,0.0) / len(values)

def simulate():
    numarray = range(40)
    particles = [Particle() for i in numarray]
    timearray = [dt*i for i in range(20000)]
    avemomentum = []
    avex = []
    avey = []
    for i in timearray:
        for j in particles:
            j.update()
        avemomentum.append(average([norm(i.v) for i in particles]))
        avex.append(average([i.x[0] for i in particles]))
        avey.append(average([i.x[1] for i in particles]))

    momentumarray = [particle.vradial() for particle in particles]
    deviationarray = [particle.deviation()[0] for particle in particles]
    #subplot(211)
    scatter(deviationarray,momentumarray)
    #subplot(212)
    #plot(numarray,momentumarray)
    #plot(avex,avey)
    #plot(timearray,avemomentum)
    axis('equal')
    show()
simulate()
