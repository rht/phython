#xmax,ymax 20, number of particles 25, distance R, velocity v, energy E,U,T,K
#standard euler algorithm
from pylab import *
#import matplotlib.animation as animation

class Particles:
    def __init__(self):
        #self.deltat = 0.1
        self.deltat = 1
        self.num = 10

    def Upair(self,i,j):
        return 0.001 / sqrt((self.r[0][i] - self.r[0][j]) ** 2 + \
               (self.r[1][i] - self.r[1][j]) ** 2)

    def Uone(self,foo):
        return sum([self.Upair(foo,j) for j in filter(lambda x: x!= foo, range(self.num))])

    def Uencalc(self):
        self.Uen = array([self.Uone(i) for i in range(self.num)])

    def generate(self):
        self.v = array([random(self.num), random(self.num)])
        self.r = array([random(self.num), random(self.num)]) * 20
        self.Uencalc()
        self.Utotal = sum(self.Uen)

    def kinetic(self):
        return 0.04 * sum(self.v ** 2)

    def vmag(self):
        return self.v[0]**2 + self.v[1]**2

    def aveU(self):
        return 0.02 * self.Utotal


    def propagate(self):
        Ubefore = self.Uen
        deltar = self.v * self.deltat
        self.r = (self.r + deltar) % 20
        self.Uencalc()
        self.v += .02 * (Ubefore - self.Uen) / deltar

#http://stackoverflow.com/questions/9401658/matplotlib-animating-a-scatter-plot
#ion()
box = Particles()
box.generate()
#scat = scatter(*box.r)

#numframes = 10
numframes = 100000
#def update_plot(i,data,scat):
#def update_plot():
    #box.propagate()
    #scat = scatter(*box.r)
    #show()

    #return scat.set_array(*box.r),

for i in range(numframes):
    box.propagate()
    #scatter(*box.r)
    #show()
figure()
hist(box.v[0]**2 + box.v[1]**2)
show()
#ani = animation.FuncAnimation(figure(), update_plot, frames=range(numframes),
                              #fargs=(,scat))

#show()
