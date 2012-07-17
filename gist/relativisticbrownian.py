from pylab import *

N = 1000
def sran(x=None): return 2*random(x)-1

positions = [sran(2), sran(2)]
positions2 = list(positions)
def getnextpos(pos, pos2):
    deltaangle = sran() * pi
    def getnewpos(pos, frac=1):
        r1 = pos[-2]
        r2 = pos[-1]
        r = r2 - r1
        angle = arctan(r[1] / r[0])
        # implement "relativistic spotlight effect"
        newangle = angle + frac * deltaangle
        return r2 + norm(r) * array([cos(newangle), sin(newangle)])
    return getnewpos(pos, frac=.5), getnewpos(pos2)

for i in range(N):
    new1, new2 = getnextpos(positions, positions2)
    positions.append(new1)
    positions2.append(new2)


plot(*array(positions).T)
plot(*array(positions2).T)
plot(positions[0][0], positions[0][1], 'ro')
plot(positions[1][0], positions[1][1], 'go')
show()
