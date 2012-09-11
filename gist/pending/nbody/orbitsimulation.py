from pylab import plot, show
from math import sqrt

vel, pos = [0.0, 6.3], [1, 0]
posx, posy = [pos[0]], [pos[1]]
t, dt = 0, 0.03

def updateprop(vel, pos, accel, posx, posy):
    for i in [0, 1]:
        vel[i] += accel[i] * dt
        pos[i] += vel[i] * dt
    posx.append(pos[0])
    posy.append(pos[1])
    return vel, pos, posx, posy

def acceleration(pos):
    a = [0, 0]
    r2 = pos[0] ** 2 + pos[1] ** 2
    for i in [0, 1]:
        a[i] = -40 * pos[i] / (r2 ** 1.5) * (1 - 0.04 / sqrt(r2))
    return a

while t <= 10:
    accel = acceleration(pos)
    vel, pos, posx, posy = updateprop(vel, pos, accel, posx, posy)
    t += dt

plot(posx, posy)
show()
