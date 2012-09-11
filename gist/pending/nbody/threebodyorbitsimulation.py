from pylab import plot, show

vel1, vel2 = [0.0, 6.3], [0.0, 2.1]
pos1, pos2 = [1, 0], [3, 0]
posx1, posy1 = [pos1[0]], [pos1[1]]
posx2, posy2 = [pos2[0]], [pos2[1]]
t, dt = 0, 0.02

def updateprop(vel, pos, accel, posx, posy):
    for i in [0, 1]:
        vel[i] += accel[i] * dt
        pos[i] += vel[i] * dt
    posx.append(pos[0])
    posy.append(pos[1])
    return vel, pos, posx, posy

def acceleration(pos1, pos2):
    a = [0, 0]
    r2 = pos1[0] ** 2 + pos1[1] ** 2
    r2rel = (pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2
    for i in [0, 1]:
        a[i] = (-40 * pos1[i] / (r2 ** 1.5)) \
               + 0.04 * (pos2[i] - pos1[i]) / (r2rel ** 1.5)
    return a

while t <= 10:
    accel1, accel2 = acceleration(pos1, pos2), acceleration(pos2, pos1)
    (vel1, pos1, posx1, posy1), (vel2, pos2, posx2, posy2) = updateprop(vel1, pos1, accel1, posx1, posy1), \
                                                         updateprop(vel2, pos2, accel2, posx2, posy2)
    t += dt

plot(posx1, posy1, 'b', posx2, posy2, 'r')
show()
