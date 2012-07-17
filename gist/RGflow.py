from pylab import *
from mpl_toolkits.mplot3d import Axes3D
# http://arxiv.org/abs/1202.5984
# http://code.google.com/p/pydmrg/
# about pycool http://arxiv.org/pdf/1201.5029.pdf

param0s = [array([1,2,.2])]
params = [param0s]
params.append([array([1,2.1,0])])
params.append([array([.5,2.1,0])])
params.append([array([.5,2.1,0])])
d = 4.01
b = .9
coef = array([b ** -2, b ** (d-4), b ** (2 * d - 6)])

for i in range(40):
    for ps in params:
        ps.append(coef * ps[-1])

fig = figure()
ax = fig.gca(projection='3d')
for ps in params:
    ax.plot(*array(ps[1:]).T)
xlabel("M")
ylabel("$\lambda$")
#zlabel("C")
show()
