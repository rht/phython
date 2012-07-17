from pylab import *
from scipy.stats import norm
from scipy.linalg.matfuncs import expm2
import os

def phi(x):
    k = 5
    #return norm.pdf(x-4) * exp(1j * k * x) + norm.pdf(x-40) * exp(-1j * k * x)
    return norm.pdf(x-4) + norm.pdf(x-40)

def phi2(x):
    k = 5
    #return norm.pdf(x-4.1) * exp(1j * k * x) + norm.pdf(x-39.9) * exp(-1j * k * x)
    return norm.pdf(x-4.1) + norm.pdf(x-39.9)



X = arange(-200,200,.5)
phis = [phi(X),phi2(X)]

#def evolve(initial,H,tfinal,timestep):
    #def propagator(H,t):
        #return expm2(-1j*H*t) #expm2: compute the matrix exponential using eigenvalue decomposition
    #time = arange(0,tfinal,timestep)
    #states = [initial]
    #U = propagator(H,timestep)
    ##evolve the system
    #for i in range(len(time)-1):
        #states.append(array(dot(U,states[-1])))
    #return time, states


def laplacematrix(dx,nx):
    return matrix((diag(ones(nx-1),1) + diag(ones(nx-1),-1) + diag(-2*ones(nx))) / (dx*dx))

def evolve(phis):
    lastphi = phis[-1]
    size = len(lastphi)
    #newphi = array(dot(laplacematrix(.1, size), lastphi)).flatten() + 2 * lastphi - phis[-2]
    newphi = .1 * .1 *array(dot(laplacematrix(.1, size), lastphi)).flatten() + 2 * lastphi - phis[-2]
    return array(newphi)

for i in range(30):
    phis.append(evolve(phis))

def createvideo(spectrums, plotter):
        #http://dawes.wordpress.com/2007/12/04/animating-png-files/
        #http://stackoverflow.com/questions/4092927/generating-movie-from-python-without-saving-individual-frames-to-files
        #http://www.scipy.org/Cookbook/Matplotlib/Animations
    import tempfile
    directory = tempfile.gettempdir()
    command = ('ffmpeg','-i', directory + '/%03d.png', 'out.avi', '-r', '5')
    #convert -delay 50 Th*.JPG anim.mpg

    ymin, ymax = min(spectrums[0]), max(spectrums[0])
    for i in range(len(spectrums)):
        figure()
        ylim(ymin, ymax)
        plotter(spectrums[i])
        filename = directory + '/%03d.png'%i
        savefig(filename)
        print('Wrote file '+ filename)
        clf()
    os.spawnvp(os.P_WAIT, 'ffmpeg', command)

createvideo(phis,plot)
#show()

#for i in phis:
    #figure()
    #plot(X,i)

#show()
