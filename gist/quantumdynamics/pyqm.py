#references
#Time evolution
##Numerical approaches to time evolution of complex quantum systems
##introduces Crank-Nicolshon scheme
##http://www.sciencedirect.com/science/article/pii/S0375960109004927

from pylab import *
#from fipy import *
#from mpl_toolkits.mplot3d import Axes3D

##Pauli spin matrices
sx = matrix([[0, 1],[ 1, 0]])
sy = matrix([[0, -1j],[1j, 0]])
sz = matrix([[1, 0],[0, -1]])
smin = matrix([[0,1],[0,0]])
smax = matrix([[0,0],[1,0]])

#helper functions
def renormalize(x): return x/norm(x)

def create1Dmesh(dx,nx):
    L = nx * dx
    #return linspace(-L,L,nx)
    return arange(-L/2.,L/2.,dx)

def create2Dmesh(dx,nx):
    L = nx * dx
    X = arange(-(L-1)/2, (L+1)/2,dx)
    return meshgrid(X,X) 

def create2DmeshK(dx,nx):
    Kx = linspace(-pi, pi, num = nx)
    return Kx, Kx

def gaussian(x, sigma):
    return  exp(-x*x/(2.*sigma**2))/sqrt(2.*pi*sigma**2)

def xtop(wavefn):
    '''convert position space wavefn into momentum space wavefn'''
    #http://docs.scipy.org/doc/numpy/reference/routines.fft.html#background-information
    return fftshift(abs(fftn(wavefn)))

def ptox(wavefn):
    return fftshift(abs(ifftn(wavefn)))

def laplacematrix(dx,nx):
    return matrix((diag(ones(nx-1),1) + diag(ones(nx-1),-1) + diag(-2*ones(nx))) / (dx*dx))
def laplace(x,dx):
    ##from scipy.ndimage.filters import laplace
    #the line above doesn't work because our wavefunction is complex
    return convolve(x,[1,-2,1],mode='same') /(dx*dx)


#laplacian matrix
#http://stackoverflow.com/questions/1764859/how-to-compute-laplacian-of-a-field
#http://en.wikipedia.org/wiki/Laplacian_matrix


def plot2Dspectrum(Hamiltonian, K):
    spectrum = array([ sort(real(eigvals(Hamiltonian(i)))) for i in K])
    figure()
    for i in range(shape(spectrum)[1]):
        plot(K, spectrum[:,i])

def animate1D(X,Z):
    fig = figure()
    #line,foo = plot(X,Z[0],X,potential(X))
    line, = plot(X,Z[0])
    for i in Z[1:]:
        waitforbuttonpress()
        line.set_ydata(i)
        draw()

def plotsurfaceanimate(X,Y,Z):
    fig = figure()
    ax = Axes3D(fig)
    ax.set_zlim3d(0,.12)
    for i in Z:
        ax.plot_surface(X,Y,i)
        draw()

def createvideo(spectrums,plotter):
        #http://dawes.wordpress.com/2007/12/04/animating-png-files/
        #http://stackoverflow.com/questions/4092927/generating-movie-from-python-without-saving-individual-frames-to-files
        #http://www.scipy.org/Cookbook/Matplotlib/Animations
    command = ('ffmpeg','-i', '%03d.png', 'out.mp4', '-r', '25')
    for i in range(len(spectrums)):
        plotter(spectrums[i])
        filename = '%03d'%i + '.png'
        savefig(filename, dpi=100)
        print 'Wrote file', filename
        clf()
    os.spawnvp(os.P_WAIT, 'ffmpeg', command)


