from pylab import *
# http://www.pma.caltech.edu/Courses/ph136/yr2004/book03/chap07/0207.1.pdf
# ugh title gets cropped, here is the solution
# http://stackoverflow.com/questions/8802918/my-matplotlib-title-gets-cropped

N = 1000
D = 3.  # size of aperture
L = 20  # distance to screen
k = 6.  # wavenumber


def plot_cornu(N, size, distance, k):
    X = zeros(N)  # horizontal axis of our phasors
    Y = zeros(N)  # vertical axis of our phasors
    for i in range(-int(N/2), int(N/2)):
    #for i in range(N - 1):
        phase = pow(size * i, 2) / 2 / distance
        X[i+1] = X[i] + cos(phase)
        Y[i+1] = Y[i] + sin(phase)
    plot(X[N/2+1:],Y[N/2+1:], 'b')
    plot(X[:N/2+1],Y[:N/2+1], 'b')

#plot_cornu(N, D, 20, k)
#show();exit()
for l in linspace(4,200,100):
    clf()
    plot_cornu(N, D, l, k)
    title('"sumoverhistories" at distance %dm \nsize of aperture %d k %d'
            % (l, D, k), fontsize=12)
    savefig("cornu%d.png" %(l))
#show()

