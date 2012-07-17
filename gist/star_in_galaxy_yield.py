from pylab import *
#by fsz, tweaked

def orbit():

    '''numerically solves and plots the solution of  the three coupled diff eqs
    derived in problem 3-5 of Pset 8 to plot the orbit of stars in
    a galaxy with a spherically symettric galactic potential, with va/vc =.451'''
    
    
    #setting up some initial values
    delta_zeta =  0.001 #sufficiently small enough to counterbalance the errors
                        #of Euler's method
    k = 0.451
    zeta = 0.0  # dimensionless time 
    theta = 0.0 
    xi = 1.0 #dimensionless radius
    phi = 0.0 #proportional to radial velocity

    #for the purpose of graphing the results
    zetamax = 20*pi
    time = arange(zeta,zetamax,delta_zeta)
    theta_r_v = [array([theta,xi,phi])]

     def propagator(state):
        theta, xi, phi = state
        return state + array([k/xi**2, phi, (k**2/xi**3)-(1/xi)]) *delta_zeta

   
    #loop, continue integrating as long as the dimensionless time is < 20 pi
    for t in time:
        # no explicit time dependence -> the system is conservative!
        theta_r_v.append(propagator(theta_r_v[-1]))


    #graphing
    title('Plot of stellar orbit for k= 0.451 for theta = 0 to 20 pi')
    polar(*array(theta_r_v).T[:-1])
    show()
    



orbit()
