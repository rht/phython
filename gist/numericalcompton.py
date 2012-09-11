from physics import *
#jul25 3.15PM - 7.05PM

mass = .1


def photon_pol(k, pol):
    k0, k1, k2, k3 = k
    vec_k = array([k1, k2, k3])
    vec_pol1 = vec_k - dot(vec_k, array([1, 0, 0]))
    vec_pol1 = vec_pol1 / norm(vec_pol1)
    vec_pol2 = cross(vec_k, vec_pol1)
    vec_pol2 /= norm(vec_pol2)

    # http://www.physics.indiana.edu/~dermisek/QFT_09/qft-I-11-4p.pdf
    # see also tom banks 43
    # circular polarization  -> not really necessary
#    if pol == 1:  # left-handed
        #return array([0] + list(1 / sqrt(2) * (vec_pol1 - 1j * vec_pol2)))
    #else:
        #return array([0] + list(1 / sqrt(2) * (vec_pol1 + 1j * vec_pol2)))
    if pol == 1:
        return array([0] + list(vec_pol1))
    else:
        return array([0] + list(vec_pol2))


def scattering_matrixoldest((p, pprime, k, kprime), mu, nu, spin1, spin2, pol1, pol2):
    # Peskin 5.74
    m = mass
    propagator = ( ((gamma(mu) * slash(k) * gamma(nu)) + 2. * gamma(mu) * p[nu])
                / (2. * fourdot(p, k)) +
                (-gamma(nu) * slash(kprime) * gamma(mu) + 2. * gamma(nu) * p[mu])
                / (-2. * fourdot(p, kprime))
                )
    return (-1j * dot(photon_pol(kprime, pol1).conj(), photon_pol(k, pol2)) *
            dot(matrix(fermion_ubar(pprime, spin1, m)),
                dot(propagator, fermion_u(p, spin2, m)).T ))


#alternative implementation, without mu, nu
def scattering_matrixold((p, pprime, k, kprime), spin1, spin2, pol1, pol2):
    # Peskin 5.74
    m = mass
    eps_star = photon_pol(kprime, pol1).conj()
    eps = photon_pol(k, pol2)
    sl_eps_star = slash(eps_star)
    sl_eps = slash(eps)
    propagator = ((sl_eps_star * slash(k) * sl_eps + 2. *
                   sl_eps_star * fourdot(eps, p))
                / (2. * fourdot(p, k)) +
                (-sl_eps * slash(kprime) * sl_eps_star +
                    2. * sl_eps * fourdot(eps_star, p))
                / (-2. * fourdot(p, kprime))
                )
    return (-1j * dot(matrix(fermion_ubar(pprime, spin1, m)),
                dot(propagator, fermion_u(p, spin2, m)).T ))


#newwer alternative implementation, without mu, nu
def scattering_matrix((p, pprime, k, kprime), spin1, spin2, pol1, pol2):
    # Peskin before 5.74
    m = mass
    eps_star = photon_pol(kprime, pol1).conj()
    eps = photon_pol(k, pol2)
    sl_eps_star = slash(eps_star)
    sl_eps = slash(eps)

    propagator = ((sl_eps_star * (slash(p) + slash(k) + m) * sl_eps)
                / (2. * fourdot(p, k)) +
                (sl_eps * (slash(p) - slash(kprime) + m) * sl_eps_star)
                / (-2. * fourdot(p, kprime))
                )
    return (-1j * dot(matrix(fermion_ubar(pprime, spin1, m)),
                dot(propagator, fermion_u(p, spin2, m)).T ))



def scattering_amplitudeold((p, pprime, k, kprime)):
    indices = [1, 2]
    indices4vector = [0, 1, 2, 3]
    amplitudelist = []
    for mu in indices4vector:
        for nu in indices4vector:
            for spin1 in indices:
                for spin2 in indices:
                    for alpha in indices4vector:
                        for beta in indices4vector:
                            for pol1 in indices:
                                for pol2 in indices:
                                    amp = scattering_matrix((p,
                                        pprime, k, kprime), mu,
                                        nu, spin1, spin2, pol1, pol2) * scattering_matrix((p,
                                            pprime, k, kprime),
                                            alpha, beta, spin2, spin1,
                                            pol2, pol1)
                                    amplitudelist.append(amp)
    return sum(amplitudelist) / 4.


# alternative implementation, without mu, nu
def scattering_amplitude((p, pprime, k, kprime)):
    indices = [1, 2]
    amplitudelist = []
    for spin1 in indices:
        for spin2 in indices:
            for pol1 in indices:
                for pol2 in indices:
                    M = scattering_matrix((p,
                        pprime, k, kprime), spin1, spin2, pol1, pol2)
                    amp = M * M.conj()
                    #print amp
                    amplitudelist.append(amp)
    return sum(amplitudelist) / 4.



def scattering_amplitude_theoretical((p, pprime, k, kprime)):
    m = mass
    # peskin 5.87
    term1 = fourdot(p, kprime) / fourdot(p, k)
    term2 = fourdot(p, k) / fourdot(p, kprime)
    term3 = 2 * m * m * (1. / fourdot(p, k) - 1. / fourdot(p, kprime))
    term4 = m ** 4 * (1. / fourdot(p, k) - 1. / fourdot(p, kprime)) ** 2
    return 2 * (term1 + term2 + term3 + term4)

#print gamma(1)
p = fourvector([0,0,0], mass)
pprime = fourvector([1,2,6], mass)
k = fourvector([4,2,0], 0)
kprime = fourvector([6,2,6], 0)
mu = 2
nu = 1
spin1 = 1
spin2 = 2
pol1, pol2 = 1, 1


print scattering_amplitude((p, pprime, k, kprime))
print scattering_amplitude_theoretical((p, pprime, k, kprime))
