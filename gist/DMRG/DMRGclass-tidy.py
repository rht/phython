#!/usr/bin/python
# -*- coding: utf-8 -*-
from mpmath import *

# from random import uniform, sample  #imported, but unused

from pylab import *
from mpl_toolkits.mplot3d import Axes3D

# import os, sys

import scipy.linalg as ln
import scipy.sparse as sparse
import scipy.sparse.linalg
from time import time

mp.dps = 60

begintime = time()
sp = sparse.csr_matrix([[0, 1], [0, 0]])
sm = sp.H
sz = .5 * sparse.csr_matrix([[1, 0], [0, -1]])


def I(m):
    # return matrix(identity(m))
    return sparse.csr_matrix(identity(m))


def MatConj(A, B):
    # return dot(dot(A,sparse.csr_matrix(B)),sparse.csr_matrix(A).H)
    return A * B * A.H


def kronecker(A, B, C, D):
    return sparse.kron(sparse.kron(sparse.kron(A, B), C), D)


def bigkronecker(A, m):
    return (A if m is 0 else sparse.kron(bigkronecker(A, m - 1), A))


def localop(Operator, site, chainlength, dimpersite):
    II = I(dimpersite)
    sitesbefore = site - 1
    sitesafter = chainlength - site

    O1 = bigkronecker(II, sitesbefore)
    O2 = bigkronecker(II, sitesafter)
    return sparse.kron(sparse.kron(O1, site), O2)


## DMRG Class

class DMRG:

    def __init__(self, B):
        self.spl = sparse.csr_matrix([[0, 1], [0, 0]])
        self.spr = self.spl
        self.sml = sparse.csr_matrix([[0, 0], [1, 0]])
        self.smr = self.sml
        self.szl = .5 * sparse.csr_matrix([[1, 0], [0, -1]])
        self.szr = self.szl
        self.GS = []
        self.GSwfn = 0
        self.HB = 0. * sparse.csr_matrix([[1, 0], [0, -1]])
        self.steprenormalization = 0
        self.entropy = []

    def renormalize(self, m):

        self.steprenormalization += 1

        global sp
        global sm
        global sz

        #print shape(self.spr), shape(self.smr), shape(self.szr), shape(self.HB)

        B = shape(self.HB)[0]

        HBB = kronecker(self.HB, I(2), I(2), I(B)) + kronecker(I(B),
                I(2), I(2), self.HB) + kronecker(self.szl, sz, I(2),
                I(B)) + .5 * kronecker(self.spl, sm, I(2), I(B)) + .5 \
            * kronecker(self.sml, sp, I(2), I(B)) + kronecker(I(B), sz,
                sz, I(B)) + .5 * kronecker(I(B), sp, sm, I(B)) + .5 \
            * kronecker(I(B), sm, sp, I(B)) + kronecker(I(B), I(2), sz,
                self.szr) + .5 * kronecker(I(B), I(2), sp, self.smr) \
            + .5 * kronecker(I(B), I(2), sm, self.spr)

        Hbb = sparse.kron(self.HB, I(2)) + sparse.kron(self.szr, sz) \
            + .5 * sparse.kron(self.spr, sm) + .5 \
            * sparse.kron(self.smr, sp)

        # gs, gsenergy = ln.eigh(HBB)[1][:, 0].reshape(2*B, 2*B), ln.eigh(HBB)[0][0]

        result = sparse.linalg.eigsh(HBB)
        (gs, gsenergy) = (result[1][:, 0].reshape(2 * B, 2 * B),
                          result[0][0])
        self.GS.append(gsenergy / (self.steprenormalization + 1.))
        self.GSwfn = gs

        U = ln.svd(gs)

        # U = sparse.linalg.svds(gs,k=4)

        rho = U[1] * U[1].T
        print (shape(U[2]), shape(rho))
        self.entropy.append(sum(log(rho)))

        T = sparse.csr_matrix(U[0])[:, 0:m].H

        #print shape(T), shape(Hbb)

        (self.HB, self.szl, self.spl, self.szr, self.spr) = (MatConj(T,
                Hbb), MatConj(T, sparse.kron(self.szl, I(2))),
                MatConj(T, sparse.kron(self.spl, I(2))), MatConj(T,
                sparse.kron(I(2), self.szr)), MatConj(T,
                sparse.kron(I(2), self.spr)))

        self.sml = self.spl.H
        self.smr = self.spr.H

        return gs

    def pplot(self):
        length = len(self.GS)
        plot(arange(length), array(self.GS))
        xlabel('Iterations')
        ylabel('GS eenrgy density')

        show()


spin12 = DMRG(0)
spin12.renormalize(4)
spin12.renormalize(8)
spin12.renormalize(16)
spin12.renormalize(32)
spin12.renormalize(64)
spin12.renormalize(128)
spin12.renormalize(128)
spin12.renormalize(128)
spin12.renormalize(128)
for i in range(10):
    spin12.renormalize(4)

print time() - begintime

spin12.pplot()
