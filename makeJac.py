# makeJac.py
import sys

from os.path import dirname, join

from time import time

from numpy import r_, zeros, array
from numpy import complex128

from scipy.sparse import hstack, vstack
from scipy.linalg import inv
from cmath import rect
from math import sqrt

from pypower.dSbus_dV import dSbus_dV

from pypower.bustypes import bustypes
from pypower.ext2int import ext2int
from pypower.makeYbus import makeYbus


def makeJacV(ppc, vref, angles):

    ppc = ext2int(ppc)
    baseMVA, bus, gen, branch = \
        ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"]

    numbus = len(bus)

    # get bus index lists of each type of bus
    ref, pv, pq = bustypes(bus, gen)

    Ybus, Tf, Ty = makeYbus(baseMVA, bus, branch)

    # set up indexing for updating V
    pvpq = r_[pv, pq]

    V = zeros(numbus, dtype=complex128)
    for b in range(numbus):
        V[b] = rect(sqrt(vref[b+1]), angles[b+1])

    # evaluate Jacobian
    dS_dVm, dS_dVa = dSbus_dV(Ybus, V)

    J11 = dS_dVa[array([pvpq]).T, pvpq].real
    J12 = dS_dVm[array([pvpq]).T, pq].real
    J21 = dS_dVa[array([pq]).T, pvpq].imag
    J22 = dS_dVm[array([pq]).T, pq].imag

    J = vstack([
            hstack([J11, J12]),
            hstack([J21, J22])
        ], format="csr")

    J = J.todense()

    Jinv = inv(J)

    return J, pv, pq, pvpq, J11, J12, J21, J22, dS_dVm, dS_dVa
