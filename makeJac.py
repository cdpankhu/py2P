# makeJac.py
import sys

from os.path import dirname, join

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


def makeJacVSC(ppc, vref, angles):

    busext = ppc["bus"]
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

    npvpq = len(pvpq)
    npq = len(pq)

    j1 = 0
    j2 = npvpq
    j3 = j2
    j4 = j3 + npq

    # Voltage angle sensitivity at bus i for P injection at bus j
    va_p = Jinv[j1:j2, j1:j2]
    # Voltage angle sensitivity at bus i for Q injection at bus j
    va_q = Jinv[j1:j2, j3:j4]
    # Voltage magnitude sensitivity at bus i for P injection at bus j
    vm_p = Jinv[j3:j4, j1:j2]
    # Voltage magnitude sensitivity at bus i for Q injection at bus j
    vm_q = Jinv[j3:j4, j3:j4]

    va_p_sc = {}
    va_q_sc = {}
    vm_p_sc = {}
    vm_q_sc = {}
    for i in busext[:, 0]:
        for j in busext[:, 0]:
            va_p_sc[i, j] = 0
            va_q_sc[i, j] = 0
            vm_p_sc[i, j] = 0
            vm_q_sc[i, j] = 0

    for i in range(npvpq):
        for j in range(npvpq):
            va_p_sc[pvpq[i]+1, pvpq[j]+1] = va_p[i, j]

    for i in range(npvpq):
        for j in range(npq):
            va_q_sc[pvpq[i]+1, pq[j]+1] = va_q[i, j]

    for i in range(npq):
        for j in range(npvpq):
            vm_p_sc[pq[i]+1, pvpq[j]+1] = vm_p[i, j]

    for i in range(npq):
        for j in range(npq):
            vm_q_sc[pq[i]+1, pq[j]+1] = vm_q[i, j]

    return vm_p_sc, vm_q_sc, va_p_sc, va_q_sc
