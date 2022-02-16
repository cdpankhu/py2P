# Coefficients.py
from numpy import diag, zeros, dot, transpose, delete, absolute
from numpy.linalg import inv
from cmath import phase, rect
from math import sqrt


def calculateptdf(lines, buses, root):
    numlines = len(lines)
    numbus = len(buses)

    lineb = []
    linea = zeros([numlines, numbus])
    for li in lines:
        lineb.append(1/lines[li].x)
        linea[li-1, lines[li].tbus-1] = 1
        linea[li-1, lines[li].fbus-1] = -1

    lineat = transpose(linea)
    bdiag = diag(lineb)
    bmat = dot(dot(lineat, bdiag), linea)
    bmat[root-1, :] = 0
    bmat[:, root-1] = 0
    bmat[root-1, root-1] = 1
    bmatinv = inv(bmat)
    bmatinv[root-1, root-1] = 0
    isf = dot(dot(bdiag, linea), bmatinv)

    ptdf = {}
    for li in lines:
        for s in buses:
            for b in buses:
                ptdf[li, s, b] = isf[int(li-1), int(s-1)] \
                    - isf[int(li-1), int(b-1)]

    return ptdf


def anglerecovery(lines, buses, fp, fq, v):
    numlines = len(lines)
    numbus = len(buses)

    c = zeros([numbus, numlines])
    for li in lines:
        c[lines[li].tbus-1, li-1] = 1
        c[lines[li].fbus-1, li-1] = -1

    # Delete first row assuming it is the root bus, need to fix for different
    # root bus
    bt = transpose(delete(c, 0, 0))

    # The original reference paper considers potential of a mesh network
    # Assuming that the network is radial, a spanning tree includes all edges
    # As such it is not necessary to calculate a spanning tree otherwise

    # Calculate beta(i,j) from the relaxed optimal solution
    # beta(i,j) = angle(vi - (zij*)x(Sij))
    # vi is the voltage magnitude squared
    betat = zeros([numlines, 1])
    for li in lines:
        z = complex(lines[li].r, lines[li].x)
        betat[li-1] = phase(v[lines[li].tbus]
                            - z.conjugate()*complex(-fp[li], -fq[li]))

    # Under a meshed network, need to ensure that the conditions for a unique
    # solution are met such that the sum of angles in any cycle is 0.
    # Since there are no cycles, beta(perpendicular) is 0, no cycles exist, and
    # the condition is guranteed.
    theta = dot(inv(bt), betat)

    thetaout = {}
    thetaout[1] = 0
    for i in range(numbus-1):
        thetaout[i+2] = theta[i, 0]

    return thetaout


def makeYbus(buses, lines):

    Ybus = {}
    for i in buses:
        for j in buses:
            Ybus[i, j] = 0

    for li in lines:
        tbus = lines[li].tbus
        fbus = lines[li].fbus
        Ybus[tbus, fbus] = -1/complex(lines[li].r, lines[li].x)
        Ybus[fbus, tbus] = Ybus[tbus, fbus]

    for b in buses:
        Ybus[b, b] = complex(buses[b].Gs, buses[b].Bs)
        for li in lines:
            tbus = lines[li].tbus
            fbus = lines[li].fbus
            if tbus == b or fbus == b:
                Ybus[b, b] += -Ybus[fbus, tbus] + 1j*lines[li].b/2

    return Ybus


def calculatevsc(buses, lines, v, angles, root):
    Ybus = makeYbus(buses, lines)

    numbus = len(buses)
    mat = zeros([2*numbus, 2*numbus])

    vrect = {}
    for b in buses:
        vrect[b] = rect(sqrt(v[b]), angles[b])

    for i in buses:
        # Adding to real equations of matrix
        for j in buses:
            if i != root:
                temp = Ybus[i, j]*vrect[j]
                row = int((i-1)*2)
                col = int((i-1)*2)
                mat[row, col+1] += temp.real
                mat[row+1, col+1] += temp.imag
            if j != root:
                temp = Ybus[i, j]*vrect[i].conjugate()
                row = int((i-1)*2)
                col = int((j-1)*2)
                mat[row, col] += temp.real
                mat[row+1, col] += temp.imag

    mat = delete(mat, [0, 1], 0)
    mat = delete(mat, [0, 1], 1)

    invmat = inv(mat)

    vscprime = zeros([2*(numbus-1), numbus-1])
    vsc = {}

    for b in range(numbus-1):
        temp = zeros([2*(numbus-1), 1])
        temp[b*2] = 1
        # temp[b*2 + 1] = 1
        vscprime[:, b] = dot(invmat, temp)[:, 0]

    for i in buses:
        for j in buses:
            vsc[i, j] = 0
            if i != root and j != root:
                row = int((i-2)*2)
                col = int(j-2)
                temp = (vrect[j].conjugate()
                        * complex(vscprime[row, col],
                                  vscprime[row+1, col]))
                vsc[i, j] = (1/absolute(vrect[j]))*temp.real

    return vsc, vscprime, vrect
