# Coefficients.py
from numpy import diag, zeros, dot, transpose
from numpy.linalg import inv


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
                ptdf[li, s, b] = isf[li-1, s-1] - isf[li-1, b-1]

    return ptdf
