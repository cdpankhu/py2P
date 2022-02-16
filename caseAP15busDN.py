# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Power flow data for 15 bus, 3 generator case.
"""

from numpy import array


def caseAP15busDN():
    """Power flow data for 9 bus, 3 generator case.
    Please see L{caseformat} for details on the case file format.
    Based on data from Joe H. Chow's book, p. 70.
    @return: Power flow data for 9 bus, 3 generator case.
    """
    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 1.0

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc["bus"] = array([
        [1,  3, 0,      0,      0, 0,      1, 1, 0, 12.66, 1, 1.1, 0.9],
        [2,  1, 0.7936, 0.1855, 0, 0.0011, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [3,  1, 0,      0,      0, 0.0028, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [4,  1, 0.0201, 0.0084, 0, 0.0024, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [5,  1, 0.0173, 0.0043, 0, 0.0004, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [6,  1, 0.0291, 0.0073, 0, 0.0008, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [7,  1, 0.0219, 0.0055, 0, 0.0006, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [8,  1, 0.0219, 0.0019, 0, 0.0006, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [9,  1, 0.0235, 0.0059, 0, 0.0012, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [10, 1, 0.0229,  0.0142, 0, 0.0004, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [11, 1, 0.0217, 0.0065, 0, 0.0004, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [12, 1, 0,      0,      0, 0.0001, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [13, 1, 0.6219, 0.1291, 0, 0.0001, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [14, 1, 0.0014, 0.0008, 0, 0.0002, 1, 1, 0, 12.66, 1, 1.1, 0.9],
        [15, 1, 0.0224, 0.0083, 0, 0.0001, 1, 1, 0, 12.66, 1, 1.1, 0.9]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    ppc["gen"] = array([
        [1,  2, 0, 5, -5, 1, 1, 1, 5, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1,  1, 0, 2,    -2,    1, 1, 1, 2,    0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0],
        [12, 0.256, 0, 0.4,  -0.4,  1, 1, 1, 0.4,  0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0]
    ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([
        [1,  2,  0.001,  0.12,   0, 2, 250, 250, 0, 0, 1, -360, 360],
        [2,  3,  0.0883, 0.1262, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [3,  4,  0.1384, 0.1978, 0, 0.256, 150, 150, 0, 0, 1, -360, 360],
        [4,  5,  0.0191, 0.0273, 0, 0.256, 300, 300, 0, 0, 1, -360, 360],
        [5,  6,  0.0175, 0.0251, 0, 0.256, 150, 150, 0, 0, 1, -360, 360],
        [6,  7,  0.0482, 0.0689, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [4,  8,  0.0407, 0.0582, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [8,  9,  0.0523, 0.0747, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [8,  10, 0.01,   0.0143, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [10, 11, 0.0241, 0.0345, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [11, 12, 0.0103, 0.0148, 0, 0.256, 250, 250, 0, 0, 1, -360, 360],
        [1,  13, 0.001,  0.12,   0, 1, 250, 250, 0, 0, 1, -360, 360],
        [13, 14, 0.1559, 0.1119, 0, 0.204, 250, 250, 0, 0, 1, -360, 360],
        [14, 15, 0.0953, 0.0684, 0, 0.204, 250, 250, 0, 0, 1, -360, 360]
    ])

    ##-----  OPF Data  -----##
    ## area data
    # area refbus
    ppc["areas"] = array([
        [1, 5]
    ])

    ## generator cost data
    # 1 startup shutdown n x1 y1 ... xn yn
    # 2 startup shutdown n c(n-1) ... c0
    ppc["gencost"] = array([
        [2, 0, 0, 3, 0, 50, 0],
        [2, 0, 0, 3, 0, 20, 0],
        [2, 0, 0, 3, 0, 10, 0]
    ])

    return ppc
