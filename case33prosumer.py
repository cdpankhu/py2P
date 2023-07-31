# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Power flow data for 15 bus, 3 generator case.
"""

from numpy import array


def case33prosumer():
    """Power flow data for 9 bus, 3 generator case.
    Please see L{caseformat} for details on the case file format.
    Based on data from Joe H. Chow's book, p. 70.
    @return: Power flow data for 9 bus, 3 generator case.
    """
    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 10.0

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    # Pd and Qd are in kVA, convert below to MVA
    ppc["bus"] = array([
        [1,    3,  0,      0,      0,  0,  1,  1,  0,  12.66,  1,  1,      1],
        [2,    1,  100,    60,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [3,    1,  90,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [4,    1,  120,    80,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [5,    1,  60,     30,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [6,    1,  60,     20,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [7,    1,  200,    100,    0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [8,    1,  200,    100,    0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [9,    1,  60,     20,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [10,   1,  60,     20,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [11,   1,  45,     30,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [12,   1,  60,     35,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [13,   1,  60,     35,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [14,   1,  120,    80,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [15,   1,  60,     10,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [16,   1,  60,     20,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [17,   1,  60,     20,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [18,   1,  90,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [19,   1,  90,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [20,   1,  90,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [21,   1,  90,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [22,   1,  90,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [23,   1,  90,     50,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [24,   1,  420,    200,    0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [25,   1,  420,    200,    0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [26,   1,  60,     25,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [27,   1,  60,     25,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [28,   1,  60,     20,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [29,   1,  120,    70,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [30,   1,  200,    600,    0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [31,   1,  150,    70,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [32,   1,  210,    100,    0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9],
        [33,   1,  60,     40,     0,  0,  1,  1,  0,  12.66,  1,  1.1,    0.9]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    # Limits in MVA
    ppc["gen"] = array([
        [1, 0,  0,  10, -10,    1,  100,    1,  10,
         -10,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        # Placeholder agent for bus 1
        [1, 0,  0,  0,   0,  1,  100,    1,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        [2, 2.5,  0,  2,  -2,     1,  100,    1,  2.5,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        [6, 1.04,  0,  3,  -3,     1,  100,    1,  3.5,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        [18, 0.72,  0,  1, -1,     1,  100,    1,  1.6,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        [22, 0.36,  0,  2, -2,     1,  100,    1,  2.3,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        [25, 1.999,  0,  2.5, -2.5,     1,  100,    1,  2.9,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
        [31, 1.3,  0,  2, -2,     1,  100,    1,  2.5,
         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0]
    ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    # branch r and x in ohms, converted to pu below
    ppc["branch"] = array([
        [1,     2,  0.0922, 0.0470, 0,  6,  0,  0,  0,  0,  1,  -360,   360],
        [2,     3,  0.4930, 0.2511, 0,  6,  0,  0,  0,  0,  1,  -360,   360],
        [3,     4,  0.3660, 0.1864, 0,  6,  0,  0,  0,  0,  1,  -360,   360],
        [4,     5,  0.3811, 0.1941, 0,  6,  0,  0,  0,  0,  1,  -360,   360],
        [5,     6,  0.8190, 0.7070, 0,  6,  0,  0,  0,  0,  1,  -360,   360],
        [6,     7,  0.1872, 0.6188, 0,  1.5,  0,  0,  0,  0,  1,  -360,   360],
        [7,     8,  0.7114, 0.2351, 0,  1.5,  0,  0,  0,  0,  1,  -360,   360],
        [8,     9,  1.0300, 0.7400, 0,  1.5,  0,  0,  0,  0,  1,  -360,   360],
        [9,     10, 1.0440, 0.7400, 0,  1.5,  0,  0,  0,  0,  1,  -360,   360],
        [10,    11, 0.1966, 0.0650, 0,  1.5,  0,  0,  0,  0,  1,  -360,   360],
        [11,    12, 0.3744, 0.1238, 0,  1.5,  0,  0,  0,  0,  1,  -360,   360],
        [12,    13, 1.4680, 1.1550, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [13,    14, 0.5416, 0.7129, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [14,    15, 0.5910, 0.5260, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [15,    16, 0.7463, 0.5450, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [16,    17, 1.2890, 1.7210, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [17,    18, 0.7320, 0.5740, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [2,     19, 0.1640, 0.1565, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [19,    20, 1.5042, 1.3554, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [20,    21, 0.4095, 0.4784, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [21,    22, 0.7089, 0.9373, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [3,     23, 0.4512, 0.3083, 0,  2,  0,  0,  0,  0,  1,  -360,   360],
        [23,    24, 0.8980, 0.7091, 0,  2,  0,  0,  0,  0,  1,  -360,   360],
        [24,    25, 0.8960, 0.7011, 0,  2,  0,  0,  0,  0,  1,  -360,   360],
        [6,     26, 0.2030, 0.1034, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [26,    27, 0.2842, 0.1447, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [27,    28, 1.0590, 0.9337, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [28,    29, 0.8042, 0.7006, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [29,    30, 0.5075, 0.2585, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [30,    31, 0.9744, 0.9630, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [31,    32, 0.3105, 0.3619, 0,  1,  0,  0,  0,  0,  1,  -360,   360],
        [32,    33, 0.3410, 0.5302, 0,  1,  0,  0,  0,  0,  1,  -360,   360]
        #Tie lines
        # [21	8	2.0000	2.0000	0	0	0	0	0	0	0 -360	360],
        # [9	15	2.0000	2.0000	0	0	0	0	0	0	0 -360	360],
        # [12	22	2.0000	2.0000	0	0	0	0	0	0	0 -360	360],
        # [18	33	0.5000	0.5000	0	0	0	0	0	0	0 -360	360],
        # [25	29	0.5000	0.5000	0	0	0	0	0	0	0 -360	360]
    ])

    ## generator cost data
    # 1 startup shutdown n x1 y1 ... xn yn
    # 2 startup shutdown n c(n-1) ... c0
    ppc["gencost"] = array([
        [2, 0, 0, 3, 0, 7.65, 0],
        [2, 0, 0, 3, 0, 7.65, 0],
        [2, 0, 0, 3, 0.05, 7.2, 0],
        [2, 0, 0, 3, 0.02, 7.6, 0],
        [2, 0, 0, 3, 0.03, 7.1, 0],
        [2, 0, 0, 3, 0.08, 7.4, 0],
        [2, 0, 0, 3, 0.03, 6.9, 0],
        [2, 0, 0, 3, 0.04, 7.4, 0]
    ])

    Vbase = ppc["bus"][0, 9]
    sBase = ppc["baseMVA"]
    zBase = (Vbase*Vbase)/sBase
    for i in range(len(ppc["branch"])):
        ppc["branch"][i, 2] = ppc["branch"][i, 2]/zBase
        ppc["branch"][i, 3] = ppc["branch"][i, 3]/zBase

    # Demand scalar
    demandScale = 2
    for i in range(len(ppc["bus"])):
        # Double base demand
        ppc["bus"][i, 2] = demandScale*ppc["bus"][i, 2] / 1e3
        ppc["bus"][i, 3] = demandScale*ppc["bus"][i, 3] / 1e3

    return ppc
