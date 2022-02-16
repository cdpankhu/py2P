# SC.py
from pandas import DataFrame
from py2P.NetworkLoad import networkload
from py2P.DLMP import calculatedlmp
from py2P.Coefficients import calculateptdf, makeYbus
from py2P.makeModel import makeModel
from math import pow, sqrt, floor
from gurobipy import Model, GRB
from time import strftime
from os import makedirs


def sc(testsystem):
    buses, lines, generators, datamat, ppc = networkload(testsystem)
    # buses[1].Vmax = 1
    # buses[1].Vmin = 1
    windset = {}
    windex = 1

    for g in generators:
        if generators[g].gtype == "Wind":
            windset[windex] = generators[g].gindex
            windex += 1

    root = -1
    for g in generators:
        if generators[g].gtype == "Root":
            root = generators[g].location

    if root == -1:
        print("No generator at root node")
        generators[1].gtype = "Root"
        root = 1

    gensetP = []
    gensetU = [root]
    for g in generators:
        if not g == root:
            gensetP.append(g)

    # Defining the set of generator buses and demand buses
    B_g = []
    for g in generators:
        if not (generators[g].location in B_g):
            B_g.append(generators[g].location)
            # For generation bus, set demand to 0
            # This needs to be changed with updated model for storage, DR, etc.
            buses[generators[g].location].Pd = 0
    B_d = [i for i in buses if i not in B_g]

    # Defining matrix of generators at each bus
    B_gn = {}
    for b in buses:
        B_gn[b] = []
    for g in generators:
        B_gn[generators[g].location].append(g)

    # Share of demand from utility versus peers
    partlevel = 1.0
    gam = {}
    for b in buses:
        gam[b] = partlevel

    pbd = {}
    for b in buses:
        for k in buses:
            pbd[b, k] = 0

    sbase = generators[1].mBase

    model = makeModel(buses, generators, lines,
                      gensetP, gensetU, SC=1, gam=gam)

    model.optimize()

    status = model.Status

    var = model.getVars()
    constrs = model.getConstrs()

    # Extracting variable results into arrays
    obj_value = model.ObjVal
    varCount = 0
    v = {}
    for b in buses:
        v[b] = var[varCount].x
        varCount += 1
    a = {}
    for li in lines:
        a[li] = var[varCount].x
        varCount += 1
    fp = {}
    for li in lines:
        fp[li] = var[varCount].x*sbase
        varCount += 1
    fq = {}
    for li in lines:
        fq[li] = var[varCount].x*sbase
        varCount += 1
    pg = {}
    for g in generators:
        pg[g] = var[varCount].x*sbase
        varCount += 1
    qg = {}
    for g in generators:
        qg[g] = var[varCount].x*sbase
        varCount += 1
    oc = var[varCount].x
    varCount += 1
    ocgen = {}
    for g in generators:
        ocgen[g] = var[varCount].x
        varCount += 1
    p = {}
    for b in buses:
        p[b] = var[varCount].x*sbase
        varCount += 1
    pnm = {}
    for i in buses:
        for j in buses:
            if abs(var[varCount].x) > 1e-6:
                pnm[i, j] = var[varCount].x*sbase
            else:
                pnm[i, j] = 0
            varCount += 1

    pgb = {}
    qgb = {}
    for i in buses:
        pgb[i] = 0
        qgb[i] = 0
    for b in buses:
        if B_gn[b]:
            pgb[b] = sum(pg[g] for g in B_gn[b])
            qgb[b] = sum(qg[g] for g in B_gn[b])

    constrCount = 0 + len(buses)*len(buses) + len(buses)*len(buses) + \
        len(buses) + len(buses) + 1 + len(generators) + len(lines) + 1
    dual_pbalance = {}
    for b in buses:
        dual_pbalance[b] = constrs[constrCount].Pi
        constrCount += 1

    DLMPInfo = DataFrame()
    NodeInfo = DataFrame()
    LineInfo = DataFrame()
    GenInfo = DataFrame()

    dispatch = pg
    print("dispatch1=:", dispatch)
    # trade_scale is kW
    # for i in dispatch:
    #     dispatch[i] = floor(dispatch[i]*10**3)/(10**3)
    SMP = generators[root].cost[1]
    status1, dlmp, pg, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
        calculatedlmp(
            dispatch, buses, generators, lines, SMP, gensetP, gensetU)

    cn = {}
    for i in buses:
        for j in buses:
            cn[i, j] = (dlmp[j] - dlmp[i])/2

    line_loading = {}
    loading = {}
    flow = {}
    for li in lines:
        line_loading[li] = 1 - round(
            (lines[li].u - sqrt(
                (pow(LineInfo[li]['fp'], 2)
                 + pow(LineInfo[li]['fq'], 2))))/lines[li].u, 3)
        loading[li] = LineInfo[li]['flow'] / lines[li].u
        flow[li] = LineInfo[li]['flow']

    v_marginrate = {}
    voltage = {}
    for b in buses:
        v_marginrate[b] = round(0.1 - sqrt(NodeInfo[b]['v']), 3)/0.1
        voltage[b] = sqrt(NodeInfo[b]['v'])

    bus_frame = DataFrame([line_loading, v_marginrate, dlmp],
                          index=['line_loading', 'v_marginrate', 'dlmp'])
    # loading_frame = DataFrame([line_loading], index='line_loading')
    # vmargin_frame = DataFrame([v_marginrate], index='v_marginrate')

    ep_u = {}
    ep_p = {}
    nuc = {}
    for b in buses:
        ep_u[b] = 0
        ep_p[b] = 0
        nuc[b] = 0

    for m in B_d:
        ep_p[m] = sum(sum(generators[g].cost[1]*pnm[n, m] for g in B_gn[n]
                          if not (g in gensetU)) for n in B_g if
                      sum(i for i in B_gn[n] if not (i in gensetU)))
        ep_u[m] = (1-gam[m])*buses[m].Pd*dlmp[m]
        nuc[m] = sum(pnm[n, m]*(dlmp[m]-dlmp[n]) for n in B_g)

    ur = {}
    cp = {}
    for b in buses:
        ur[b] = ep_u[b] + nuc[b]
        cp[b] = ep_u[b] + ep_p[b] + nuc[b]

    gc = ocgen
    sum_gc_p = sum(gc[g] for g in gensetP)
    sum_gc_u = sum(gc[g] for g in gensetU)
    sum_dgr = sum(gc[g] for g in gensetP)
    sum_dgp = sum_dgr - sum_gc_p
    up = sum(ur[b] for b in buses) - sum_gc_u
    sum_cp = sum(cp[b] for b in buses)
    sum_ep_p = sum(ep_p[b] for b in buses)
    sum_ep_u = sum(ep_u[b] for b in buses)
    sum_nuc = sum(nuc[b] for b in buses)

    cashflow_mat = DataFrame([ep_p, ep_u, nuc, ur, cp],
                             index=['ep_p', 'ep_u', 'nuc', 'ur', 'cp'])
    cashflow_sum = DataFrame([sum_cp, sum_ep_p, sum_ep_u, sum_nuc,
                              sum_gc_p, sum_gc_u, up],
                             index=['sum_cp', 'sum_ep_p', 'sum_ep_u', 'sum_nuc',
                                    'sum_ocp', 'sum_ocu', 'up'])
    data_mat = DataFrame([pnm], index=["trades_dis"])
    data_mat = data_mat.transpose()

    print("partlevel = ", partlevel)
    print("status = ", status)
    print("max loading = ", max(loading[li] for li in lines)*100, "%")
    print("max voltage = ", round(max(voltage[b] for b in buses), 3))
    print("min voltage = ", round(min(voltage[b] for b in buses), 3))
    print("max dlmp = ", round(max(dlmp[b] for b in buses), 3))
    print("min dlmp = ", round(min(dlmp[b] for b in buses), 3))
    method = "sc"

    basefile = "C:\\Users\\Colton\\github\\py2P\\results\\"
    dirname = basefile+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)

    busfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)+"_busframe.csv"
    nodefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)+"_nodeframe.csv"
    linefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)+"_lineframe.csv"
    dlmpfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)+"_dlmpframe.csv"
    genfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)+"_genframe.csv"
    cashflowmatfile = dirname+"\\"+method+strftime("%Y-%m-%d-%H-%M-%S_")+"_" + \
        testsystem+"_"+str(100*partlevel)+"_cfmat.csv"
    cashflowsumfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+str(100*partlevel)+"_cfsum.csv"
    datafile = dirname+"\\"+method+strftime("%Y-%m-%d-%H-%M-%S_")+"_" + \
        testsystem+"_"+str(100*partlevel)+"_data.csv"

    makedirs(dirname, exist_ok=True)
    bus_frame.to_csv(busfile)
    cashflow_mat.to_csv(cashflowmatfile)
    cashflow_sum.to_csv(cashflowsumfile)
    NodeInfo.to_csv(nodefile)
    LineInfo.to_csv(linefile)
    DLMPInfo.to_csv(dlmpfile)
    GenInfo.to_csv(genfile)
    data_mat.to_csv(datafile)

    return pnm, dlmp
