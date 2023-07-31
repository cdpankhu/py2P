# bSC.py
from pandas import DataFrame
from py2P.py2Pbase.bP2PTradeType import Agent
from py2P.py2Pbase.bNetworkLoad import networkload
from py2P.py2Pbase.bTradeFunction import generatetrade, selecttrade
from py2P.py2Pbase.bDLMP import calculatedlmp
from math import pow, sqrt
from gurobipy import Model, GRB
from time import strftime
from os import makedirs


def sc(testsystem):
    buses, lines, generators, datamat = networkload(testsystem)
    buses[1].Vmax = 1
    buses[1].Vmin = 1
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

    m = Model()

    v = m.addVars((b for b in buses), lb=0, name="v")
    # Variable for current square
    a = m.addVars((li for li in lines), lb=0, name="a")
    fp = m.addVars((li for li in lines), name="fp", lb=float('-inf'))
    fq = m.addVars((li for li in lines), name="fq", lb=float('-inf'))
    pg = m.addVars((g for g in generators), name="pg", lb=float('-inf'))
    qg = m.addVars((g for g in generators), name="qg", lb=float('-inf'))
    oc = m.addVar(lb=0, name="oc")
    ocgen = m.addVars((g for g in generators), lb=0, name="ocgen")
    # P2P variables
    p = m.addVars((b for b in buses), name='p', lb=float('-inf'))
    pnm = m.addVars((b for b in buses), (b for b in buses), name='pnm',
                    lb=float('-inf'))

    m.update()

    # Set objective functions
    obj = -oc
    m.setObjective(obj, GRB.MAXIMIZE)

    # Create operating constraints
    # Constraints for bidirectional trading
    m.addConstrs((pnm[i, j] + pnm[j, i] == 0 for i in buses for j in buses),
                 name='bidirection')

    m.addConstrs((pnm[i, j] >= 0 for i in B_g for j in buses), name='signg')
    m.addConstrs((pnm[i, j] <= 0 for i in B_d for j in buses), name='signd')
    m.addConstrs((p[i] == sum(pnm[i, j] for j in buses) for i in buses),
                 name='nodalbalance')

    m.addConstrs((sum(pg[g] for g in B_gn[n]
                      if not (g in gensetU)) == p[n] for n in B_g),
                 name='genpower')
    m.addConstrs(((-gam[n]) * buses[n].Pd == p[n] for n in B_d), name='demand')

    m.addConstr((oc == sum(ocgen[g] for g in generators)),
                name='operatingcost')

    # Polynomial generation cost, single order
    m.addConstrs(
        (ocgen[g] == (generators[g].cost[1]*pg[g] + generators[g].cost[2])
            for g in generators), name="gencost"
    )

    # Apparent power flow limits on the receiving node
    m.addConstrs(
        ((fp[li]*fp[li] + fq[li]*fq[li])
            <= lines[li].u*lines[li].u for li in lines),
        name="linecapfw"
    )
    # Apparent power flow limits on the sending node
    m.addConstrs(
        ((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r)
            + (fq[li] - a[li]*lines[li].x)*(fq[li] - a[li]*lines[li].x)
            <= (lines[li].u)*(lines[li].u) for li in lines),
        name="linecapbw"
    )
    # Constraint relating the line flows and voltages
    m.addConstrs(
        ((v[lines[li].fbus]
            - 2*(lines[li].r*fp[li] + lines[li].x*fq[li])
            + a[li]*(lines[li].r*lines[li].r + lines[li].x*lines[li].x))
            == v[lines[li].tbus] for li in lines),
        name="betweennodes"
    )
    # Second order conic constraint
    m.addConstrs(
        ((fp[li]*fp[li] + fq[li]*fq[li])
            <= v[lines[li].fbus]*a[li] for li in lines),
        name="socp"
    )
    # Voltage contraint for root node
    for g in generators:
        if generators[g].gtype == "Root":
            m.addConstr(v[generators[g].location] == 1, name="rootvoltage")
    # Active power balance constraint
    m.addConstrs(
        ((
            sum(fp[li] for li in buses[b].outline)
            - sum(fp[li] - lines[li].r*a[li] for li in buses[b].inline)
            - sum(pg[g] for g in B_gn[b]) + buses[b].Pd + v[b]*buses[b].Gs
            ) == 0 for b in buses),
        name="pbalance"
    )
    # Reactive power balance constraint
    m.addConstrs(
        ((
            sum(fq[li] for li in buses[b].outline)
            - sum(fq[li] - lines[li].x*a[li] for li in buses[b].inline)
            - sum(qg[g] for g in B_gn[b]) + buses[b].Qd - v[b]*buses[b].Bs
            ) == 0 for b in buses),
        name="qbalance"
    )
    m.addConstrs(
        (v[b] <= buses[b].Vmax*buses[b].Vmax for b in buses),
        name="vmax"
    )
    m.addConstrs(
        (v[b] >= buses[b].Vmin*buses[b].Vmin for b in buses),
        name="vmin"
    )
    m.addConstrs(
        (pg[g] >= generators[g].Pmin for g in generators),
        name="pgmin"
    )
    m.addConstrs(
        (pg[g] <= generators[g].Pmax for g in generators),
        name="pgmax"
    )
    m.addConstrs(
        (qg[g] >= generators[g].Qmin for g in generators),
        name="qgmin"
    )
    m.addConstrs(
        (qg[g] <= generators[g].Qmax for g in generators),
        name="qgmax"
    )
    m.Params.QCPDual = 1
    m.update()

    m.optimize()

    status = m.Status

    var = m.getVars()

    # Extracting variable results into arrays
    obj_value = m.ObjVal
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
        fp[li] = var[varCount].x
        varCount += 1
    fq = {}
    for li in lines:
        fq[li] = var[varCount].x
        varCount += 1
    pg = {}
    for g in generators:
        pg[g] = var[varCount].x
        varCount += 1
    qg = {}
    for g in generators:
        qg[g] = var[varCount].x
        varCount += 1
    oc = var[varCount].x
    varCount += 1
    ocgen = {}
    for g in generators:
        ocgen[g] = var[varCount].x
        varCount += 1
    p = {}
    for b in buses:
        p[b] = var[varCount].x
        varCount += 1
    pnm = {}
    for i in buses:
        for j in buses:
            if abs(var[varCount].x) > 1e-3:
                pnm[i, j] = var[varCount].x
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

    DLMPInfo = DataFrame()
    NodeInfo = DataFrame()
    LineInfo = DataFrame()
    GenInfo = DataFrame()

    dispatch = pg
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
                              sum_dgr, sum_dgp, up],
                             index=['sum_cp', 'sum_ep_p', 'sum_ep_u', 'sum_nuc',
                                    'sum_dgr', 'sum_dgp', 'up'])

    print("partlevel = ", partlevel)
    print("status = ", status)
    print("max loading = ", max(loading[li] for li in lines)*100, "%")
    print("max voltage = ", round(max(voltage[b] for b in buses), 3))
    print("min voltage = ", round(min(voltage[b] for b in buses), 3))
    print("max dlmp = ", round(max(dlmp[b] for b in buses), 3))
    print("min dlmp = ", round(min(dlmp[b] for b in buses), 3))
    method = "basesc"

    basefile = "N:\\Python\\py2P\\results\\"
    dirname = basefile+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem
    busfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem \
        + "_busframe.csv"
    nodefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem \
        + "_nodeframe.csv"
    linefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem \
        + "_lineframe.csv"
    dlmpfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem \
        + "_dlmpframe.csv"
    genfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem \
        + "_genframe.csv"
    cashflowmatfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" \
        + testsystem+"_cfmat.csv"
    cashflowsumfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_cfsum.csv"
    makedirs(dirname, exist_ok=True)
    bus_frame.to_csv(busfile)
    cashflow_mat.to_csv(cashflowmatfile)
    cashflow_sum.to_csv(cashflowsumfile)
    NodeInfo.to_csv(nodefile)
    LineInfo.to_csv(linefile)
    DLMPInfo.to_csv(dlmpfile)
    GenInfo.to_csv(genfile)

    return dispatch, dlmp, pnm
