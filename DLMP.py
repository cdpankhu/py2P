# DLMP.py
from gurobipy import Model, GRB
from math import pow


def calculatedlmp(
        dispatch, buses, generators, lines, SMP, gensetP, gensetU):

    B_gn = {}
    for i in buses:
        B_gn[i] = []
    for g in generators:
        B_gn[generators[g].location].append(g)

    m = Model()
    # Create variables
    # Variable for voltage square
    v = m.addVars((b for b in buses), lb=0, name="v")
    # Variable for current square
    a = m.addVars((li for li in lines), lb=0, name="a")
    fp = m.addVars((li for li in lines), name="fp", lb=float('-inf'))
    fq = m.addVars((li for li in lines), name="fq", lb=float('-inf'))
    pg = m.addVars((g for g in generators), name="pg")
    qg = m.addVars((g for g in generators), name="qg")
    oc = m.addVar(lb=0, name="oc")
    ocgen = m.addVars((g for g in generators), lb=0, name="ocgen")

    m.update()

    # Set objective functions
    obj = -oc
    m.setObjective(obj, GRB.MAXIMIZE)

    m.update()

    # Create constraints
    # Define operating cost
    m.addConstr(oc == sum(ocgen[g] for g in generators), name="oc")
    # Define linear generation cost
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
    m.addConstrs(
        (dispatch[g] == pg[g] for g in gensetP),
        name="genpower"
    )
    m.Params.QCPDual = 1
    m.update()

    m.optimize()

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
    pq = {}
    for g in generators:
        pq[g] = var[varCount].x
        varCount += 1
    oc = var[varCount].x
    varCount += 1
    ocgen = {}
    for g in generators:
        ocgen[g] = var[varCount].x
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

    # Getting constraint duals
    constrs = m.getConstrs()
    qconstrs = m.getQConstrs()

    # Getting quadratic constraint duals
    qconstrCount = 0
    dual_linecapfw = {}
    for li in lines:
        dual_linecapfw[li] = qconstrs[qconstrCount].QCPi
        qconstrCount += 1
    dual_linecapbw = {}
    for li in lines:
        dual_linecapbw[li] = qconstrs[qconstrCount].QCPi
        qconstrCount += 1
    dual_socp = {}
    for li in lines:
        dual_socp[li] = qconstrs[qconstrCount].QCPi
        qconstrCount += 1

    # Getting linear constraint duals
    # Move count past oc, and ocgen
    constrCount = 0 + 1 + len(generators)
    dual_oc = constrs[0].Pi
    dual_betweennodes = {}
    for li in lines:
        dual_betweennodes[li] = constrs[constrCount].Pi
        constrCount += 1

    # Move count past root generator voltage constraint
    constrCount += 1

    dual_pbalance = {}
    for b in buses:
        dual_pbalance[b] = constrs[constrCount].Pi
        constrCount += 1
    dual_qbalance = {}
    for b in buses:
        dual_qbalance[b] = constrs[constrCount].Pi
        constrCount += 1
    dual_vmax = {}
    for b in buses:
        dual_vmax[b] = constrs[constrCount].Pi
        constrCount += 1
    dual_vmin = {}
    for b in buses:
        dual_vmin[b] = constrs[constrCount].Pi
        constrCount += 1
    dual_pgmin = {}
    for g in generators:
        dual_pgmin[g] = constrs[constrCount].Pi
        constrCount += 1
    dual_pgmax = {}
    for g in generators:
        dual_pgmax[g] = constrs[constrCount].Pi
        constrCount += 1
    dual_qgmin = {}
    for g in generators:
        dual_qgmin[g] = constrs[constrCount].Pi
        constrCount += 1
    dual_qgmax = {}
    for g in generators:
        dual_qgmax[g] = constrs[constrCount].Pi
        constrCount += 1

        # Deriving DLMP coefficients from model results
    a1 = {}
    a2 = {}
    a3 = {}
    a4 = {}
    a5 = {}
    for li in lines:
        a1[li] = ((
            (pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
            + a[li]*fq[li]*(pow(lines[li].r, 2) - pow(lines[li].x, 2))
            - 2*a[li]*fp[li]*lines[li].r*lines[li.x])
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2)))
            )

        a2[li] = ((
            (pow(fp[li], 2) + pow(fq[li], 2))*lines[li].r
            - a[li]*fp[li]*(pow(lines[li].r, 2) + pow(lines[li].x, 2)))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2)))
            )

        a3[li] = ((
            -(pow(fp[li], 2) + pow(fq[li], 2))*lines[li].r
            + a[li]*fp[li]*(pow(lines[li].r, 2) - pow(lines[li].x, 2))
            + 2*a[li]*fq[li]*lines[li].r*lines[li].x)
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2)))
            )

        a4[li] = ((
            2*(pow(fq[li], 3)*lines[li].r - pow(fp[li], 3)*lines[li].x)
            + 2*fp[li]*fq[li]*(fp[li]*lines[li].r - fq[li]*lines[li].x))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2)))
            )

        a51 = ((
            2*(pow(fq[li], 3)*lines[li].r - 2*pow(fp[li], 3)*lines[li].x)
            + 2*fp[li]*fq[li]*(fp[li]*lines[li].r - fq[li]*lines[li].x))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        a52 = ((
            2*pow(a[li], 2)*(fq[li]*pow(lines[li].r, 3)
                             - fp[li]*pow(lines[li].x, 3))
            - 4*a[li]*fp[li]*(pow(lines[li].r, 2) - pow(lines[li].x, 2)))
               / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                  - a[li]*fq[li]
                  * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # Complete a53 and a54

    return m, B_gn
