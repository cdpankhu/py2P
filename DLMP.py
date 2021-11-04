# DLMP.py
from gurobipy import Model, GRB
from math import pow, sqrt
from pandas import DataFrame


def calculatedlmp(
        dispatch, buses, generators, lines, SMP, gensetP, gensetU):

    # Consider reducing cyclomatic complexity by moving shadow price extraction
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

    status = m.Status
    if not (status == GRB.OPTIMAL):
        return status, {}, {}, {}, {}, {}, {}

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
        # ((fp[l]^2+fq[l]^2)*lines[l].x+a[l]*fq[l]*(lines[l].r^2-lines[l].x^2)
        # -2*a[l]*fp[l]*lines[l].r*lines[l].x)
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a1[li] = ((
            (pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
            + a[li]*fq[li]*(pow(lines[li].r, 2) - pow(lines[li].x, 2))
            - 2*a[li]*fp[li]*lines[li].r*lines[li].x)
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # ((fp[l]^2+fq[l]^2)*lines[l].r-a[l]*fp[l]*(lines[l].r^2+lines[l].x^2))
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a2[li] = ((
            (pow(fp[li], 2) + pow(fq[li], 2))*lines[li].r
            - a[li]*fp[li]*(pow(lines[li].r, 2) + pow(lines[li].x, 2)))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # (-(fp[l]^2+fq[l]^2)*lines[l].r+a[l]*fp[l]*(lines[l].r^2-lines[l].x^2)
        # +2*a[l]*fq[l]*lines[l].r*lines[l].x)
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a3[li] = ((
            -(pow(fp[li], 2) + pow(fq[li], 2))*lines[li].r
            + a[li]*fp[li]*(pow(lines[li].r, 2) - pow(lines[li].x, 2))
            + 2*a[li]*fq[li]*lines[li].r*lines[li].x)
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # (2*(fq[l]^3*lines[l].r-fp[l]^3*lines[l].x)
        # +2*fp[l]*fq[l]*(fp[l]*lines[l].r-fq[l]*lines[l].x))
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a4[li] = ((
            2*(pow(fq[li], 3)*lines[li].r - pow(fp[li], 3)*lines[li].x)
            + 2*fp[li]*fq[li]*(fp[li]*lines[li].r - fq[li]*lines[li].x))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # (2*(fq[l]^3*lines[l].r-2*fp[l]^3*lines[l].x)
        # +2*fp[l]*fq[l]*(fp[l]*lines[l].r-fq[l]*lines[l].x))
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a51 = ((
            2*(pow(fq[li], 3)*lines[li].r - 2*pow(fp[li], 3)*lines[li].x)
            + 2*fp[li]*fq[li]*(fp[li]*lines[li].r - fq[li]*lines[li].x))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # (2*a[l]^2*(fq[l]*lines[l].r^3-fp[l]*lines[l].x^3)
        # -4*a[l]*fp[l]*fq[l]*(lines[l].r^2-lines[l].x^2))
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a52 = ((
            2*pow(a[li], 2)*(fq[li]*pow(lines[li].r, 3)
                             - fp[li]*pow(lines[li].x, 3))
            - 4*a[li]*fp[li]*fq[li]
               * (pow(lines[li].r, 2) - pow(lines[li].x, 2)))
               / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                  - a[li]*fq[li]
                  * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # (4*a[l]*lines[l].r*lines[l].x*(fp[l]^2-fq[l]^2))
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a53 = ((
            4*a[li]*lines[li].r*lines[li].x*(pow(fp[li], 2) - pow(fq[li], 2)))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        # (-2*a[l]^2*lines[l].r*lines[l].x*(fp[l]*lines[l].r-fq[l]*lines[l].x))
        # /((fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
        a54 = ((
            -2*pow(a[li], 2)*lines[li].r*lines[li].x
            * (fp[li]*lines[li].r - fq[li]*lines[li].x))
            / ((pow(fp[li], 2) + pow(fq[li], 2))*lines[li].x
                - a[li]*fq[li]
                * (pow(lines[li].r, 2) + pow(lines[li].x, 2))))

        a5[li] = a51 + a52 + a53 + a54

    dlmp = {}
    dlmp[1] = dual_pbalance[1]
    eq40 = {}
    eq41 = {}
    eq42 = {}
    eq43 = {}
    eq44 = {}
    for li in lines:
        eq40[li] = 0
        eq41[li] = 0
        eq42[li] = 0
        eq43[li] = 0
        eq44[li] = 0

    for li in lines:
        # The real power price at the ancestor node
        eq40[lines[li].fbus] = a1[li]*dual_pbalance[lines[li].tbus]
        # The reactive power price at the current node
        eq41[lines[li].fbus] = a2[li]*dual_qbalance[lines[li].fbus]
        # The reactive power price at the ancestor node
        eq42[lines[li].fbus] = a3[li]*dual_qbalance[lines[li].tbus]
        # The contribution of the first complex power constraint
        eq43[lines[li].fbus] = a4[li]*dual_linecapfw[li]
        # The contribution of the second complex power constraint
        eq44[lines[li].fbus] = a5[li]*dual_linecapbw[li]
        dlmp[lines[li].fbus] = (eq40[lines[li].fbus] + eq41[lines[li].fbus]
                                + eq42[lines[li].fbus] + eq43[lines[li].fbus]
                                + eq44[lines[li].fbus])

    print("::Network Info::\n")
    for b in buses:
        pgb[b] = round(pgb[b], 3)
        qgb[b] = round(qgb[b], 3)
        v[b] = round(v[b], 3)
    NodeInfo = DataFrame([pgb, qgb, v], index=['pg', 'qg', 'v'])
    print("::NodeInfo::\n", NodeInfo, "\n\n")
    flow = {}
    for li in lines:
        flow[li] = round(sqrt(pow(fp[li], 2) + pow(fq[li], 2)), 3)
        a[li] = round(a[li], 3)
        fp[li] = round(fp[li], 3)
        fq[li] = round(fq[li], 3)
    LineInfo = DataFrame([flow, a, fp, fq], index=['flow', 'i', 'fp', 'fq'])
    print("::LineInfo::\n", LineInfo, "\n\n")

    dlmpp = {}
    for b in buses:
        dlmpp[b] = round(dlmp[b], 2)
        eq40[b] = round(eq40[b], 2)
        eq41[b] = round(eq41[b], 2)
        eq42[b] = round(eq42[b], 2)
        eq43[b] = round(eq43[b], 2)
        eq44[b] = round(eq44[b], 2)
    DLMPInfo = DataFrame([dlmp, eq40, eq41, eq42, eq43, eq44],
                         index=["\u03BB_i", "eq40", "eq41", "eq42", "eq43",
                                "eq44"])
    print("::DLMPInfo::\n", DLMPInfo, "\n\n")

    GenInfo = DataFrame([pg, qg, ocgen], index=['pg', 'qg', 'ocgen'])

    return status, dlmp, pg, NodeInfo, LineInfo, DLMPInfo, GenInfo
