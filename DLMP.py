# DLMP.py
from gurobipy import Model, GRB


def calculateDLMP(
        dispatch, buses, generators, B_gn, lines, SMP, gensetP, gensetU):
    m = Model()

    # Create variables
    # Variable for voltage square
    v = m.addVars((b for b in buses), lb=0, name="v")
    # Variable for current square
    a = m.addVars((li for li in lines), lb=0, name="a")
    fp = m.addVars((li for li in lines), name="fp")
    fq = m.addVars((li for li in lines), name="fq")
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
    m.addConstrs(
        ((fp[li]*fp[li] + fq[li]*fq[li])
            <= lines[li].u*lines[li].u for li in lines),
        name="linecapfw"
    )
    m.addConstrs(
        ((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r)
            + (fq[li] - a[li]*lines[li].x)*(fq[li] - a[li]*lines[li].x)
            <= (lines[li].u)*(lines[li].u) for li in lines),
        name="linecapbw"
    )
    m.addConstrs(
        ((v[lines[li].fbus] - 2*(
                lines[li].r*fp[li] + lines[li].x*fq[li]
                + a[li]*(lines[li].r*lines[li].r + lines[li].x*lines[li].x)))
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
            m.addConstr(v[generators[g].location] == 1)

    m.addConstrs(
        ((
            sum(fp[li] for li in buses[b].outline)
            - sum(fp[li] - lines[li].r*a[li] for li in buses[b].inline)
            - sum(pg[g] for g in B_gn[b]) + buses[b].Pd + v[b]*buses[b].Gs
            ) == 0 for b in buses),
        name="pbalance"
    )
    m.addConstrs(
        ((
            sum(fq[li] for li in buses[b].outline)
            - sum(fp[li] - lines[li].x*a[li] for li in buses[b].inline)
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
    # m.addConstrs(
    #     (dispatch[g] == pg[g] for g in gensetP),
    #     name="genpower"
    # )

    m.update()

    m.optimize()

    return m
