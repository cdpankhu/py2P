# DLMP.py
from gurobipy import Model, GRB


def calculateDLMP(
        dispatch, buses, generators, B_gn, lines, SMP, gensetP, gensetU):
    m = Model()

    # Create variables
    # Variable for voltage square
    v = m.addVars(len(buses), lb=0, name="v")
    # Variable for current square
    a = m.addVars(len(lines), lb=0, name="a")
    fp = m.addVars(len(lines), name="fp")
    fq = m.addVars(len(lines), name="fq")
    pg = m.addVars(len(generators), name="pg")
    qg = m.addVars(len(generators), name="qg")
    oc = m.addVar(lb=0, name="oc")
    ocgen = m.addVars(len(generators), lb=0, name="ocgen")

    m.update()

    # Set objective functions
    obj = -oc
    m.setObjective(obj, GRB.MAXIMIZE)

    # Create constraints
    # Define operating cost
    m.addConstr(oc == sum(ocgen[g] for g in range(len(generators))), name="oc")
    # Define linear generation cost
    m.addConstrs(
        (ocgen[g-1] == (generators[g].cost[1]*pg[g-1] + generators[g].cost[2])
            for g in generators), name="gencost"
    )
    m.addConstrs(
        ((fp[li-1]*fp[li-1] + fq[li-1]*fq[li-1])
            <= lines[li].u*lines[li].u for li in lines),
        name="linecapfw"
    )
    m.addConstrs(
        ((fp[li-1] - a[li-1]*lines[li].r)*(fp[li-1] - a[li-1]*lines[li].r)
            + (fq[li-1] - a[li-1]*lines[li].x)*(fq[li-1] - a[li-1]*lines[li].x)
            <= (lines[li].u)*(lines[li].u) for li in lines),
        name="linecapbw"
    )
    m.addConstrs(
        ((v[lines[li].fbus-1] - 2*(
                lines[li].r*fp[li-1] + lines[li].x*fq[li-1]
                + a[li-1]*(lines[li].r*lines[li].r + lines[li].x*lines[li].x)))
            == v[lines[li].tbus-1] for li in lines),
        name="betweennodes"
    )
    m.addConstrs(
        ((fp[li-1]*fp[li-1] + fq[li-1]*fq[li-1])
            <= v[lines[li].fbus-1]*a[li-1] for li in lines),
        name="socp"
    )

    # Voltage contraint for root node
    for g in generators:
        if generators[g].type == "Root":
            m.addConstr(v[generators[g].location-1] == 1)

    m.addConstrs(
        ((
            sum(fp[li-1] for li in buses[b].outline)
            - sum(fp[li-1] - lines[li].r*a[li-1] for li in buses[b].inline)
            - sum(pg[g-1] for g in B_gn[b]) + buses[b].Pd + v[b-1]*buses[b].Gs
            ) == 0 for b in buses),
        name="pbalance"
    )
    m.addConstrs(
        ((
            sum(fq[li-1] for li in buses[b].outline)
            - sum(fp[li-1] - lines[li].x*a[li-1] for li in buses[b].inline)
            - sum(qg[g-1] for g in B_gn[b]) + buses[b].Qd - v[b-1]*buses[b].Bs
            ) == 0 for b in buses),
        name="qbalance"
    )
    m.addConstrs(
        (v[b-1] <= buses[b].Vmax*buses[b].Vmax for b in buses),
        name="vmax"
    )
    m.addConstrs(
        (v[b-1] >= buses[b].Vmin*buses[b].Vmin for b in buses),
        name="vmin"
    )
    m.addConstrs(
        (pg[g-1] >= generators[g].Pmin for g in generators),
        name="pgmin"
    )
    m.addConstrs(
        (pg[g-1] <= generators[g].Pmax for g in generators),
        name="pgmax"
    )
    m.addConstrs(
        (qg[g-1] >= generators[g].Qmin for g in generators),
        name="qgmin"
    )
    m.addConstrs(
        (qg[g-1] <= generators[g].Qmax for g in generators),
        name="qgmax"
    )
    m.addConstrs(
        (dispatch[g] == pg[g-1] for g in generators),
        name="genpower"
    )

    m.optimize()
