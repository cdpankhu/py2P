# makeModel.py
# Module containing code to make various nodels for use
from gurobipy import Model, GRB


def makeModel(buses, generators, lines, gensetP, gensetU, **optional):

    # Defining the set of generator buses and demand buses
    B_g = []
    for g in generators:
        if not (generators[g].location in B_g):
            B_g.append(generators[g].location)
            # For generation bus, set demand to 0
            # This needs to be changed with updated model for storage, DR, etc.
            buses[generators[g].location].Pd = 0
    B_d = [i for i in buses if i not in B_g]

    B_gn = {}
    for i in buses:
        B_gn[i] = []
    for g in generators:
        B_gn[generators[g].location].append(g)

    sbase = generators[1].mBase

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

    # Add coordination variables if SC
    if "SC" in optional:
        # P2P variables
        p = m.addVars((b for b in buses), name='p', lb=float('-inf'))
        pnm = m.addVars((b for b in buses), (b for b in buses), name='pnm',
                        lb=float('-inf'))

    m.update()

    # Set objective functions
    obj = -oc
    m.setObjective(obj, GRB.MAXIMIZE)

    m.update()

    # Create constraints
    # Add additional coordination constraints if SC
    if "SC" in optional:
        # Constraints for bidirectional trading
        m.addConstrs((pnm[i, j] + pnm[j, i] == 0 for i in buses for j in buses),
                     name='bidirection')

        m.addConstrs(
            (pnm[i, j] >= 0 for i in B_g for j in buses), name='signg')
        m.addConstrs(
            (pnm[i, j] <= 0 for i in B_d for j in buses), name='signd')
        m.addConstrs((p[i] == sum(pnm[i, j] for j in buses) for i in buses),
                     name='nodalbalance')

        m.addConstrs((sum(pg[g] for g in B_gn[n]
                          if not (g in gensetU)) == p[n] for n in B_g),
                     name='genpower')
        m.addConstrs(((-optional['gam'][n]) * buses[n].Pd/sbase == p[n]
                      for n in B_d), name='demand')

    # Define operating cost
    m.addConstr(oc == sum(ocgen[g] for g in generators), name="oc")
    # Define linear generation cost
    m.addConstrs(
        (ocgen[g] == (generators[g].cost[1]*pg[g]*sbase + generators[g].cost[2])
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
            m.addConstr(v[generators[g].location] == pow(generators[g].Vg, 2),
                        name="rootvoltage")
    # Active power balance constraint
    m.addConstrs(
        ((
            sum(fp[li] for li in buses[b].outline)
            - sum(fp[li] - lines[li].r*a[li] for li in buses[b].inline)
            - sum(pg[g] for g in B_gn[b])
            + buses[b].Pd/sbase + v[b]*buses[b].Gs
            ) == 0 for b in buses),
        name="pbalance"
    )
    # Reactive power balance constraint with branch admittance
    m.addConstrs(
        ((
            sum(fq[li] - lines[li].b*v[b]/2 for li in buses[b].outline)
            - sum(fq[li] - lines[li].x*a[li]
                  + lines[li].b*v[b]/2 for li in buses[b].inline)
            - sum(qg[g] for g in B_gn[b])
            + buses[b].Qd/sbase - v[b]*buses[b].Bs
            ) == 0 for b in buses),
        name="qbalance"
    )
    if "NodeInfo" in optional:
        for b in buses:
            if optional['NodeInfo'][b] == 1:
                m.addConstr(v[b] <= float('inf'), name="vmax["+str(b)+"]")
            else:
                m.addConstr(v[b] <= buses[b].Vmax*buses[b].Vmax,
                            name="vmax["+str(b)+"]")
            if optional['NodeInfo'][b] == -1:
                m.addConstr(v[b] >= 0, name="vmin["+str(b)+"]")
            else:
                m.addConstr(v[b] >= buses[b].Vmin*buses[b].Vmin,
                            name="vmin["+str(b)+"]")
    else:
        m.addConstrs(
            (v[b] <= buses[b].Vmax*buses[b].Vmax for b in buses),
            name="vmax"
        )
        m.addConstrs(
            (v[b] >= buses[b].Vmin*buses[b].Vmin for b in buses),
            name="vmin"
        )
    m.addConstrs(
        (pg[g] >= generators[g].Pmin/sbase for g in generators),
        name="pgmin"
    )
    m.addConstrs(
        (pg[g] <= generators[g].Pmax/sbase for g in generators),
        name="pgmax"
    )
    m.addConstrs(
        (qg[g] >= generators[g].Qmin/sbase for g in generators),
        name="qgmin"
    )
    m.addConstrs(
        (qg[g] <= generators[g].Qmax/sbase for g in generators),
        name="qgmax"
    )

    if 'dispatch' in optional:
        m.addConstrs(
            (optional['dispatch'][g]/sbase == pg[g] for g in gensetP),
            name="genpower"
            )
    m.Params.QCPDual = 1
    m.Params.BarQCPConvTol = 1e-7
    m.update()

    return m
