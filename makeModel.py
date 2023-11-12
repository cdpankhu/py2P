# makeModel.py
# Module containing code to make various nodels for use
from gurobipy import Model, GRB


def makeModel(buses, generators, lines, sBase, gensetP, gensetU, conic, **optional):

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

    m = Model()
    # Create variables
    # Variable for voltage square
    v = m.addVars((b for b in buses), lb=0, name="v")
    # Variable for current square
    a = m.addVars((li for li in lines), lb=0, name="a")
    fp = m.addVars((li for li in lines), name="fp", lb=float('-inf'))
    fq = m.addVars((li for li in lines), name="fq", lb=float('-inf'))
    pg = m.addVars((g for g in generators), name="pg", lb=float('-inf'))
    qg = m.addVars((g for g in generators), name="qg", lb=float('-inf'))
    oc = m.addVar(lb=0, name="oc")
    ocgen = m.addVars((g for g in generators), lb=float('-inf'), name="ocgen")

    # Add coordination variables if SC
    if "SC" in optional:
        # P2P variables
        # Total injection at bus b
        p = m.addVars((b for b in buses), name='p', lb=float('-inf'))
        # trade from bus i to bus j
        pnm = m.addVars((b for b in buses), (b for b in buses), name='pnm',
                        lb=float('-inf'))
        # trade from bus i to utility
        pnu = m.addVars((b for b in buses), name='pnm',
                        lb=float('-inf'))

    if "LineInfo" in optional:
        fslack = m.addVars((li for li in lines), name='fslack', lb = 0)
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
        # Injection at bus i is equal to the sum of transactions at bus i
        m.addConstrs((p[i] == sum(pnm[i, j] for j in buses) + pnu[i] for i in buses),
                     name='nodalbalance')
        # The sum of individual generators at bus n equal to total injection
        m.addConstrs((sum(pg[g] for g in B_gn[n]
                          if not (g in gensetU)) == p[n] for n in B_g),
                     name='genpower')
        # Demand buses consume specified demand
        m.addConstrs((-buses[n].Pd/sBase == p[n]
                      for n in B_d), name='demand')

    # Define operating cost
    pfactor = 10000
    if "LineInfo" in optional:
        m.addConstr(oc == sum(ocgen[g] for g in generators)+ sum(fslack[li]*fslack[li] for li in lines)*pfactor, name="oc")
    else:
        m.addConstr(oc == sum(ocgen[g] for g in generators), name="oc")
    # Define linear generation cost
    m.addConstrs(
        (ocgen[g] >= ((pg[g]*sBase)*(pg[g]*sBase)*generators[g].cost[0] + generators[g].cost[1]*pg[g]*sBase + generators[g].cost[2])
            for g in generators), name="gencost"
    )
    if "LineInfo" in optional:
        M = 1000
        for li in lines:
            # Apparent power flow limits on the receiving node
            if optional['LineInfo'][li] == -1 or optional['LineInfo'][li] == 2:
                # m.addConstr(((fp[li]*fp[li] + fq[li]*fq[li])
                #              <= M*(lines[li].u*lines[li].u)/(sBase*sBase)),
                #             name="linecapfw["+str(li)+"]")
                # Attempting to add slack variable with penalty for soft relaxed lines
                m.addConstr(((fp[li]*fp[li] + fq[li]*fq[li])
                             <= (lines[li].u*lines[li].u)/(sBase*sBase) + fslack[li]),
                            name="linecapfw["+str(li)+"]")
            else:
                m.addConstr(
                    ((fp[li]*fp[li] + fq[li]*fq[li])
                        <= (lines[li].u*lines[li].u)/(sBase*sBase)),
                    name="linecapfw["+str(li)+"]")
            # Apparent power flow limits on the sending node
            if optional['LineInfo'][li] == 1 or optional['LineInfo'][li] == 2:
                # m.addConstr(
                #     ((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r)
                #         + (fq[li] - a[li]*lines[li].x)
                #         * (fq[li] - a[li]*lines[li].x)
                #         <= M*(lines[li].u*lines[li].u)/(sBase*sBase)),
                #     name="linecapbw["+str(li)+"]")
                # Attempting to add slack variable with penalty for soft relaxed lines
                m.addConstr(
                    ((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r)
                        + (fq[li] - a[li]*lines[li].x)
                        * (fq[li] - a[li]*lines[li].x)
                        <= (lines[li].u*lines[li].u)/(sBase*sBase) + fslack[li]),
                    name="linecapbw["+str(li)+"]")
            else:
                m.addConstr(
                    ((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r)
                        + (fq[li] - a[li]*lines[li].x)
                        * (fq[li] - a[li]*lines[li].x)
                        <= (lines[li].u*lines[li].u)/(sBase*sBase)),
                    name="linecapbw["+str(li)+"]")
    else:
        # Apparent power flow limits on the receiving node
        m.addConstrs(
            ((fp[li]*fp[li] + fq[li]*fq[li])
                <= (lines[li].u*lines[li].u)/(sBase*sBase) for li in lines),
            name="linecapfw"
        )
        # Apparent power flow limits on the sending node
        m.addConstrs(
            ((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r)
                + (fq[li] - a[li]*lines[li].x)*(fq[li] - a[li]*lines[li].x)
                <= (lines[li].u*lines[li].u)/(sBase*sBase) for li in lines),
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
    if conic == 1:
        m.addConstrs(
            ((fp[li]*fp[li] + fq[li]*fq[li])
                <= v[lines[li].fbus]*a[li] for li in lines),
            name="socp"
        )
    else:
        m.addConstrs(
            ((fp[li]*fp[li] + fq[li]*fq[li])
                == v[lines[li].fbus]*a[li] for li in lines),
            name="socp"
        )

    # Voltage constraint for root node
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
            + buses[b].Pd/sBase + v[b]*buses[b].Gs
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
            + buses[b].Qd/sBase - v[b]*buses[b].Bs
            ) == 0 for b in buses),
        name="qbalance"
    )
    if "NodeInfo" in optional:
        #m.addConstrs(
        #    (v[b] <= float('inf') for b in buses),
        #    name="vmax"
        #)
        #m.addConstrs(
        #    (v[b] >= 0 for b in buses),
        #    name="vmin"
        #)
        for b in buses:
            if optional['NodeInfo'][b] == 1 or optional['NodeInfo'][b] == 2:
                m.addConstr(v[b] <= float('inf'), name="vmax["+str(b)+"]")
            else:
                m.addConstr(v[b] <= buses[b].Vmax*buses[b].Vmax,
                            name="vmax["+str(b)+"]")
        for b in buses:
            if optional['NodeInfo'][b] == -1 or optional['NodeInfo'][b] == 2:
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
        (pg[g] >= generators[g].Pmin/sBase for g in generators),
        name="pgmin"
    )
    m.addConstrs(
        (pg[g] <= generators[g].Pmax/sBase for g in generators),
        name="pgmax"
    )
    m.addConstrs(
        (qg[g] >= generators[g].Qmin/sBase for g in generators),
        name="qgmin"
    )
    m.addConstrs(
        (qg[g] <= generators[g].Qmax/sBase for g in generators),
        name="qgmax"
    )

    if 'dispatch' in optional:
        m.addConstrs(
            (optional['dispatch'][g]/sBase == pg[g] for g in gensetP),
            name="genpower"
            )

    # if 'SC' in optional:
    #     m.addConstr(pnm[12, 11] == 0.0217)
    #     m.addConstr(pnm[12, 10] == 0.0229)
    #     m.addConstr(pnm[12, 9] == 0.0235)
    #     m.addConstr(pnm[12, 8] == 0.0219)
    #     m.addConstr(pnm[12, 7] == 0.0219)
    #     m.addConstr(pnm[12, 6] == 0.0291)
    #     m.addConstr(pnm[12, 5] == 0.0173)
    #     m.addConstr(pnm[12, 4] == 0.0201)
    #     m.addConstr(pnm[12, 2] == 0.0756)
    #     m.addConstr(pnm[12, 2] == 0.254)

    m.Params.QCPDual = 1
    m.Params.BarQCPConvTol = 1e-7
    m.Params.BarHomogeneous = 1
    m.Params.DualReductions = 0
    m.Params.OptimalityTol = 1e-5
    # m.Params.InfUnbdInfo = 1
    # m.Params.Presolve = 0
    if conic != 1:
        m.Params.NonConvex = 2
    m.update()

    return m
