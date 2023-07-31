# DLMP.py
from gurobipy import Model, GRB
from math import pow, sqrt
from pandas import DataFrame
from py2P.Coefficients import anglerecovery
from py2P.makeModel import makeModel


def calculatedlmp(
        dispatch, buses, generators, lines, sBase, SMP, gensetP, gensetU, conic, **optional):

    print("dispatch:", dispatch)
    # Consider reducing cyclomatic complexity by moving shadow price extraction
    B_gn = {}
    for i in buses:
        B_gn[i] = []
    for g in generators:
        B_gn[generators[g].location].append(g)

    m = Model()
    if 'NodeInfo' in optional and 'LineInfo' in optional:
        m = makeModel(buses, generators, lines, sBase, gensetP, gensetU, conic,
                      dispatch=dispatch, NodeInfo=optional['NodeInfo'],
                      LineInfo=optional['LineInfo'])
    elif 'NodeInfo' in optional:
        m = makeModel(buses, generators, lines, sBase, gensetP, gensetU, conic,
                      dispatch=dispatch, NodeInfo=optional['NodeInfo'])
    elif 'LineInfo' in optional:
        m = makeModel(buses, generators, lines, sBase, gensetP, gensetU, conic,
                      dispatch=dispatch, LineInfo=optional['LineInfo'])
    else:
        m = makeModel(buses, generators, lines, sBase, gensetP,
                      gensetU, conic, dispatch=dispatch)

    m.optimize()

    status = m.Status
    print("status: ", status)

    # IIS is currently disabled
    # if (status == GRB.INFEASIBLE):
    #     nodestatus, linestatus = calculatebinding(
    #         m, len(buses), len(lines), len(generators))
    #     NodeInfo = DataFrame([nodestatus], index=['status'])
    #     LineInfo = DataFrame([linestatus], index=['status'])
    #     print(NodeInfo.loc['status'], LineInfo.loc['status'])
    #     return status, {}, {}, NodeInfo, LineInfo, {}, {}
    if (status == GRB.INFEASIBLE):
        nodestatus = {} 
        linestatus = {}
        for b in buses:
            nodestatus[b] = 0
        for li in lines:
            linestatus[li] = 0
        NodeInfo = DataFrame([nodestatus], index=['status'])
        LineInfo = DataFrame([linestatus], index=['status'])
        return status, m, {}, NodeInfo, LineInfo, {}, {}

    var = m.getVars()
    # Getting constraint duals
    constrs = m.getConstrs()
    qconstrs = m.getQConstrs()

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

    theta = anglerecovery(lines, buses, fp, fq, v)

    state_linecapfw = {}
    state_linecapbw = {}
    socp_lhs = {}
    socp_rhs = {}
    for li in lines:
        state_linecapfw[li] = int((fp[li]*fp[li] + fq[li]*fq[li]) > (lines[li].u*lines[li].u)/(sBase*sBase))
        state_linecapbw[li] = int((fp[li] - a[li]*lines[li].r)*(fp[li] - a[li]*lines[li].r) + (fq[li] - a[li]*lines[li].x) * (fq[li] - a[li]*lines[li].x) > (lines[li].u*lines[li].u)/(sBase*sBase))
        socp_lhs[li] = (fp[li]*fp[li] + fq[li]*fq[li])
        socp_rhs[li] = v[lines[li].fbus]*a[li]

    print("::Network Info::\n")
    for b in buses:
        pgb[b] = round(pgb[b]*sBase, 4)
        qgb[b] = round(qgb[b]*sBase, 4)
        v[b] = round(v[b], 4)
        theta[b] = round(theta[b], 5)
    NodeInfo = DataFrame([pgb, qgb, v, theta], index=[
                         'pg', 'qg', 'v', 'theta'])
    print("::NodeInfo::\n", NodeInfo, "\n\n")
    flow = {}
    a_real = {}
    fp_real = {}
    fq_real = {}
    for li in lines:
        flow[li] = round(sqrt(pow(fp[li]*sBase, 2) + pow(fq[li]*sBase, 2)), 3)
        a_real[li] = round(a[li], 3)
        fp_real[li] = round(fp[li]*sBase, 3)
        fq_real[li] = round(fq[li]*sBase, 3)
    LineInfo = DataFrame([flow, a_real, fp_real, fq_real, state_linecapfw, state_linecapbw, socp_lhs, socp_rhs], index=['flow', 'i', 'fp', 'fq', 'linecapfw', 'linecapbw', 'socp_lhs', 'socp_rhs'])
    print("::LineInfo::\n", LineInfo, "\n\n")

    for g in generators:
        pg[g] = pg[g]*sBase
        qg[g] = qg[g]*sBase

    GenInfo = DataFrame([pg, qg, ocgen], index=['pg', 'qg', 'ocgen'])

    if ('LineInfo' in optional) or ('NodeInfo' in optional) or conic != 1:
        return status, {}, pg, NodeInfo, LineInfo, {}, GenInfo
    

    # Getting quadratic constraint duals
    # Skip quadratic ocgen constraints
    qconstrCount = 0 + len(generators)
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
    # If ocgen is linear, need to skip len(generators)
    constrCount = 0 + 1
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
        #print("lines: ", li, lines[li].x, lines[li].r)
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

    dlmpp = {}
    for b in buses:
        dlmpp[b] = round(dlmp[b]/sBase, 2)
        eq40[b] = round(eq40[b]/sBase, 2)
        eq41[b] = round(eq41[b]/sBase, 2)
        eq42[b] = round(eq42[b]/sBase, 2)
        eq43[b] = round(eq43[b]/sBase, 2)
        eq44[b] = round(eq44[b]/sBase, 2)
    DLMPInfo = DataFrame([dlmpp, eq40, eq41, eq42, eq43, eq44],
                         index=["\u03BB_i", "eq40", "eq41", "eq42", "eq43",
                                "eq44"])
    print("::DLMPInfo::\n", DLMPInfo, "\n\n")

    return status, dlmpp, pg, NodeInfo, LineInfo, DLMPInfo, GenInfo


def calculatebinding(model, buscount, linecount, gencount):
    # constrs = model.getConstrs()
    # qconstrs = model.getQConstrs()
    # constrCount = 0
    # constrCount = 0 + 1
    # constrs[0].IISConstrForce = 0
    # for li in range(linecount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1

    # # Move count past root generator voltage constraint
    # constrs[constrCount].IISConstrForce = 1
    # constrCount += 1

    # for b in range(buscount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1
    # for b in range(buscount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1
    # #vmax
    # for b in range(buscount):
    #     constrCount += 1
    # #vmin
    # for b in range(buscount):
    #     constrCount += 1
    # for g in range(gencount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1
    # for g in range(gencount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1
    # for g in range(gencount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1
    # for g in range(gencount):
    #     constrs[constrCount].IISConstrForce = 1
    #     constrCount += 1
    
    # qconstrCount = 0
    # for g in range(gencount):
    #     qconstrs[qconstrCount].IISQConstrForce = 1
    #     qconstrCount += 1
    # #linecapfw
    # for li in range(linecount):
    #     qconstrCount += 1
    
    # #linecapbw
    # for li in range(linecount):
    #     qconstrCount += 1
    
    # for li in range(linecount):
    #     qconstrs[qconstrCount].IISQConstrForce = 1
    #     qconstrCount += 1

    model.computeIIS()
    iisconstr = model.IISConstr
    iisqconstr = model.IISQConstr

    # Getting base index for flow and voltage constraints
    linecapfwindex = 0 + gencount
    linecapbwindex = linecapfwindex + linecount
    vmaxindex = 0 + 1 + linecount + 1 + buscount + buscount
    vminindex = vmaxindex + buscount

    nodestatus = {}
    linestatus = {}

    # Assign nodestatus[b] == 1 if vmax violation, -1 if vmin violation
    for i in range(buscount):
        nodestatus[i+1] = 0
        if iisconstr[vmaxindex + i] == 1:
            nodestatus[i+1] = 1
        elif iisconstr[vminindex + i] == 1:
            nodestatus[i+1] = -1

    # Assign linestatus[li] == 1 for linecapfw violation, -1 if linecapbw
    for i in range(linecount):
        linestatus[i+1] = 0
        if iisqconstr[linecapfwindex + i] == 1:
            linestatus[i+1] = -1
        elif iisqconstr[linecapbwindex + i] == 1:
            linestatus[i+1] = 1

    return nodestatus, linestatus
