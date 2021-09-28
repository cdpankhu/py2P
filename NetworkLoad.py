# NetworkLoad.py
import os
import pandas
import NetworkDataType
import pypower.api


def networkload(testsystem):
    if testsystem == "AP15busDN":
        usingMatPower = 0
    elif testsystem == "IEEE33busDN" or testsystem == "IEEE33busDN2" \
            or testsystem == "ISONE8busTN":
        usingMatPower = 0
    elif testsystem == "R2-25.00-1":
        usingMatPower = 0
    else:
        usingMatPower = 1
        PUcorrection = 100
        if testsystem == "case141_OPF":
            PUcorrection = 10

    if usingMatPower == 0:
        return networkload_nomatpower(testsystem)
    elif usingMatPower == 1:
        return networkload_matpower(testsystem, PUcorrection)


def networkload_nomatpower(testsystem):
    filename_Node = os.getcwd()+"/data/"+testsystem+"/Node.csv"
    filename_Generator = os.getcwd()+"/data"+testsystem+"/Generator.csv"
    filename_Line = os.getcwd()+"/data/"+testsystem+"/Line.csv"

    busmat = pandas.read_csv(filename_Node)
    buses = []
    for i in range(1, len(busmat)):
        bindex = busmat[i, 1]
        btype = 0
        Pd = busmat[i, 2]
        Qd = busmat[i, 3]
        Gs = busmat[i, 6]
        Bs = busmat[i, 7]
        area = 0
        Vm = 0
        Va = 0
        baseKV = busmat[i, 8]
        bzone = 0
        Vmax = busmat[i, 4]
        Vmin = busmat[i, 5]
        buses.append(NetworkDataType.Bus(
            bindex, btype, Pd, Qd, Gs, Bs,
            area, Vm, Va, baseKV, bzone, Vmax, Vmin))

    genmat = pandas.read_csv(filename_Generator)
    generators = []
    for i in range(1, len(genmat)):
        gindex = genmat[i, 1]
        gtype = genmat[i, 2]
        location = genmat[i, 3]
        Pg = genmat[i, 4]
        Qg = genmat[i, 5]
        Qmax = genmat[i, 6]
        Qmin = genmat[i, 7]
        Vg = genmat[i, 8]
        mBase = genmat[i, 9]
        status = genmat[i, 10]
        Pmax = genmat[i, 11]
        Pmin = genmat[i, 12]
        cost = [genmat[i, 13], genmat[i, 14], genmat[i, 15]]
        SUcost = genmat[i, 16]
        SDcost = genmat[i, 17]
        RU = genmat[i, 18]
        RD = genmat[i, 19]
        UPtime = genmat[i, 20]
        DNtime = genmat[i, 21]
        generators.append(NetworkDataType.Generator(
            gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status,
            Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime))

    branchmat = pandas.read_csv(filename_Line)
    lines = []
    for i in range(1, len(branchmat)):
        lindex = i
        fbus = int(branchmat[i, 2])
        tbus = int(branchmat[i, 1])
        r = branchmat[i, 3]
        x = branchmat[i, 4]
        b = branchmat[i, 5]
        u = branchmat[i, 6]
        buses[tbus].children.append(fbus)
        buses[fbus].ancestor.append(tbus)
        buses[tbus].inline.append(lindex)
        buses[fbus].outline.append(lindex)
        lines.append(NetworkDataType.Line(
            lindex, fbus, tbus, r, x, b, u))

    for j in range(1, len(lines)):
        for k in range(1, len(lines)):
            if lines[j].fbus == lines[k].tbus:
                lines[j].ancestor.append(k)
            elif lines[j].tbus == lines[k].fbus:
                lines[j].children.append(k)

    datamat = dict()
    datamat["bus"] = busmat
    datamat["branch"] = branchmat
    datamat["gen"] = genmat
    return buses, lines, generators, datamat


def networkload_matpower(testsystem, PUcorrection):
    mpc = matpowercase(testsystem)
    datamat = mpc

    buses = []
    for i in range(1, len(mpc["bus"])):
        bindex = mpc["bus"][i, 1]
        btype = mpc["bus"][i, 2]
        Pd = mpc["bus"][i, 3]
        Qd = mpc["bus"][i, 4]
        Gs = mpc["bus"][i, 5]
        Bs = mpc["bus"][i, 6]
        area = mpc["bus"][i, 7]
        Vm = mpc["bus"][i, 8]
        Va = mpc["bus"][i, 9]
        baseKV = mpc["bus"][i, 10]
        bzone = mpc["bus"][i, 11]
        Vmax = mpc["bus"][i, 12]
        Vmin = mpc["bus"][i, 13]
        buses.append(NetworkDataType.Bus(
            bindex, btype, Pd, Qd, Gs, Bs,
            area, Vm, Va, baseKV, bzone, Vmax, Vmin))

    generators = []
    for i in range(1, len(mpc["gen"])):
        gindex = i
        gtype = "NotDefined"
        location = mpc["gen"][i, 1]
        Pg = mpc["gen"][i, 2]
        Qg = mpc["gen"][i, 3]
        Qmax = mpc["gen"][i, 4]
        Qmin = mpc["gen"][i, 5]
        Vg = mpc["gen"][i, 6]
        mBase = mpc["gen"][i, 7]
        status = mpc["gen"][i, 8]
        Pmax = mpc["gen"][i, 9]
        Pmin = mpc["gen"][i, 10]
        cost = [mpc["gencost"][i, 5],
                mpc["gencost"][i, 6], mpc["gencost"][i, 7]]
        SUcost = mpc["gencost"][i, 2]
        SDcost = mpc["gencost"][i, 3]
        RU = 0
        RD = 0
        UPtime = 0
        DNtime = 0
        generators.append(NetworkDataType.Generator(
            gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status,
            Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime))

    lines = []
    for i in range(1, len(mpc["branch"])):
        lindex = i
        fbus = int(mpc["branch"][i, 2])
        tbus = int(mpc["branch"][i, 1])
        r = mpc["branch"][i, 3]/PUcorrection
        x = mpc["branch"][i, 4]/PUcorrection
        b = mpc["branch"][i, 5]
        u = 10000
        buses[tbus].children.append(fbus)
        buses[fbus].ancestor.append(tbus)
        buses[tbus].inline.append(lindex)
        buses[fbus].outline.append(lindex)
        lines.append(NetworkDataType.Line(
            lindex, fbus, tbus, r, x, b, u))

    if testsystem == "case141_OPF":
        for j in range(1, 6):
            lines[j].u = 50
        for j in range(7, len(lines)):
            lines[j].u = 5.0

    for j in range(1, len(lines)):
        for k in range(1, len(lines)):
            if lines[j].fbus == lines[k].tbus:
                lines[j].ancestor.append(k)
            elif lines[j].tbus == lines[k].fbus:
                lines[j].children.append(k)

    return buses, lines, generators, datamat


def matpowercase(testsystem):
    if testsystem == "case4gs":
        return pypower.api.case4gs()
    elif testsystem == "case6ww":
        return pypower.api.case6ww()
    elif testsystem == "case9":
        return pypower.api.case9()
    elif testsystem == "case9Q":
        return pypower.api.case9Q()
    elif testsystem == "case9target":
        return pypower.api.case9target()
    elif testsystem == "case14":
        return pypower.api.case14()
    elif testsystem == "case24_ieee_rts":
        return pypower.api.case24_ieee_rts()
    elif testsystem == "case30":
        return pypower.api.case30()
    elif testsystem == "case30pwl":
        return pypower.api.case30pwl()
    elif testsystem == "case30Q":
        return pypower.api.case30Q()
    elif testsystem == "case39":
        return pypower.api.case39()
    elif testsystem == "case57":
        return pypower.api.case57()
    elif testsystem == "case118":
        return pypower.api.case118()
    elif testsystem == "case300":
        return pypower.api.case300()
