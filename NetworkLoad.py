# NetworkLoad.py
from os import getcwd
from pandas import read_csv
from py2P.NetworkDataType import Bus, Line, Generator
from pypower.api import \
    case118, case14, case24_ieee_rts, case30, case300, case30Q, case30pwl, \
    case39, case4gs, case57, case6ww, case9, case9Q, case9target
from py2P.case9DN import case9DN
from py2P.caseAP15busDN import caseAP15busDN
from py2P.caseAP15busDNtest import caseAP15busDNtest
from py2P.caseAP15busDN3gen import caseAP15busDN3gen
from py2P.case33prosumer import case33prosumer
from py2P.case141caracas import case141caracas


def networkload(testsystem):
    if testsystem == "AP15busDN" or testsystem == "AP15busDN_M":
        usingMatPower = 0
    elif testsystem == "IEEE33busDN" or testsystem == "IEEE33busDN2" \
            or testsystem == "ISONE8busTN":
        usingMatPower = 0
    elif testsystem == "R2-25.00-1":
        usingMatPower = 0
    elif testsystem == "Test2Bus":
        usingMatPower = 0
    elif testsystem == "6BusWiley":
        usingMatPower = 0
    else:
        usingMatPower = 1
        # PUcorrection = 100
        # if testsystem == "case141_OPF":
        #     PUcorrection = 10

    if usingMatPower == 0:
        return networkload_nomatpower(testsystem)
    elif usingMatPower == 1:
        # return networkload_matpower(testsystem, PUcorrection)
        return networkload_matpower(testsystem)


def networkload_nomatpower(testsystem):
    # filename_Node = getcwd()+"\\data\\"+testsystem+"\\Node.csv"
    # filename_Generator = getcwd()+"\\data\\"+testsystem+"\\Generator.csv"
    # filename_Line = getcwd()+"\\data\\"+testsystem+"\\Line.csv"
    filename_Node = \
        "N:\\Python\\py2P\\data\\"+testsystem+"\\Node.csv"
    filename_Generator = \
        "N:\\Python\\py2P\\data\\"+testsystem+"\\Generator.csv"
    filename_Line = \
        "N:\\Python\\py2P\\data\\"+testsystem+"\\Line.csv"

    sBase = 1

    busmat = read_csv(filename_Node)
    buses = {}
    for i in range(len(busmat)):
        bindex = busmat["Node"][i]
        btype = 0
        Pd = busmat["Pd"][i]
        Qd = busmat["Qd"][i]
        Gs = busmat["Gs"][i]
        Bs = busmat["Bs"][i]
        area = 0
        Vm = 0
        Va = 0
        baseKV = busmat["baseKV"][i]
        bzone = 0
        Vmax = busmat["Vmax"][i]
        Vmin = busmat["Vmin"][i]
        buses[bindex] = Bus(
            bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, bzone,
            Vmax, Vmin)
        # buses.append(Bus(
        #    bindex, btype, Pd, Qd, Gs, Bs,
        #    area, Vm, Va, baseKV, bzone, Vmax, Vmin))

    genmat = read_csv(filename_Generator)
    generators = {}
    for i in range(len(genmat)):
        gindex = genmat["gindex"][i]
        gtype = genmat["gtype"][i]
        location = genmat["location"][i]
        Pg = genmat["Pg"][i]
        Qg = genmat["Qg"][i]
        Qmax = genmat["Qmax"][i]
        Qmin = genmat["Qmin"][i]
        Vg = genmat["Vg"][i]
        mBase = genmat["mBase"][i]
        status = genmat["status"][i]
        Pmax = genmat["Pmax"][i]
        Pmin = genmat["Pmin"][i]
        cost = [
            genmat["cost-2nd"][i], genmat["cost-1st"][i], genmat["cost-nl"][i]]
        SUcost = genmat["SUcost"][i]
        SDcost = genmat["SDcost"][i]
        RU = genmat["RU"][i]
        RD = genmat["RD"][i]
        UPtime = genmat["UPtime"][i]
        DNtime = genmat["Dntime"][i]
        generators[gindex] = Generator(
            gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status,
            Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
        # generators.append(Generator(
        #    gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status,
        #    Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime))

    branchmat = read_csv(filename_Line)
    lines = {}
    for i in range(len(branchmat)):
        lindex = i+1
        # fbus = int(branchmat["to"][i])
        # tbus = int(branchmat["from"][i])
        # Definition of fbus is based on {Low, 2013} which is opposite matpower
        fbus = int(branchmat["from"][i])
        tbus = int(branchmat["to"][i])
        r = branchmat["r"][i]
        x = branchmat["x"][i]
        b = branchmat["b"][i]
        u = branchmat["s"][i]
        # buses[fbus].children.append(tbus)
        # buses[tbus].ancestor.append(fbus)
        # buses[tbus].inline.append(lindex)
        # buses[fbus].outline.append(lindex)
        buses[tbus].children.append(fbus)
        buses[fbus].ancestor.append(tbus)
        buses[tbus].inline.append(lindex)
        buses[fbus].outline.append(lindex)
        lines[lindex] = Line(
            lindex, fbus, tbus, r, x, b, u)
        # lines.append(Line(
        #    lindex, fbus, tbus, r, x, b, u))

    for j in lines:
        for k in lines:
            if lines[j].fbus == lines[k].tbus:
                lines[j].ancestor.append(k)
            elif lines[j].tbus == lines[k].fbus:
                lines[j].children.append(k)

    datamat = dict()
    datamat["bus"] = busmat
    datamat["branch"] = branchmat
    datamat["gen"] = genmat
    return buses, lines, generators, sBase, datamat, {}


def networkload_matpower(testsystem):
    mpc = matpowercase(testsystem)
    datamat = mpc

    sBase = mpc["baseMVA"]

    buses = {}
    for i in range(len(mpc["bus"])):
        bindex = mpc["bus"][i, 0]
        btype = mpc["bus"][i, 1]
        Pd = mpc["bus"][i, 2]
        Qd = mpc["bus"][i, 3]
        Gs = mpc["bus"][i, 4]
        Bs = mpc["bus"][i, 5]
        area = mpc["bus"][i, 6]
        Vm = mpc["bus"][i, 7]
        Va = mpc["bus"][i, 8]
        baseKV = mpc["bus"][i, 9]
        bzone = mpc["bus"][i, 10]
        Vmax = mpc["bus"][i, 11]
        Vmin = mpc["bus"][i, 12]
        buses[bindex] = Bus(
            bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, bzone,
            Vmax, Vmin)
        # buses.append(Bus(
        #    bindex, btype, Pd, Qd, Gs, Bs,
        #    area, Vm, Va, baseKV, bzone, Vmax, Vmin))

    generators = {}
    for i in range(len(mpc["gen"])):
        gindex = i+1
        gtype = "NotDefined"
        location = mpc["gen"][i, 0]
        Pg = mpc["gen"][i, 1]
        Qg = mpc["gen"][i, 2]
        Qmax = mpc["gen"][i, 3]
        Qmin = mpc["gen"][i, 4]
        Vg = mpc["gen"][i, 5]
        mBase = mpc["gen"][i, 6]
        status = mpc["gen"][i, 7]
        Pmax = mpc["gen"][i, 8]
        Pmin = mpc["gen"][i, 9]
        cost = [mpc["gencost"][i, 4],
                mpc["gencost"][i, 5], mpc["gencost"][i, 6]]
        SUcost = mpc["gencost"][i, 1]
        SDcost = mpc["gencost"][i, 2]
        RU = 0
        RD = 0
        UPtime = 0
        DNtime = 0
        generators[gindex] = Generator(
            gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status,
            Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime)
        # generators.append(Generator(
        #    gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase, status,
        #    Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime))

    lines = {}
    for i in range(len(mpc["branch"])):
        lindex = i+1
        # fbus = int(mpc["branch"][i, 0])
        # tbus = int(mpc["branch"][i, 1])
        # Definition of fbus is based on {Low, 2013} which is opposite matpower
        fbus = int(mpc["branch"][i, 1])
        tbus = int(mpc["branch"][i, 0])
        # r = mpc["branch"][i, 2]/PUcorrection
        # x = mpc["branch"][i, 3]/PUcorrection
        r = mpc["branch"][i, 2]
        x = mpc["branch"][i, 3]
        b = mpc["branch"][i, 4]
        u = mpc["branch"][i, 5]
        # For limit of 0, set limit high since no limit considered
        if u == 0:
            u = 1000
        # buses[fbus].children.append(tbus)
        # buses[tbus].ancestor.append(fbus)
        # buses[tbus].inline.append(lindex)
        # buses[fbus].outline.append(lindex)
        buses[tbus].children.append(fbus)
        buses[fbus].ancestor.append(tbus)
        buses[tbus].inline.append(lindex)
        buses[fbus].outline.append(lindex)
        lines[lindex] = Line(
            lindex, fbus, tbus, r, x, b, u)

    if testsystem == "case141_OPF":
        for j in lines:
            if lines[j].lindex < 7:
                lines[j].u = 50
            else:
                lines[j].u = 5.0

    for j in lines:
        for k in lines:
            # if lines[j].fbus == lines[k].tbus:
            #     lines[j].ancestor.append(k)
            # elif lines[j].tbus == lines[k].fbus:
            #     lines[j].children.append(k)
            if lines[j].fbus == lines[k].tbus:
                lines[j].ancestor.append(k)
            elif lines[j].tbus == lines[k].fbus:
                lines[j].children.append(k)

    return buses, lines, generators, sBase, datamat, mpc


def matpowercase(testsystem):
    if testsystem == "case4gs":
        return case4gs()
    elif testsystem == "case6ww":
        return case6ww()
    elif testsystem == "case9":
        return case9()
    elif testsystem == "case9Q":
        return case9Q()
    elif testsystem == "case9target":
        return case9target()
    elif testsystem == "case14":
        return case14()
    elif testsystem == "case24_ieee_rts":
        return case24_ieee_rts()
    elif testsystem == "case30":
        return case30()
    elif testsystem == "case30pwl":
        return case30pwl()
    elif testsystem == "case30Q":
        return case30Q()
    elif testsystem == "case39":
        return case39()
    elif testsystem == "case57":
        return case57()
    elif testsystem == "case118":
        return case118()
    elif testsystem == "case300":
        return case300()
    elif testsystem == "case9DN":
        return case9DN()
    elif testsystem == "caseAP15busDN":
        return caseAP15busDN()
    elif testsystem == "caseAP15busDNtest":
        return caseAP15busDNtest()
    elif testsystem == "caseAP15busDN3gen":
        return caseAP15busDN3gen()
    elif testsystem == "case33prosumer":
        return case33prosumer()
    elif testsystem == "case141caracas":
        return case141caracas()
