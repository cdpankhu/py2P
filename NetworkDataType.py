# NetworkDataType.py

class Bus:
    def __init__(
            self, bindex, btype, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV,
            bzone, Vmax, Vmin):
        self.bindex = bindex
        self.btpye = btype
        self.Pd = Pd
        self.Qd = Qd
        self.Gs = Gs
        self.Bs = Bs
        self.area = area
        self.Vm = Vm
        self.Va = Va
        self.baseKV = baseKV
        self.bzone = bzone
        self.Vmax = Vmax
        self.Vmin = Vmin
        self.children = []
        self.ancestor = []
        self.inline = []
        self.outline = []


class Line:
    def __init__(self, lindex, fbus, tbus, r, x, b, u):
        self.lindex = lindex
        self.fbus = fbus
        self.tbus = tbus
        self.r = r
        self.x = x
        self.b = b
        self.u = u
        self.children = []
        self.ancestor = []


class Generator:
    def __init__(
            self, gindex, gtype, location, Pg, Qg, Qmax, Qmin, Vg, mBase,
            status, Pmax, Pmin, cost, SUcost, SDcost, RU, RD, UPtime, DNtime):
        self.gindex = gindex
        self.gtype = gtype
        self.location = location
        self.Pg = Pg
        self.Qg = Qg
        self.Qmax = Qmax
        self.Qmin = Qmin
        self.Vg = Vg
        self.mBase = mBase
        self.status = status
        self.Pmax = Pmax
        self.Pmin = Pmin
        self.cost = cost
        self.SUCost = SUcost
        self.SDCost = SDcost
        self.RU = RU
        self.RD = RD
        self.UPtime = UPtime
        self.DNtime = DNtime
