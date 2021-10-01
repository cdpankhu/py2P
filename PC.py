# PC.py
import gurobipy
import pandas
from py2P.NetworkDataType import Bus, Line, Generator
from py2P.P2PTradeType import Trade, Agent
from py2P.NetworkLoad import networkload


testsystem = "AP15busDN"


def pc(testsystem):
    buses, lines, generators, datamat = networkload(testsystem)
    lineset = list(range(len(lines)))
    busset = list(range(len(buses)))
    genset = list(range(len(generators)))
    buses[0].Vmax = 1
    buses[0].Vmin = 1
    windset = []

    for g in genset:
        if generators[g].gtype == "Wind":
            windset.append(generators[g].gindex)

    root = -1
    for g in genset:
        if generators[g].gtype == "Root":
            root = g

    if root == -1:
        print("No generator at root node")
        generators[0].gtype = "Root"
        root = 0

    gensetU = root
    gensetP = genset
    del gensetP[root]

    # Defining the set of generator buses and demand buses
    B_g = []
    for g in genset:
        if not (generators[g].location in B_g):
            B_g.append(generators[g].location)
            buses[generators[g].location].Pd = 0
    B_d = [i for i in busset if i not in B_g]

    # Defining matrix of generators at each bus
    B_gn = []
    for b in busset:
        B_gn.append([])
    for g in genset:
        B_gn[generators[g].location-1].append(g)

    # Share of demand from utility versus peers
    partlevel = 0.0
    gam = partlevel*[1 for i in range(len(busset))]

    # Scale of each trade in MWh
    trade_scale = 1e-3

    # Defining trading agents at non-root buses
    agents = []
    for b in busset:
        if not B_gn[b]:
            temp_agent = Agent(
                b+1, b+1, 1, 0, 0, [], gam[b]*buses[b].Pd, gam[b].buses[b].Qd)
            agents.append(temp_agent)
        else:
            for g in B_gn[b]:
                if generators[g].gtype != "Root":
                    temp_agent = Agent(
                        b+1, b+1, 1, generators[g].Pmin, generators[g].Pmax,
                        generators[g].cost, gam[b]*buses[b].Pd,
                        gam[b]*buses[b].Qd)
                    agents.append(temp_agent)

    # Defining trades
