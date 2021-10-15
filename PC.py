# PC.py
import gurobipy
from pandas import DataFrame
from py2P.NetworkDataType import Bus, Line, Generator
from py2P.P2PTradeType import Trade, Agent
from py2P.NetworkLoad import networkload
from py2P.TradeFunction import generatetrade, selecttrade
from math import ceil


testsystem = "AP15busDN"


def pc(testsystem):
    buses, lines, generators, datamat = networkload(testsystem)
    buses[1].Vmax = 1
    buses[1].Vmin = 1
    windset = {}
    windex = 1

    for g in generators:
        if generators[g].gtype == "Wind":
            windset[windex] = generators[g].gindex
            windex += 1

    root = -1
    for g in generators:
        if generators[g].gtype == "Root":
            root = generators[g].location

    if root == -1:
        print("No generator at root node")
        generators[1].gtype = "Root"
        root = 1

    gensetP = []
    gensetU = root
    for g in generators:
        if not g == root:
            gensetP.append(g)

    # Defining the set of generator buses and demand buses
    B_g = []
    for g in generators:
        if not (generators[g].location in B_g):
            B_g.append(generators[g].location)
            # For generation bus, set demand to 0
            # This needs to be changed with updated model for storage, DR, etc.
            buses[generators[g].location].Pd = 0
    B_d = [i for i in buses if i not in B_g]

    # Defining matrix of generators at each bus
    B_gn = {}
    for b in buses:
        B_gn[b] = []
    for g in generators:
        B_gn[generators[g].location].append(g)

    # Share of demand from utility versus peers
    partlevel = 1.0
    gam = {}
    for b in buses:
        gam[b] = partlevel

    # Scale of each trade in MWh
    trade_scale = 1e-3  # trade in kWh
    #trade_scale = 1  # trade in MWh

    # Defining trading agents at non-root buses
    # Current definition defines a single agent of type 1 at each demand bus.
    # Current definition defines duplicate generator agents for each gen at bus
    # May need to refactor with addition of individual assets at single bus
    # multiple assets at single bus could act as aggregation as well???
    agents = {}
    agentID = 1
    for b in buses:
        if not B_gn[b]:
            agents[agentID] = Agent(
                b, b, 1, 0, 0, [0, 0, 0], gam[b]*buses[b].Pd, gam[b]*buses[b].Qd)
            agentID += 1
        else:
            for g in B_gn[b]:
                if generators[g].gtype != "Root":
                    agents[agentID] = Agent(
                        b, b, 1, generators[g].Pmin, generators[g].Pmax,
                        generators[g].cost, gam[b]*buses[b].Pd,
                        gam[b]*buses[b].Qd)
                    agentID += 1

    # Defining trades and result matrices
    trades = generatetrade(agents, trade_scale)
    trades_old = generatetrade(agents, trade_scale)

    trades_mat = {}
    trades_mat_buy = {}
    trades_mat_sell = {}
    for i in agents:
        for j in agents:
            trades_mat[i, j] = []
            trades_mat_buy[i, j] = []
            trades_mat_sell[i, j] = []
    trade_selected_mat = {}
    trade_selected_mat_buy = {}
    trade_selected_mat_sell = {}
    for i in buses:
        for j in buses:
            trade_selected_mat[i, j] = []
            trade_selected_mat_buy[i, j] = []
            trade_selected_mat_sell[i, j] = []
    for w in trades:
        trades_mat[trades[w].As, trades[w].Ab].append(trades[w].tradeID)

    trades_selected = []
    trades_rest = []
    trades_selected_old = []
    trades_rest_old = []
    dispatch = {}
    GenCon = {}
    Generation = {}
    Consumption = {}
    GenCon_old = {}
    dispatch_old = {}
    for i in buses:
        dispatch[i] = 0
        GenCon[i] = 0
        Generation[i] = 0
        Consumption[i] = 0
        GenCon_old[i] = 0
        dispatch_old[i] = 0
    status = 0
    DLMP = {}
    DLMP_stack = []
    dispatch_stack = []
    DLMPInfo = DataFrame()
    NodeInfo = DataFrame()
    LineInfo = DataFrame()
    dispatch_peerG = {}
    for i in generators:
        dispatch_peerG[i] = 0
    param_delta = 1.0

    # need to add timer code for performance assessment
    while True:
        accepted = {}
        rejected = {}
        pg = {}
        Revenue_Sell = {}
        Cost_Buy = {}
        Cost_Network = {}
        Cost_DG = {}
        trades_old = trades
        trades_selected_old = trades_selected
        for w in trades:
            trades[w].Pes = 0
            trades[w].Peb = 0
        while True:
            accepted = {}
            rejected = {}
            pg = {}

            gap = 1e-3
            for agentID in agents:
                temp_accepted, temp_rejected, temp_pg, Revenue_Sell, Cost_Buy, \
                    Cost_Network, Cost_DG = \
                    selecttrade(trades, agents, agentID, gap, trade_scale)
                accepted[agentID] = temp_accepted
                rejected[agentID] = temp_rejected
                pg[agentID] = temp_pg

            # price delta = $5/MWh = $0.005/kWh
            deltaP = 5*trade_scale
            for w in trades:
                if (w in accepted[trades[w].Ab]) \
                        and not (w in accepted[trades[w].As]):
                    if trades[w].Peb > trades[w].Pes:
                        trades[w].Pes += deltaP
                    else:
                        trades[w].Peb += deltaP
            diff = -1
            if not trades:
                diff = 0
            else:
                diff = sum((w in accepted[trades[w].Ab]) and not (
                    w in accepted[trades[w].As]) for w in trades)

            if diff == 0:
                break

        trades_selected = []
        trades_rest = []
        for w in trades:
            if (w in accepted[trades[w].Ab]) and (w in accepted[trades[w].As]):
                trades_selected.append(w)
            else:
                trades_rest.append(w)

        # for i in agents:
        #     dispatch[i] = ceil(agents[i].Pd/trade_scale)  # MW to kW conversion

        # Determine nodal operating points according to incremental trades
        for w in trades_selected:
            dispatch[agents[trades[w].Ab].location] -= 1 * trade_scale
            dispatch[agents[trades[w].As].location] += 1 * trade_scale
            GenCon[agents[trades[w].Ab].location] -= 1 * trade_scale
            GenCon[agents[trades[w].As].location] += 1 * trade_scale
            Consumption[agents[trades[w].Ab].location] -= 1 * trade_scale
            Generation[agents[trades[w].As].location] += 1 * trade_scale

        # Set generator dispatch values
        for g in gensetP:
            dispatch_peerG[g] = Generation[generators[g].location]

        SMP = generators[root].cost[1]

        return dispatch_peerG, buses,\
            generators, B_gn, lines, SMP, gensetP, gensetU
