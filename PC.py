# PC.py
import gurobipy
from pandas import DataFrame
from py2P.NetworkDataType import Bus, Line, Generator
from py2P.P2PTradeType import Trade, Agent
from py2P.NetworkLoad import networkload
from py2P.TradeFunction import generatetrade, selecttrade
from py2P.DLMP import calculatedlmp
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
    # trade_scale = 1  # trade in MWh

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
                b, b, 1, 0, 0, [0, 0, 0], gam[b]*buses[b].Pd,
                gam[b]*buses[b].Qd)
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
    trades_selected_mat = {}
    trades_selected_mat_buy = {}
    trades_selected_mat_sell = {}
    for i in buses:
        for j in buses:
            trades_selected_mat[i, j] = []
            trades_selected_mat_buy[i, j] = []
            trades_selected_mat_sell[i, j] = []
    for w in trades:
        trades_mat[trades[w].As, trades[w].Ab].append(trades[w].tradeID)

    trades_selected = []
    trades_rest = []
    trades_selected_old = []
    trades_rest_old = []
    dispatch = {}
    gencon = {}
    generation = {}
    consumption = {}
    gencon_old = {}
    dispatch_old = {}
    for i in buses:
        dispatch[i] = 0
        gencon[i] = 0
        generation[i] = 0
        consumption[i] = 0
        gencon_old[i] = 0
        dispatch_old[i] = 0
    status = 0
    dlmp = {}
    dlmp_stack = []
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
        print("Outer loop")
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
            print("Inner loop")
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

        # Determine nodal operating points according to incremental trades
        for w in trades_selected:
            dispatch[agents[trades[w].Ab].location] -= 1 * trade_scale
            dispatch[agents[trades[w].As].location] += 1 * trade_scale
            gencon[agents[trades[w].Ab].location] -= 1 * trade_scale
            gencon[agents[trades[w].As].location] += 1 * trade_scale
            consumption[agents[trades[w].Ab].location] -= 1 * trade_scale
            generation[agents[trades[w].As].location] += 1 * trade_scale

        # Set generator dispatch values
        for g in gensetP:
            dispatch_peerG[g] = generation[generators[g].location]

        # Calculate network charge
        SMP = generators[root].cost[1]
        print(generation, consumption)
        status, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo = \
            calculatedlmp(
                dispatch_peerG, buses, generators, lines, SMP, gensetP, gensetU
                )

        # Update network charge
        deltaNw = trade_scale
        if status == 2:
            for w in trades:
                if w in trades_selected:
                    trades[w].costNw = (
                        dlmp[agents[trades[w].Ab].location]
                        - dlmp[agents[trades[w].As].location])/2
        else:
            for w in trades_selected:
                trades[w].costNw += param_delta*deltaNw/trade_scale

        dispatch_stack.append(dispatch)
        dlmp_stack.append(dlmp)

        if gencon_old == gencon and status == 2:
            break
        else:
            for d in dispatch:
                dispatch_old[d] = dispatch[d]
            for g in gencon:
                gencon_old[g] = gencon[g]

    trades_dis = {}
    nwcharge_dis = {}
    echarge_dis = {}
    for i in buses:
        for j in buses:
            trades_dis[i, j] = 0
            nwcharge_dis[i, j] = 0
            echarge_dis[i, j] = 0
    temp = {}
    for b in buses:
        temp[b] = []

    for w in trades_selected:
        trades_dis[agents[trades[w].Ab].location,
                   agents[trades[w].As].location] -= trade_scale
        trades_dis[agents[trades[w].As].location,
                   agents[trades[w].Ab].location] += trade_scale
        echarge_dis[agents[trades[w].Ab].location,
                    agents[trades[w].As].location] -= trades[w].Peb
        echarge_dis[agents[trades[w].As].location,
                    agents[trades[w].Ab].location] += trades[w].Pes
        nwcharge_dis[agents[trades[w].Ab].location,
                     agents[trades[w].As].location] \
            -= round(trades[w].costNw, 4)
        nwcharge_dis[agents[trades[w].As].location,
                     agents[trades[w].Ab].location] \
            += round(trades[w].costNW, 4)
        temp[agents[trades[w].Ab].location].append(w)
        temp[agents[trades[w].As].location].append(w)

    price_mat = {}
    nprice_mat = {}
    pricemin = {}
    pricemax = {}
    nprice_mat_min = {}
    nprice_mat_max = {}
    eprice_mat = {}
    eprice_mat_min = {}
    eprice_mat_max = {}
    for i in buses:
        for j in buses:
            price_mat[i, j] = 0
            nprice_mat[i, j] = 0
            pricemin[i, j] = 0
            pricemax[i, j] = 0
            nprice_mat_min[i, j] = 0
            nprice_mat_max[i, j] = 0
            eprice_mat[i, j] = 0
            eprice_mat_min[i, j] = 0
            eprice_mat_max[i, j] = 0

    for w in trades_selected:
        trades_selected_mat[agents[trades[w].As].location,
                            agents[trades[w].Ab].location].append(
                            trades[w].tradeID)

    for i in buses:
        for j in buses:
            temp_mat = []
            for w in trades_selected_mat[i, j]:
                temp_mat.append(trades[w].Pes)
                pricemin[i, j] = min(temp_mat)
                pricemin[j, i] = -min(temp_mat)
                pricemax[i, j] = max(temp_mat)
                pricemax[j, i] = max(temp_mat)
                if pricemin[i, j] == pricemax[i, j]:
                    price_mat[i, j] = pricemin[i, j]
                    price_mat[j, i] = -pricemin[i, j]
                else:
                    print("error: pricemin != pricemax")

    for w in trades_selected:
        nprice_mat[agents[trades[w].Ab].location,
                   agents[trades[w].As].location] \
            = 0.5*round(trades[w].costNw, 5)
        nprice_mat[agents[trades[w].As].location,
                   agents[trades[w].Ab].location] \
            = 0.5*round(trades[w].costNw, 5)

    for w in trades:
        nprice_mat_min[agents[trades[w].Ab].location,
                       agents[trades[w].As].location] \
            = min(round(trades[w].costNW, 5),
                  nprice_mat_min[agents[trades[w].Ab].location,
                                 agents[trades[w].As].location])
        nprice_mat_min[agents[trades[w].As].location,
                       agents[trades[w].Ab].location] \
            = min(round(trades[w].costNw, 5),
                  nprice_mat_min[agents[trades[w].As].location,
                                 agents[trades[w].Ab].location])
        nprice_mat_max[agents[trades[w].Ab].location,
                       agents[trades[w].As].location] \
            = max(round(trades[w].costNW, 5),
                  nprice_mat_max[agents[trades[w].Ab].location,
                                 agents[trades[w].As].location])
        nprice_mat_max[agents[trades[w].As].location,
                       agents[trades[w].Ab].location] \
            = max(round(trades[w].costNw, 5),
                  nprice_mat_max[agents[trades[w].As].location,
                                 agents[trades[w].Ab].location])

    for w in trades_selected:
        eprice_mat[agents[trades[w].Ab].location,
                   agents[trades[w].As].location] = round(trades[w].Pes, 5)
        eprice_mat[agents[trades[w].As].location,
                   agents[trades[w].Ab].location] = round(trades[w].Peb, 5)

    for w in trades:
        eprice_mat_min[agents[trades[w].Ab].location,
                       agents[trades[w].As].location] \
            = min(round(trades[w].Pes, 5),
                  eprice_mat_min[agents[trades[w].Ab].location,
                                 agents[trades[w].As].location])
        eprice_mat_min[agents[trades[w].As].location,
                       agents[trades[w].Ab].location] \
            = min(round(trades[w].Pes, 5),
                  eprice_mat_min[agents[trades[w].As].location,
                                 agents[trades[w].Ab].location])
        eprice_mat_max[agents[trades[w].Ab].location,
                       agents[trades[w].As].location] \
            = max(round(trades[w].Pes, 5),
                  eprice_mat_max[agents[trades[w].Ab].location,
                                 agents[trades[w].As].location])
        eprice_mat_max[agents[trades[w].As].location,
                       agents[trades[w].Ab].location] \
            = max(round(trades[w].Pes, 5),
                  eprice_mat_max[agents[trades[w].As].location,
                                 agents[trades[w].Ab].location])

    peer_revenue = {}
    # Base code may multiply by 2, to confirm
    # I would think this will return a revenue of 0 since Pes = Peb
    for i in buses:
        peer_revenue[i] = 0
        for j in buses:
            peer_revenue[i] += echarge_dis[i, j]

    consumer_cost = 0
    generator_revenue = 0
    for b in B_d:
        consumer_cost -= peer_revenue[b]
    for b in B_g:
        generator_revenue += peer_revenue[b]

    return consumer_cost, generator_revenue
