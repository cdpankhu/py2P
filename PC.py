# PC.py
from pandas import DataFrame
from py2P.P2PTradeType import Agent
from py2P.NetworkLoad import networkload
from py2P.TradeFunction import generatetrade, selecttrade
from py2P.DLMP import calculatedlmp
from py2P.Coefficients import calculateptdf
from py2P.makeJac import makeJacVSC
from py2P.penaltyutil import vm_p_sum
from py2P.SC import sc
from math import pow, sqrt
from time import strftime
from os import makedirs


def pc(testsystem, trade_scale, deltaP, pen_scale, pen_start, clearing):
    buses, lines, generators, sBase, datamat, ppc = networkload(testsystem)
    # buses[1].Vmax = 1
    # buses[1].Vmin = 1
    windset = {}
    windex = 1

    # get reference dlmp for nuc of tradings
    pnm, dlmp_reference = sc(testsystem)

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
    gensetU = [root]
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

    # # Scale of each trade in MWh
    # trade_scale = 1e-3  # trade in kWh
    # # trade_scale = 1  # trade in MWh

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

    ptdf = calculateptdf(lines, buses, root)

    # Defining trades and result matrices
    trades = generatetrade(agents, trade_scale)
    trades_old = generatetrade(agents, trade_scale)

    #Set trade nuc based on reference system usage charges
    for w in trades:
        trades[w].costNw = (
             dlmp_reference[agents[trades[w].Ab].location]
             - dlmp_reference[agents[trades[w].As].location])/2
        if agents[trades[w].Ab].location == 2 and agents[trades[w].As].location == 1:
            trades[w].costNw = 1/2
        elif agents[trades[w].Ab].location == 4 and agents[trades[w].As].location == 1:
            trades[w].costNw = 3/2
        elif agents[trades[w].Ab].location == 5 and agents[trades[w].As].location == 1:
            trades[w].costNw = 4/2
        elif agents[trades[w].Ab].location == 6 and agents[trades[w].As].location == 1:
            trades[w].costNw = 5/2
        elif agents[trades[w].Ab].location == 7 and agents[trades[w].As].location == 1:
            trades[w].costNw = 6/2
        elif agents[trades[w].Ab].location == 8 and agents[trades[w].As].location == 1:
            trades[w].costNw = 4/2
        elif agents[trades[w].Ab].location == 9 and agents[trades[w].As].location == 1:
            trades[w].costNw = 5/2
        elif agents[trades[w].Ab].location == 10 and agents[trades[w].As].location == 1:
            trades[w].costNw = 5/2
        elif agents[trades[w].Ab].location == 11 and agents[trades[w].As].location == 1:
            trades[w].costNw = 6/2
        elif agents[trades[w].Ab].location == 13 and agents[trades[w].As].location == 1:
            trades[w].costNw = 1/2
        elif agents[trades[w].Ab].location == 14 and agents[trades[w].As].location == 1:
            trades[w].costNw = 2/2
        elif agents[trades[w].Ab].location == 15 and agents[trades[w].As].location == 1:
            trades[w].costNw = 3/2
        elif agents[trades[w].Ab].location == 2 and agents[trades[w].As].location == 12:
            trades[w].costNw = 6/2
        elif agents[trades[w].Ab].location == 4 and agents[trades[w].As].location == 12:
            trades[w].costNw = 4/2
        elif agents[trades[w].Ab].location == 5 and agents[trades[w].As].location == 12:
            trades[w].costNw = 5/2
        elif agents[trades[w].Ab].location == 6 and agents[trades[w].As].location == 12:
            trades[w].costNw = 6/2
        elif agents[trades[w].Ab].location == 7 and agents[trades[w].As].location == 12:
            trades[w].costNw = 7/2
        elif agents[trades[w].Ab].location == 8 and agents[trades[w].As].location == 12:
            trades[w].costNw = 3/2
        elif agents[trades[w].Ab].location == 9 and agents[trades[w].As].location == 12:
            trades[w].costNw = 4/2
        elif agents[trades[w].Ab].location == 10 and agents[trades[w].As].location == 12:
            trades[w].costNw = 2/2
        elif agents[trades[w].Ab].location == 11 and agents[trades[w].As].location == 12:
            trades[w].costNw = 1/2
        elif agents[trades[w].Ab].location == 13 and agents[trades[w].As].location == 12:
            trades[w].costNw = 8/2
        elif agents[trades[w].Ab].location == 14 and agents[trades[w].As].location == 12:
            trades[w].costNw = 9/2
        elif agents[trades[w].Ab].location == 15 and agents[trades[w].As].location == 12:
            trades[w].costNw = 10/2

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
    # If more than one peer per bus, agent>buses
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
    for i in agents:
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
    pnm_stack = []
    DLMPInfo = DataFrame()
    NodeInfo = DataFrame()
    LineInfo = DataFrame()
    GenInfo = DataFrame()
    dispatch_peerG = {}
    max_penalty_delta = 0
    max_nuc_delta = 0
    min_trade_price = 0
    penalty_old = {}
    for w in trades:
        penalty_old[w] = 0
    for i in generators:
        dispatch_peerG[i] = 0
    penalty_iter = 0
    outer_iter = 0

    # need to add timer code for performance assessment
    # Outer loop of iterative process
    while True:
        # print("Outer loop")
        outer_iter += 1
        accepted = {}
        rejected = {}
        pg = {}
        Revenue_Sell = {}
        Cost_Buy = {}
        Cost_Network = {}
        Cost_Penalty = {}
        Cost_DG = {}
        trades_old = trades
        trades_selected_old = trades_selected
        trade_prices = []
        trade_penalty_deltas = []
        for w in trades:
            trade_prices.append(trades[w].Pes)
            trade_penalty_deltas.append(
                abs(trades[w].penalty - penalty_old[w]))
        min_trade_price = min(trade_prices)
        max_penalty_delta = max(trade_penalty_deltas)
        # set trade prices as the minimum final price minus max delta
        for w in trades:
            # trades[w].Pes = min_trade_price - max_penalty_delta - max_nuc_delta
            # trades[w].Peb = min_trade_price - max_penalty_delta - max_nuc_delta
            trades[w].Pes = 0
            trades[w].Peb = 0
            if trades[w].Pes < 0:
                trades[w].Pes = 0
            if trades[w].Peb < 0:
                trades[w].Peb = 0
            penalty_old[w] = trades[w].penalty
            trades[w].cleared = 0
        max_penalty_delta = 0
        max_nuc_delta = 0
        # Inner loop for negotiation of trades
        while True:
            # print("Inner loop")
            accepted = {}
            rejected = {}
            pg = {}

            gap = 1e-3
            for agentID in agents:
                print("agent: " + str(agentID) + " iter: " + str(outer_iter))
                status, temp_accepted, temp_rejected, temp_pg, Revenue_Sell, \
                    Cost_Buy, Cost_Network, Cost_Penalty, Cost_DG, model = \
                    selecttrade(trades, agents, agentID, gap, trade_scale)
                accepted[agentID] = temp_accepted
                rejected[agentID] = temp_rejected
                pg[agentID] = temp_pg

                if status != 2:
                    print("error with select trade")
                    return trades, agents, model, agentID, trade_scale

            # price delta = $5/MWh = $0.005/kWh
            # deltaP = 1
            # Update trade prices according to presence in optimal sets
            # Note price only increases and so may overshoot optimal
            changestate = 0
            for w in trades:
                if (w in accepted[trades[w].Ab]) \
                        and not (w in accepted[trades[w].As]):
                    changestate = 1
                    if trades[w].Peb > trades[w].Pes:
                        trades[w].Pes += deltaP
                    else:
                        trades[w].Peb += deltaP
                # Clear trades as agreed.
                if (w in accepted[trades[w].Ab]) \
                        and (w in accepted[trades[w].As]) and clearing == 1:
                    trades[w].cleared = 1
                # print(trades[w].tradeID, trades[w].As, trades[w].Ab,
                #       trades[w].Pes, trades[w].Peb, trades[w].costNw,
                #       trades[w].penalty,
                #       trades[w].Peb+trades[w].penalty+trades[w].costNw,
                #       w in accepted[trades[w].As], w in accepted[trades[w].Ab],
                #       trades[w].cleared)

            diff = -1
            if not trades:
                diff = 0
            else:
                diff = sum((w in accepted[trades[w].Ab]) and not (
                    w in accepted[trades[w].As]) for w in trades)

            if changestate == 0:
                break

        trades_selected = []
        trades_rest = []
        for w in trades:
            if (w in accepted[trades[w].Ab]) and (w in accepted[trades[w].As]):
                trades_selected.append(w)
            else:
                trades_rest.append(w)

        for i in agents:
            dispatch[i] = 0
            gencon[i] = 0
            generation[i] = 0
            consumption[i] = 0

        # Determine nodal operating points according to incremental trades
        for w in trades_selected:
            dispatch[agents[trades[w].Ab].location] -= 1 * trade_scale
            dispatch[agents[trades[w].As].location] += 1 * trade_scale
            gencon[agents[trades[w].Ab].location] -= 1 * trade_scale
            gencon[agents[trades[w].As].location] += 1 * trade_scale
            consumption[agents[trades[w].Ab].location] -= 1 * trade_scale
            generation[agents[trades[w].As].location] += 1 * trade_scale

        for a in agents:
            dispatch[a] = round(dispatch[a], 3)
            gencon[a] = round(gencon[a], 3)
            consumption[a] = round(consumption[a], 3)
            generation[a] = round(generation[a], 3)

        # Set generator dispatch values
        # This only works with one P generator per bus
        # Need to use agent object to work properly
        for g in gensetP:
            dispatch_peerG[g] = generation[generators[g].location]

        # Calculate network charge
        SMP = generators[root].cost[1]
        print(dispatch_peerG)
        status, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
            calculatedlmp(
                dispatch_peerG, buses, generators, lines, sBase, SMP, gensetP, gensetU
                )

        # Update network charge
        # Trade state is feasible and DLMPs become network charges
        if status == 2:
            # For all trades, update network charge
            for w in trades:
                costNwold = trades[w].costNw
                # trades[w].costNw = (
                #     dlmp[agents[trades[w].Ab].location]
                #     - dlmp[agents[trades[w].As].location])/2
                costNwdelta = abs(costNwold - trades[w].costNw)
                # track largest abs nuc change
                if costNwdelta > max_nuc_delta:
                    max_nuc_delta = costNwdelta
        # Trade state infeasible, update network charge
        else:
            penalty_iter += 1
            NodeState = NodeInfo.loc['status']
            LineState = LineInfo.loc['status']
            print(NodeState, LineState)
            nodestatus = sum(abs(NodeState[b]) for b in buses)
            linestatus = sum(abs(LineState[li]) for li in lines)
            # Run relaxed power flow to get feasible voltages and DLMPs
            print("Relaxed model")
            status2, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
                calculatedlmp(dispatch_peerG, buses, generators, lines, sBase,
                              SMP, gensetP, gensetU, NodeInfo=NodeState,
                              LineInfo=LineState)
            # Update nuc with relaxed DLMPs
            # if status2 == 2:
            #     for w in trades:
            #         trades[w].costNw = (
            #             dlmp[agents[trades[w].Ab].location]
            #             - dlmp[agents[trades[w].As].location])/2
            # Update penalty for ALL trades according to proportional
            # contribution to violated flow limits via ptdfs
            for li in lines:
                if LineState[li] != 0:
                    fwptdfsum = sum(abs(round(ptdf[li,
                                                   agents[trades[w].As].location,
                                                   agents[trades[w].Ab].location],
                                              4))
                                    for w in trades if
                                    round(ptdf[li, agents[trades[w].As].location,
                                               agents[trades[w].Ab].location],
                                          4) > 0)
                    bwptdfsum = sum(abs(round(ptdf[li,
                                                   agents[trades[w].As].location,
                                                   agents[trades[w].Ab].location],
                                              4))
                                    for w in trades if
                                    round(ptdf[li, agents[trades[w].As].location,
                                               agents[trades[w].Ab].location],
                                          4) < 0)
                    #  penpool is defined as the some of ptdf for bad trades
                    penpool = 0
                    if LineState[li] == -1:
                        penpool = bwptdfsum/pen_scale
                    elif LineState[li] == 1:
                        penpool = fwptdfsum/pen_scale
                    # Testing code to set penalty at advanced state
                    # Remove for final code/model
                    if penalty_iter == 1:
                        for w in trades:
                            wptdf = round(ptdf[li, agents[trades[w].As].location,
                                               agents[trades[w].Ab].location], 4)
                            if wptdf > 0:
                                trades[w].penalty += pen_start
                            elif wptdf < 0:
                                trades[w].penalty += 0
                    # end of code to remove

                    for w in trades:
                        wptdf = round(ptdf[li, agents[trades[w].As].location,
                                           agents[trades[w].Ab].location], 4)

                        if wptdf > 0:
                            trades[w].penalty += (LineState[li] * penpool
                                                  * wptdf / fwptdfsum)
                        elif wptdf < 0:
                            trades[w].penalty += (LineState[li] * penpool
                                                  * wptdf / bwptdfsum)
                    for w in trades:
                        print(w, li, LineState[li], agents[trades[w].As].location,
                              agents[trades[w].Ab].location, trades[w].Pes,
                              trades[w].Peb, trades[w].costNw,
                              trades[w].penalty,
                              trades[w].costNw+trades[w].penalty,
                              w in trades_selected)
            if nodestatus != 0:
                vm_p_sc, vm_q_sc, va_p_sc, va_q_sc = makeJacVSC(
                    ppc, NodeInfo.loc['v'], NodeInfo.loc['theta'])
                for b in buses:
                    if NodeState[b] != 0:
                        vsc_pos_sum, vsc_neg_sum = vm_p_sum(
                            b, trades, agents, vm_p_sc)
                        penpool = 0
                        if NodeState[b] == 1:
                            penpool = vsc_pos_sum
                        elif NodeState[b] == -1:
                            penpool = vsc_neg_sum
                        for w in trades:
                            wvsc = round(
                                vm_p_sc[b, agents[trades[w].As].location]
                                - vm_p_sc[b, agents[trades[w].Ab].location], 4)
                            if wvsc > 0:
                                trades[w].penalty += NodeState[b] * \
                                    wvsc * penpool / vsc_pos_sum
                            elif wvsc < 0:
                                trades[w].penalty += NodeState[b] * \
                                    wvsc * penpool / vsc_neg_sum
        dispatch_temp = {}
        dlmp_temp = {}
        pnm_temp = {}
        for i in buses:
            for j in buses:
                pnm_temp[i, j] = 0
        for w in trades_selected:
            pnm_temp[agents[trades[w].Ab].location,
                     agents[trades[w].As].location] -= trade_scale
            pnm_temp[agents[trades[w].As].location,
                     agents[trades[w].Ab].location] += trade_scale
        for i in dispatch:
            dispatch_temp[i] = dispatch[i]
        for i in dlmp:
            dlmp_temp[i] = dlmp[i]
        dispatch_stack.append(dispatch_temp)
        dlmp_stack.append(dlmp_temp)
        pnm_stack.append(pnm_temp)

        # If gencon is unchanged, break while loop
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
            += round(trades[w].costNw, 4)
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
            = min(round(trades[w].costNw, 5),
                  nprice_mat_min[agents[trades[w].Ab].location,
                                 agents[trades[w].As].location])
        nprice_mat_min[agents[trades[w].As].location,
                       agents[trades[w].Ab].location] \
            = min(round(trades[w].costNw, 5),
                  nprice_mat_min[agents[trades[w].As].location,
                                 agents[trades[w].Ab].location])
        nprice_mat_max[agents[trades[w].Ab].location,
                       agents[trades[w].As].location] \
            = max(round(trades[w].costNw, 5),
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

    line_loading = {}
    loading = {}
    flow = {}
    for li in lines:
        line_loading[li] = 1 - round(
            (lines[li].u - sqrt(
                (pow(LineInfo[li]['fp'], 2)
                 + pow(LineInfo[li]['fq'], 2))))/lines[li].u, 3)
        loading[li] = LineInfo[li]['flow'] / lines[li].u
        flow[li] = LineInfo[li]['flow']

    v_marginrate = {}
    voltage = {}
    for b in buses:
        v_marginrate[b] = round(0.1 - sqrt(NodeInfo[b]['v']), 3)/0.1
        voltage[b] = sqrt(NodeInfo[b]['v'])

    bus_frame = DataFrame([line_loading, v_marginrate, dlmp],
                          index=['line_loading', 'v_marginrate', 'dlmp'])
    # loading_frame = DataFrame([line_loading], index='line_loading')
    # vmargin_frame = DataFrame([v_marginrate], index='v_marginrate')

    demand_P2P = {}
    ep_u = {}
    ep_p = {}
    nuc = {}
    pen = {}
    for b in buses:
        demand_P2P[b] = 0
        ep_u[b] = 0
        ep_p[b] = 0
        nuc[b] = 0
        pen[b] = 0

    for m in B_d:
        demand_P2P[m] = sum(trades_dis[n, m] for n in B_g)
        ep_u[m] = (buses[m].Pd - demand_P2P[m])*dlmp[m]*trade_scale

    for w in trades_selected:
        ep_p[agents[trades[w].Ab].location] += trades[w].Pes*trade_scale
        nuc[agents[trades[w].Ab].location] += 0.5*trades[w].costNw*trade_scale
        nuc[agents[trades[w].As].location] += 0.5*trades[w].costNw*trade_scale
        pen[agents[trades[w].Ab].location] += trades[w].penalty*trade_scale

    ur = {}
    cp = {}
    for b in buses:
        ur[b] = ep_u[b] + nuc[b]
        cp[b] = ur[b] + ep_p[b] + pen[b]

    gc = {}
    for g in generators:
        gc[g] = generators[g].cost[1] * round(pg[g], 5)

    sum_gc_p = sum(gc[g] for g in gensetP)
    sum_gc_u = sum(gc[g] for g in gensetU)
    sum_dgr = sum(ep_p[b] for b in buses)
    sum_dgp = sum_dgr - sum_gc_p
    up = sum(ur[b] for b in buses) - sum_gc_u
    sum_cp = sum(cp[b] for b in buses)
    sum_ep_p = sum(ep_p[b] for b in buses)
    sum_ep_u = sum(ep_u[b] for b in buses)
    sum_nuc = sum(nuc[b] for b in buses)
    sum_pen = sum(pen[b] for b in buses)

    trades_mat = {}
    trades_mat["tradeID"] = []
    trades_mat["As"] = []
    trades_mat["Ab"] = []
    trades_mat["Pes"] = []
    trades_mat["Peb"] = []
    trades_mat["costNw"] = []
    trades_mat["penalty"] = []
    trades_mat["cleared"] = []
    for w in trades:
        trades_mat["tradeID"].append(trades[w].tradeID)
        trades_mat["As"].append(trades[w].As)
        trades_mat["Ab"].append(trades[w].Ab)
        trades_mat["Pes"].append(trades[w].Pes)
        trades_mat["Peb"].append(trades[w].Peb)
        trades_mat["costNw"].append(trades[w].costNw)
        trades_mat["penalty"].append(trades[w].penalty)
        if w in trades_selected:
            trades_mat["cleared"].append(1)
        else:
            trades_mat["cleared"].append(0)

    trades_frame = DataFrame(data=trades_mat)

    cashflow_mat = DataFrame([ep_p, ep_u, nuc, pen, ur, cp],
                             index=['ep_p', 'ep_u', 'nuc', 'pen', 'ur', 'cp'])
    cashflow_sum = DataFrame([sum_cp, sum_ep_p, sum_ep_u,
                              sum_nuc, sum_pen, sum_dgr, sum_dgp, up],
                             index=['sum_cp', 'sum_ep_p', 'sum_ep_u',
                                    'sum_nuc', 'sum_pen', 'sum_dgr',
                                    'sum_dgp', 'up'])

    data_mat = DataFrame([trades_dis, nwcharge_dis, echarge_dis],
                         index=['trades_dis', 'nwcharge_dis', 'echarge_dis'])
    data_mat = data_mat.transpose()

    print("partlevel = ", partlevel)
    print("status = ", status)
    print("max loading = ", max(loading[li] for li in lines)*100, "%")
    print("max voltage = ", round(max(voltage[b] for b in buses), 3))
    print("min voltage = ", round(min(voltage[b] for b in buses), 3))
    print("max dlmp = ", round(max(dlmp[b] for b in buses), 3))
    print("min dlmp = ", round(min(dlmp[b] for b in buses), 3))
    method = "pc"

    options = "_"+str(trade_scale)+"_"+str(deltaP)+"_"+str(pen_scale) + \
        "_"+str(pen_start)+"_"+str(clearing)

    basefile = "C:\\Users\\Colton\\github\\py2P\\results\\"
    dirname = basefile+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+str(100*partlevel)+options
    busfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+str(100*partlevel)+options+"_busframe.csv"
    nodefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+str(100*partlevel)+options+"_nodeframe.csv"
    linefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+str(100*partlevel)+options+"_lineframe.csv"
    dlmpfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+str(100*partlevel)+options+"_dlmpframe.csv"
    genfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+str(100*partlevel)+options+"_genframe.csv"
    cashflowmatfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+str(100*partlevel)+options+"_cfmat.csv"
    cashflowsumfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+str(100*partlevel)+options+"_cfsum.csv"
    datafile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+str(100*partlevel)+options+"_data.csv"
    tradefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+str(100*partlevel)+options+"_trades.csv"

    makedirs(dirname, exist_ok=True)
    bus_frame.to_csv(busfile)
    cashflow_mat.to_csv(cashflowmatfile)
    cashflow_sum.to_csv(cashflowsumfile)
    NodeInfo.to_csv(nodefile)
    LineInfo.to_csv(linefile)
    DLMPInfo.to_csv(dlmpfile)
    GenInfo.to_csv(genfile)
    data_mat.to_csv(datafile)
    trades_frame.to_csv(tradefile)

    return dispatch_stack, dlmp_stack, pnm_stack, trades, trades_selected, agents
