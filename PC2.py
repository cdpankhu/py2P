# PC.py
from pandas import DataFrame, Series
from py2P.P2PTradeType import Agent
from py2P.NetworkLoad import networkload
from py2P.TradeFunction import generatetrade, selecttrade
from py2P.DLMP import calculatedlmp
from py2P.Coefficients import calculateptdf, calculatezth
from py2P.makeJac import makeJacVSC
from py2P.penaltyutil import vm_p_sum
from math import pow, sqrt
from time import strftime
from os import makedirs
from numpy import sign


def pc(testsystem, trade_scale, deltaP, pen_scale, pen_start, clearing):
    buses, lines, generators, datamat, ppc = networkload(testsystem)
    windset = {}
    windex = 1

    Zth = calculatezth(buses, lines)

    for g in generators:
        if generators[g].gtype == "Wind":
            windset[windex] = generators[g].gindex
            windex += 1

    # Note: the root method below only works if gen 1 is slack gen at node 1
    # Currently root was referring to gen index as opposed to root node location
    root = -1
    for g in generators:
        if generators[g].gtype == "Root":
            root = generators[g].location

    if root == -1:
        print("No generator at root node")
        generators[1].gtype = "Root"
        root = int(generators[1].location)

    ptdf = calculateptdf(lines, buses, root)

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

    # Defining trading agents at non-root buses
    # Current definition defines a single agent of type 1 at each demand bus.
    # Current definition defines duplicate generator agents for each gen at bus
    # May need to refactor with addition of individual assets at single bus
    # multiple assets at single bus could act as aggregation as well???
    # Current fixed tariff rates for buying and selling to grid
    # in practice, utility could set based on program, time, etc.
    # Tg is the feed in tariff rate, Td is the cost of utility supply
    Tg = 9
    Td = 60
    agents = {}
    agentID = 1
    for b in buses:
        if not B_gn[b]:
            agents[agentID] = Agent(
                b, b, 1, 0, 0, [0, 0, 0], buses[b].Pd, buses[b].Qd, Tg, Td)
            agentID += 1
        else:
            for g in B_gn[b]:
                if generators[g].gtype != "Root":
                    agents[agentID] = Agent(
                        b, b, 1, generators[g].Pmin, generators[g].Pmax,
                        generators[g].cost, buses[b].Pd, buses[b].Qd, Tg, Td)
                    agentID += 1

    # Defining trades and result matrices
    trades = generatetrade(agents, trade_scale)
    trades_old = generatetrade(agents, trade_scale)

    # Set trade nuc based on reference system usage charges
    Ez = 1
    for w in trades:
        trades[w].costNw = Ez*abs(Zth[trades[w].As, trades[w].Ab])

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
    status_old = 0
    dlmp = {}
    dlmp_stack = []
    dispatch_stack = []
    pnm_stack = []
    DLMPInfo = DataFrame()
    NodeInfo = DataFrame()
    LineInfo = DataFrame()
    GenInfo = DataFrame()
    LineState = Series()
    NodeState = Series()
    LineState_old = Series()
    NodeState_old = Series()
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
    changestate = 0
    pg = {}
    pd = {}
    pg_old = {}
    pd_old = {}

    # need to add timer code for performance assessment
    # Outer loop of iterative process
    while True:
        # print("Outer loop")
        outer_iter += 1
        accepted = {}
        rejected = {}
        Revenue_Sell = {}
        Cost_Buy = {}
        Cost_Network = {}
        Cost_Penalty = {}
        Cost_DG = {}
        trade_prices = []
        trade_penalty_deltas = []
        if changestate == 0:
            trades_old = trades
            trades_selected_old = trades_selected
            trades_rest_old = trades_rest
            status_old = status
            LineState_old = LineState
            NodeState_old = NodeState
            for k in pg:
                pg_old[k] = pg[k]
            for k in pd:
                pd_old[k] = pd[k]
        for w in trades:
            trade_prices.append(trades[w].Pes)
            trade_penalty_deltas.append(
                abs(trades[w].penalty - penalty_old[w]))
        min_trade_price = min(trade_prices)
        max_penalty_delta = max(trade_penalty_deltas)
        # set trade prices as the minimum final price minus max delta
        for w in trades:
            # trades[w].Pes = trades[w].Pes - max_penalty_delta - max_nuc_delta
            # trades[w].Peb = trades[w].Peb - max_penalty_delta - max_nuc_delta
            if trades[w].Pes < 0:
                trades[w].Pes = 0
            if trades[w].Peb < 0:
                trades[w].Peb = 0
            penalty_old[w] = trades[w].penalty
            # trades[w].cleared = 0
        max_penalty_delta = 0
        max_nuc_delta = 0

        # Old inner loop, select optimal set of trades
        accepted = {}
        rejected = {}
        pg = {}
        pd = {}

        gap = 1e-3
        for agentID in agents:
            print("agent: " + str(agentID) + " iter: " + str(outer_iter))
            status, temp_accepted, temp_rejected, temp_pg, temp_pgutility, \
                temp_pgpeer, temp_pd, temp_pdutility, temp_pdpeer = \
                selecttrade(trades, agents, agentID, gap, trade_scale)
            accepted[agentID] = temp_accepted
            rejected[agentID] = temp_rejected
            pg[agentID, "pg"] = temp_pg
            pg[agentID, "pgutility"] = temp_pgutility
            pg[agentID, "pgpeer"] = temp_pgpeer
            pd[agentID, "pd"] = temp_pd
            pd[agentID, "pdutility"] = temp_pdutility
            pd[agentID, "pdpeer"] = temp_pdpeer

            if status != 2:
                print("error with select trade")
                return trades, agents, 0, agentID, trade_scale, 0

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
            #       trades[w].penalty, trades[w].Pes-trades[w].costNw,
            #       trades[w].Peb+trades[w].penalty+trades[w].costNw,
            #       w in accepted[trades[w].As], w in accepted[trades[w].Ab])

        trades_selected = []
        trades_rest = []
        # Populate trades selected and order by nuc lowest to highest
        for w in trades:
            if (w in accepted[trades[w].Ab]) and (w in accepted[trades[w].As]):
                if len(trades_selected) == 0:
                    trades_selected.append(w)
                else:
                    inserted = 0
                    for i in range(len(trades_selected)):
                        # insert prior to trade which costNw is less than
                        if trades[w].costNw < trades[trades_selected[i]].costNw:
                            trades_selected.insert(i, w)
                            inserted = 1
                            break
                    # costNw is greater than all others, put at end of list
                    if inserted == 0:
                        trades_selected.append(w)
            else:
                trades_rest.append(w)

        # Currently this uses agents, not sure if it makes sense for buses
        for i in agents:
            dispatch[i] = 0
            gencon[i] = 0
            generation[i] = 0
            consumption[i] = 0

        # Determine dispatch contributions from P2P trades
        for w in trades_selected:
            dispatch[agents[trades[w].Ab].location] -= 1 * trade_scale
            dispatch[agents[trades[w].As].location] += 1 * trade_scale
            gencon[agents[trades[w].Ab].location] -= 1 * trade_scale
            gencon[agents[trades[w].As].location] += 1 * trade_scale
            consumption[agents[trades[w].Ab].location] -= 1 * trade_scale
            generation[agents[trades[w].As].location] += 1 * trade_scale

        # Add contributions from buy/sell to utility
        for a in agents:
            dispatch[a] = round(dispatch[a] - pd[a, "pdutility"]
                                + pg[a, "pgutility"], 3)
            gencon[a] = round(gencon[a] - pd[a, "pdutility"]
                              + pg[a, "pgutility"], 3)
            consumption[a] = round(consumption[a] - pd[a, "pdutility"], 3)
            generation[a] = round(generation[a] + pg[a, "pgutility"], 3)

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
                dispatch_peerG, buses, generators, lines, SMP, gensetP, gensetU
                )

        # Update network charge
        # Trade state is feasible and DLMPs become network charges
        if status == 2 and changestate == 0:
            # For all trades, update network charge
            for w in trades:
                costNwold = trades[w].costNw
                costNwdelta = abs(costNwold - trades[w].costNw)
                # track largest abs nuc change
                if costNwdelta > max_nuc_delta:
                    max_nuc_delta = costNwdelta
        # Trade state infeasible, update penalty on violated trades
        # Trades are successsively removed from the trade stack according to
        # highest network charges first
        elif changestate == 0:
            penalty_iter += 1
            NodeState = NodeInfo.loc['status']
            LineState = LineInfo.loc['status']
            print(NodeState, LineState)
            nodestatus = sum(abs(NodeState[b]) for b in buses)
            linestatus = sum(abs(LineState[li]) for li in lines)
            # Run relaxed power flow to get feasible voltages and DLMPs
            print("Relaxed model")
            status2, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
                calculatedlmp(dispatch_peerG, buses, generators, lines,
                              SMP, gensetP, gensetU, NodeInfo=NodeState,
                              LineInfo=LineState)
            # Update penalty for trades according to proportional
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
                    # Adjust congestion for P2P trades
                    for w in trades:
                        wptdf = round(ptdf[li, agents[trades[w].As].location,
                                           agents[trades[w].Ab].location], 4)

                        trades[w].penalty += LineState[li] * wptdf * pen_scale
                        # trades[w].penalty += (LineState[li] * penpool
                        #                       * wptdf / fwptdfsum)
                    # Adjust congestion for P2utility trades
                    for a in agents:
                        wptdfg = round(ptdf[li, agents[a].location, root], 4)
                        wptdfd = round(ptdf[li, root, agents[a].location], 4)

                        agents[a].gpenalty += LineState[li] * \
                            wptdfg * pen_scale
                        agents[a].dpenalty += LineState[li] * \
                            wptdfd * pen_scale

                    # for w in trades:
                    #     print(li, LineState[li], agents[trades[w].As].location,
                    #           agents[trades[w].Ab].location, trades[w].Pes,
                    #           trades[w].Peb, trades[w].costNw,
                    #           trades[w].penalty,
                    #           trades[w].costNw+trades[w].penalty,
                    #           trades[w].cleared)
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
        if changestate == 0:
            dispatch_stack.append(dispatch_temp)
            dlmp_stack.append(dlmp_temp)
            pnm_stack.append(pnm_temp)

        # !! Congestion clearing process !!
        # Stable solution shifted from feasible to infeasible after penalty
        # Revert penalty back once and lock in trades for congested line
        if status == 2 and changestate == 0 and status_old != 2:
            print(LineState, LineState_old)
            for li in lines:
                # Check if line congestion was alleviated
                if LineState_old[li] != 0:
                    line_limit = lines[li].u * LineState_old[li]
                    line_loading = 0
                    line_violation = 0
                    # First accept all P2P trades in the non-violating direction
                    for w in trades_selected_old:
                        wptdf = round(ptdf[li, agents[trades[w].As].location,
                                           agents[trades[w].Ab].location], 4)

                        if wptdf != 0 and sign(wptdf) != sign(LineState_old[li]):
                            trades[w].penalty -= LineState_old[li] * wptdf \
                                                 * pen_scale
                            line_loading += wptdf*trade_scale
                            # Lock congestion on that line
                            trades[w].cleared = 1
                    # Fix all non-violating direction P2U trades as is. This
                    # includes trades of 0 with grid
                    for a in agents:
                        wptdfg = round(ptdf[li, agents[a].location, root])
                        wptdfd = round(ptdf[li, root, agents[a].location])

                        if wptdfg != 0 and sign(wptdfg) != sign(LineState_old[li]):
                            agents[a].gpenalty -= LineState_old[li] * wptdfg \
                                                 * pen_scale
                            agents[a].utilityfixed = 1
                            agents[a].pgutility = pg_old[a, "pgutility"]
                            line_loading += wptdfg*agents[a].pgutility
                        if wptdfd != 0 and sign(wptdfd) != sign(LineState_old[li]):
                            agents[a].dpenalty -= LineState_old[li] * wptdfd \
                                                 * pen_scale
                            agents[a].utilityfixed = 1
                            agents[a].pdutility = pd_old[a, "pdutility"]
                            line_loading += wptdfd*agents[a].pdutility
                    # Accept violating P2P trades while line capacity remains
                    for w in trades_selected_old:
                        wptdf = round(ptdf[li, agents[trades[w].As].location,
                                           agents[trades[w].Ab].location], 4)

                        if wptdf != 0 and sign(wptdf) == sign(LineState_old[li]):
                            line_loading += wptdf*trade_scale
                            if line_limit > 0 and line_loading >= line_limit:
                                line_violation = 1
                            elif line_limit < 0 and line_loading <= line_limit:
                                line_violation = 1
                            if line_violation == 0:
                                trades[w].penalty -= LineState_old[li] * wptdf \
                                    * pen_scale
                                # Lock that trade as being in the optimal set
                                trades[w].cleared = 1
                            else:
                                trades[w].blocked = 1
                    # Block all non-selected P2P trades impacting congestion
                    for w in trades_rest_old:
                        wptdf = round(ptdf[li, agents[trades[w].As].location,
                                           agents[trades[w].Ab].location], 4)

                        if wptdf != 0:
                            trades[w].blocked = 1
                    # If any capacity remaining, allocate to P2U transactions
                    # This also fixes trades if 0 in the violating direction
                    # First, reduce penalties for violating tradees and order
                    # agents according to lowest NUC with root
                    nuc_agents = []
                    for a in agents:
                        wptdfg = round(ptdf[li, agents[a].location, root])
                        wptdfd = round(ptdf[li, root, agents[a].location])
                        if wptdfg != 0 or wptdfd != 0:
                            agents[a].utilityfixed = 1
                            agents[a].dpenalty -= LineState_old[li] * wptdfd \
                                * pen_scale
                            agents[a].gpenalty -= LineState_old[li] * wptdfg \
                                * pen_scale

                        if len(nuc_agents) == 0:
                            nuc_agents.append(a)
                        else:
                            inserted = 0
                            for i in range(len(nuc_agents)):
                                # insert prior to trade which Zth less than
                                if Zth[a, root] < Zth[nuc_agents[i], root]:
                                    nuc_agents.insert(i, a)
                                    inserted = 1
                                    break
                            # costNw is greater than all others, put at end of list
                            if inserted == 0:
                                nuc_agents.append(a)
                    # Next, iterate through agents in order of nuc for clearing
                    for a in nuc_agents:
                        wptdfg = round(ptdf[li, agents[a].location, root])
                        wptdfd = round(ptdf[li, root, agents[a].location])

                        if wptdfg != 0 and sign(wptdfg) == sign(LineState_old[li]):
                            if pg_old[a, "pgutility"] >= abs(line_limit - line_loading) - trade_scale:
                                agents[a].pgutility = abs(
                                    line_limit - line_loading) - trade_scale
                                line_loading += wptdfg * agents[a].pgutility
                                break
                            else:
                                agents[a].pgutility = pg_old[a, "pgutility"]
                                line_loading += wptdfg * agents[a].pgutility
                        if wptdfd != 0 and sign(wptdfd) == sign(LineState_old[li]):
                            if pd_old[a, "pdutility"] >= abs(line_limit - line_loading) - trade_scale:
                                agents[a].pdutility = abs(
                                    line_limit - line_loading) - trade_scale
                                line_loading += wptdfd * agents[a].pdutility
                                break
                            else:
                                agents[a].pdutility = pd_old[a, "pdutility"]
                                line_loading += wptdfd * agents[a].pdutility

                                # If gencon is unchanged, break while loop, feasible stable solution
        if gencon_old == gencon and status == 2 and changestate == 0:
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
    seller_rev = {}
    seller_rev_peer = {}
    seller_rev_utility = {}
    buyer_cost = {}
    buyer_cost_peer = {}
    buyer_cost_utility = {}
    buyer_cost_system = {}
    buyer_cost_nuc = {}
    buyer_cost_pen = {}
    seller_cost_nuc = {}
    seller_cost_gen = {}
    seller_cost = {}
    utility_cost = 0
    utility_revenue = 0

    for b in buses:
        demand_P2P[b] = 0
        ep_u[b] = 0
        ep_p[b] = 0
        nuc[b] = 0
        pen[b] = 0
        seller_rev[b] = 0
        seller_rev_peer[b] = 0
        seller_rev_utility[b] = 0
        buyer_cost[b] = 0
        buyer_cost_peer[b] = 0
        buyer_cost_utility[b] = 0
        buyer_cost_system[b] = 0
        buyer_cost_nuc[b] = 0
        buyer_cost_pen[b] = 0
        seller_cost_nuc[b] = 0
        seller_cost_gen[b] = 0
        seller_cost[b] = 0

    # Cost of the system supplied energy due to trade scale rounding
    for a in agents:
        b = agents[a].location
        buyer_cost_system[b] += (agents[a].Pd - pd[a, "pdpeer"]
                                 - pd[a, "pdutility"]) * agents[a].Td
        seller_rev_utility[b] += pg[a, "pgutility"] * agents[a].Tg
        buyer_cost_utility[b] += pd[a, "pdutility"] * agents[a].Tg

    for g in gensetP:
        b = generators[g].location
        seller_cost_gen[b] += generators[g].cost[1] * dispatch_peerG[g]

    utility_cost = generators[root].cost[1] * GenInfo[root]["pg"]

    for m in B_d:
        demand_P2P[m] = sum(trades_dis[n, m] for n in B_g)
        ep_u[m] = (buses[m].Pd - demand_P2P[m])*dlmp[m]*trade_scale

    for w in trades_selected:
        ep_p[agents[trades[w].Ab].location] += trades[w].Pes*trade_scale
        nuc[agents[trades[w].Ab].location] += trades[w].costNw*trade_scale
        nuc[agents[trades[w].As].location] += trades[w].costNw*trade_scale
        pen[agents[trades[w].Ab].location] += trades[w].penalty*trade_scale
        buyer_cost_peer[agents[trades[w].Ab].location] += trades[w].Pes*trade_scale
        seller_rev_peer[agents[trades[w].As].location] += trades[w].Pes*trade_scale
        buyer_cost_nuc[agents[trades[w].Ab].location] += trades[w].costNw*trade_scale
        seller_cost_nuc[agents[trades[w].As].location] += trades[w].costNw*trade_scale
        buyer_cost_pen[agents[trades[w].Ab].location] += trades[w].penalty*trade_scale

    for b in buses:
        seller_rev[b] += seller_rev_peer[b] + seller_rev_utility[b]
        buyer_cost[b] += buyer_cost_peer[b] + buyer_cost_utility[b]\
            + buyer_cost_nuc[b] + buyer_cost_pen[b] + buyer_cost_system[b]
        seller_cost[b] += seller_cost_gen[b] + seller_cost_nuc[b]
        utility_revenue += buyer_cost_utility[b] + buyer_cost_system[b]

    ur = {}
    cp = {}
    for b in buses:
        buyer_cost_utility[b] = pd[b, "pdutility"]*Td
        ur[b] = ep_u[b] + nuc[b]
        cp[b] = ur[b] + ep_p[b] + pen[b]

    gc = {}
    for g in generators:
        gc[g] = generators[g].cost[1] \
            * round(pg[generators[g].location, "pg"], 5)

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

    sum_seller_rev = sum(seller_rev[b] for b in buses)
    sum_seller_rev_peer = sum(seller_rev_peer[b] for b in buses)
    sum_seller_rev_utility = sum(seller_rev_utility[b] for b in buses)
    sum_buyer_cost = sum(buyer_cost[b] for b in buses)
    sum_buyer_cost_peer = sum(buyer_cost_peer[b] for b in buses)
    sum_buyer_cost_utility = sum(buyer_cost_utility[b] for b in buses)
    sum_buyer_cost_system = sum(buyer_cost_system[b] for b in buses)
    sum_buyer_cost_nuc = sum(buyer_cost_nuc[b] for b in buses)
    sum_buyer_cost_pen = sum(buyer_cost_pen[b] for b in buses)
    sum_seller_cost_nuc = sum(seller_cost_nuc[b] for b in buses)
    sum_seller_cost_gen = sum(seller_cost_gen[b] for b in buses)
    sum_seller_cost = sum(seller_cost[b] for b in buses)

    trades_mat = {}
    trades_mat["tradeID"] = []
    trades_mat["As"] = []
    trades_mat["Ab"] = []
    trades_mat["Pes"] = []
    trades_mat["Peb"] = []
    trades_mat["costNw"] = []
    trades_mat["penalty"] = []
    trades_mat["cleared"] = []
    trades_mat["blocked"] = []
    trades_mat["passed"] = []
    trades_mat['buyer_pay'] = []
    trades_mat['seller_rev'] = []
    for w in trades:
        trades_mat["tradeID"].append(trades[w].tradeID)
        trades_mat["As"].append(trades[w].As)
        trades_mat["Ab"].append(trades[w].Ab)
        trades_mat["Pes"].append(trades[w].Pes)
        trades_mat["Peb"].append(trades[w].Peb)
        trades_mat["costNw"].append(trades[w].costNw)
        trades_mat["penalty"].append(trades[w].penalty)
        trades_mat["cleared"].append(trades[w].cleared)
        trades_mat["blocked"].append(trades[w].blocked)
        if w in trades_selected:
            trades_mat["passed"].append(1)
        else:
            trades_mat["passed"].append(0)
        trades_mat['buyer_pay'].append(trades[w].Peb + trades[w].costNw
                                       + trades[w].penalty)
        trades_mat['seller_rev'].append(trades[w].Pes - trades[w].costNw)

    trades_frame = DataFrame(data=trades_mat)

    revcost_mat = DataFrame([seller_rev, seller_rev_peer, seller_rev_utility,
                             buyer_cost, buyer_cost_peer, buyer_cost_utility,
                             buyer_cost_system, buyer_cost_nuc,
                             buyer_cost_pen, seller_cost, seller_cost_nuc,
                             seller_cost_gen],
                            index=["seller_rev", "seller_rev_peer",
                                   "seller_rev_utility", "buyer_cost",
                                   "buyer_cost_peer", "buyer_cost_utility",
                                   "buyer_cost_system", "buyer_cost_nuc",
                                   "buyer_cost_pen", "seller_cost",
                                   "seller_cost_nuc", "seller_cost_gen"])
    revcost_sum = DataFrame([sum_seller_rev, sum_seller_rev_peer,
                             sum_seller_rev_utility, sum_buyer_cost,
                             sum_buyer_cost_peer, sum_buyer_cost_utility,
                             sum_buyer_cost_system, sum_buyer_cost_nuc,
                             sum_buyer_cost_pen, sum_seller_cost,
                             sum_seller_cost_nuc, sum_seller_cost_gen,
                             utility_revenue, utility_cost],
                            index=["sum_seller_rev", "sum_seller_rev_peer",
                                   "sum_seller_rev_utility", "sum_buyer_cost",
                                   "sum_buyer_cost_peer",
                                   "sum_buyer_cost_utility",
                                   "sum_buyer_cost_system",
                                   "sum_buyer_cost_nuc",
                                   "sum_buyer_cost_pen", "sum_seller_cost",
                                   "sum_seller_cost_nuc", "sum_seller_cost_gen",
                                   "utility_revenue", "utility_cost"])

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

    print("status = ", status)
    print("max loading = ", max(loading[li] for li in lines)*100, "%")
    print("max voltage = ", round(max(voltage[b] for b in buses), 3))
    print("min voltage = ", round(min(voltage[b] for b in buses), 3))
    print("max dlmp = ", round(max(dlmp[b] for b in buses), 3))
    print("min dlmp = ", round(min(dlmp[b] for b in buses), 3))
    method = "pc2"

    options = "_"+str(trade_scale)+"_"+str(deltaP)+"_"+str(pen_scale) + \
        "_"+str(pen_start)+"_"+str(clearing)

    basefile = "C:\\Users\\Colton\\github\\py2P\\results\\"
    dirname = basefile+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem + "_"+options
    busfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+options+"_busframe.csv"
    nodefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+options+"_nodeframe.csv"
    linefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+options+"_lineframe.csv"
    dlmpfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+options+"_dlmpframe.csv"
    genfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+options+"_genframe.csv"
    cashflowmatfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+options+"_cfmat.csv"
    cashflowsumfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+options+"_cfsum.csv"
    revcostmatfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+options+"_revcostmat.csv"
    revcostsumfile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+options+"_revcostsum.csv"
    datafile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_" + \
        testsystem+"_"+options+"_data.csv"
    tradefile = dirname+"\\"+strftime("%Y-%m-%d-%H-%M-%S_")+method+"_"+testsystem + \
        "_"+options+"_trades.csv"

    makedirs(dirname, exist_ok=True)
    bus_frame.to_csv(busfile)
    cashflow_mat.to_csv(cashflowmatfile)
    cashflow_sum.to_csv(cashflowsumfile)
    revcost_mat.to_csv(revcostmatfile)
    revcost_sum.to_csv(revcostsumfile)
    NodeInfo.to_csv(nodefile)
    LineInfo.to_csv(linefile)
    DLMPInfo.to_csv(dlmpfile)
    GenInfo.to_csv(genfile)
    data_mat.to_csv(datafile)
    trades_frame.to_csv(tradefile)
    print(dirname)

    return dispatch_stack, dlmp_stack, pnm_stack, trades, trades_selected, agents
