# TradeFunction.py
from .P2PTradeType import Trade
from math import floor, ceil
from gurobipy import Model, GRB


def generatetrade(agents, trade_scale):
    trades = {}
    tradeID = 1

    for s in agents:
        for b in agents:
            up = min(
                floor(agents[s].pgMax/trade_scale),
                ceil(agents[b].Pd/trade_scale))
            # +1 since range is < max
            for k in range(1, up+1):
                w = Trade(tradeID, 1, 0.0, 0.0, s, b, 0)
                trades[tradeID] = w
                tradeID += 1
    return trades


def selecttrade(W, A, agentID, gap, trade_scale):
    # Trade Set W
    m = Model("selecttrade")
    m.Params.OutputFlag = 0
    # Create variables
    # Note, need to add constraint/variable for utility of marginal demand
    pg = m.addVar(name="pg", lb=A[agentID].pgMin, ub=A[agentID].pgMax)
    revenuesell = m.addVar(name="revenuesell")
    costbuy = m.addVar(name="costbuy")
    costnetwork = m.addVar(name="costnetwork")
    costdg = m.addVar(name="costdg")
    u = m.addVars((w for w in W), vtype=GRB.BINARY, name="u")

    m.update()

    # Set objective functions
    obj = revenuesell
    obj -= costbuy
    obj -= costnetwork
    obj -= costdg
    m.setObjective(obj, GRB.MAXIMIZE)

    m.update()

    # Create constraints
    # Define revenue from sales to peers
    m.addConstr(
        revenuesell
        == sum(W[i].Pes * u[i] for i in W if W[i].As == agentID),
        "c0")
    # Define cost of purchasing from peers
    m.addConstr(
        costbuy == sum(W[i].Peb * u[i] for i in W if W[i].Ab == agentID),
        "c1")
    # Define cost of network on transactions
    m.addConstr(
        costnetwork == trade_scale * (
            sum(W[i].costNw * u[i] for i in W if W[i].Ab == agentID)
            + sum(W[i].costNw * u[i] for i in W if W[i].As == agentID)),
        "c2")
    # Define cost of DG generation operation
    m.addConstr(
        costdg == (
            A[agentID].costFn[0]*pg*pg + A[agentID].costFn[1]*pg
            + A[agentID].costFn[2]
            ),
        "c3")
    # Require that power from Dg is equal to the sum of selling trades
    # Need to make sure trading to oneself is okay!
    m.addConstr(
        pg/trade_scale == sum(u[i] for i in W if W[i].As == agentID),
        "c4")
    # Requires that the peer demand is equal to the sum of purchased trades
    # Currently Pd is fixed, no marginal value of additional demand. FIX!!!
    m.addConstr(
        floor(A[agentID].Pd/trade_scale)
        == sum(u[i] for i in W if W[i].Ab == agentID),
        "c5")

    # Run optimization
    for i in W:
        if not (W[i].As == agentID or W[i].Ab == agentID):
            m.addConstr(u[i] == 0)

    m.update()

    m.optimize()
    v = m.getVars()
    w = {}
    tradeID = 1
    for i in range(5, len(v)):
        w[tradeID] = v[i].x
        tradeID += 1
    accepted = [i for i in w if w[i] == 1]
    rejected = [i for i in w if w[i] == 0]
    pg = m.getVarByName("pg").x
    revenuesell = m.getVarByName("revenuesell").x
    costbuy = m.getVarByName("costbuy").x
    costnetwork = m.getVarByName("costnetwork").x
    costdg = m.getVarByName("costdg").x
    return accepted, rejected, pg, revenuesell, costbuy, costnetwork, costdg
