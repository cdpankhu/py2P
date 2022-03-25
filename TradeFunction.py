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
                w = Trade(tradeID, 1, 0.0, 0.0, s, b, 0.0, 0.0)
                trades[tradeID] = w
                tradeID += 1
    return trades


def selecttrade(W, A, agentID, gap, trade_scale):
    # Trade Set W
    m = Model("selecttrade")
    m.Params.OutputFlag = 0
    m.Params.DualReductions = 0
    # Create variables
    # Note, need to add constraint/variable for utility of marginal demand
    # Total energy sold by peer
    pg = m.addVar(name="pg", lb=A[agentID].pgMin, ub=A[agentID].pgMax)
    # Total energy sold by peer to other peers
    pgpeer = m.addVar(name="pgpeer", lb=0)
    # Total energy sold by peer to utility/grid
    pgutility = m.addVar(name="pgutility", lb=0)
    pd = m.addVar(name="pd", lb=0)
    pdpeer = m.addVar(name="pdpeer", lb=0)
    pdutility = m.addVar(name="pdutility", lb=0)
    revenuepeer = m.addVar(name="revenuepeer")
    revenueutility = m.addVar(name="revenueutility")
    revenuesell = m.addVar(name="revenuesell")
    costpeer = m.addVar(name="costpeer")
    costutility = m.addVar(name="costutility")
    costbuy = m.addVar(name="costbuy")
    costnetwork = m.addVar(name="costnetwork", lb=float('-inf'))
    costpenalty = m.addVar(name="costpenalty", lb=float('-inf'))
    costdg = m.addVar(name="costdg")
    u = m.addVars((w for w in W), vtype=GRB.BINARY, name="u")

    m.update()

    # Set objective functions
    obj = revenuesell
    obj -= costbuy
    obj -= costnetwork
    obj -= costpenalty
    obj -= costdg
    m.setObjective(obj, GRB.MAXIMIZE)

    m.update()

    # Create constraints
    # Define revenue from sales to peers
    m.addConstr(
        revenuepeer == trade_scale * (
            sum(W[i].Pes * u[i] for i in W if W[i].As == agentID)),
        "c0")
    m.addConstr(
        revenueutility == pgutility*(A[agentID].Tg - A[agentID].gpenalty),
        "c1")
    m.addConstr(
        revenuesell == revenuepeer + revenueutility,
        "c2")
    # Define cost of purchasing from peers
    m.addConstr(
        costpeer == trade_scale * (
            sum(W[i].Peb * u[i] for i in W if W[i].Ab == agentID)),
        "c3")
    m.addConstr(
        costutility == pdutility*(A[agentID].Td + A[agentID].dpenalty),
        "c4")
    m.addConstr(
        costbuy == costpeer + costutility,
        "c5")
    # Define cost of network on transactions, buyer and seller share cost
    m.addConstr(
        costnetwork == trade_scale * (
            sum(W[i].costNw * u[i] for i in W if W[i].Ab == agentID)
            + sum(W[i].costNw * u[i] for i in W if W[i].As == agentID)),
        "c6")
    # Define cost of penalty on transactions, buyer pays full penalty
    m.addConstr(
        costpenalty == trade_scale * (
            sum(W[i].penalty * u[i] for i in W if W[i].Ab == agentID)),
        "c7")
    # Define cost of DG generation operation
    m.addConstr(
        costdg == (
            A[agentID].costFn[0]*pg*pg + A[agentID].costFn[1]*pg
            + A[agentID].costFn[2]),
        "c8")
    # Define pgpeer as the sum of all energy sold to other peers
    m.addConstr(
        pgpeer == trade_scale * sum(u[i] for i in W if W[i].As == agentID),
        "c9")
    # Require that pg is equal to the sum of total energy sold to peers plus
    # total energy sold to utility
    m.addConstr(
        pg == pgpeer + pgutility,
        "c10")
    # Require that total energy purchased is equal to demand. Currently demand
    # is fixed, however value of marginal demand should be considered.
    m.addConstr(
        pd == trade_scale * floor(A[agentID].Pd/trade_scale),
        "c11")
    # Define pdpeer as the sum of all energy bought from other peers
    m.addConstr(
        pdpeer == trade_scale * sum(u[i] for i in W if W[i].Ab == agentID),
        "c12")
    # Define total energy purchased as the energy purchased from peers plus
    # the energy purchased from the utility
    m.addConstr(
        pd == pdpeer + pdutility,
        "c13")

    # Add additional trade specific constraints
    for i in W:
        # If a trade does not include the agent, set select to 0
        if not (W[i].As == agentID or W[i].Ab == agentID):
            m.addConstr(u[i] == 0, 'c14['+str(i)+"]")
        else:
            # If a trade is set to cleared, select the trade in optimal set
            if W[i].cleared == 1:
                m.addConstr(u[i] == 1, 'c15['+str(i)+"]")
            # If a trade is blocked, exclude the trade from optimal set
            if W[i].blocked == 1:
                m.addConstr(u[i] == 0, 'c16['+str(i)+"]")
    # Constraints on trading to
    if A[agentID].utilityfixed == 1:
        # If setpoint of agent is fixed, ensure demand meets that expected
        m.addConstr(
            pdutility == A[agentID].pdutility,
            'c17')
        # If setpoint of agent is fixed, ensure supply meets that expected
        m.addConstr(
            pgutility == A[agentID].pgutility,
            'c18')

    m.update()

    m.optimize()

    status = m.Status
    if not (status == GRB.OPTIMAL):
        return status, [], [], 0, 0, 0, 0, 0, 0, m
    v = m.getVars()
    w = {}
    tradeID = 1
    # 6 to skip the first 6 vars
    for i in range(15, 15 + len(W)):
        w[tradeID] = v[i].x
        tradeID += 1
    accepted = [i for i in w if w[i] == 1]
    rejected = [i for i in w if w[i] == 0]
    pg = m.getVarByName("pg").x
    pd = m.getVarByName("pd").x
    revenuesell = m.getVarByName("revenuesell").x
    costbuy = m.getVarByName("costbuy").x
    costnetwork = m.getVarByName("costnetwork").x
    costpenalty = m.getVarByName("costpenalty").x
    costdg = m.getVarByName("costdg").x
    pgutility = m.getVarByName("pgutility").x
    pgpeer = m.getVarByName("pgpeer").x
    pdutility = m.getVarByName("pdutility").x
    pdpeer = m.getVarByName("pdpeer").x

    return status, accepted, rejected, pg, pgutility, pgpeer, pd, pdutility, \
        pdpeer
