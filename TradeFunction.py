# TradeFunction.py
from .P2PTradeType import Trade, Agent
from math import floor, ceil
from gurobipy import *


def generatetrade(agents, trade_scale):
    trades = []
    tradeID = 1

    agentset = list(range(len(agents)))
    for s in agentset:
        for b in agentset:
            up = min(
                floor(agents[s].pgMax/trade_scale),
                ceil(agents[b].Pd/trade_scale))
            for k in range(1, up):
                w = Trade(tradeID, 1, 0.0, 0.0, s, b, 0)
                trades.append(w)
                tradeID += 1
    return trades


def selecttrade(W, A, agentID, gap, trade_scale):
    # Trade Set W
    m = Model("selecttrade")

    # Creat variables
    u = {}
    for w in range(len(W)):
        u[w] = m.addVar(vtype=GRB.BINARY, name="u")

    # Set objective functions

    # Create constraints

    # Run optimization
