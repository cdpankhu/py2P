# P2PTradeType.py

class Trade:
    def __init__(self, tradeID, kind, Pes, Peb, As, Ab, costNw, penalty):
        self.tradeID = tradeID
        self.kind = kind
        self.Pes = Pes
        self.Peb = Peb
        self.As = As
        self.Ab = Ab
        self.costNw = costNw
        self.penalty = penalty


class Agent:
    def __init__(self, agentID, location, kind, pgMin, pgMax, costFn, Pd, Qd):
        self.agentID = agentID
        self.location = location
        self.kind = kind
        self.pgMin = pgMin
        self.pgMax = pgMax
        self.costFn = costFn
        self.Pd = Pd
        self.Qd = Qd
