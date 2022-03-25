# P2PTradeType.py

class Trade:
    def __init__(self, tradeID, kind, Pes, Peb, As, Ab, costNw, penalty):
        self.tradeID = tradeID
        self.kind = kind
        self.Pes = Pes
        self.Peb = Peb
        self.As = As
        self.Ab = Ab
        self.cleared = 0
        self.blocked = 0
        self.costNw = costNw
        self.penalty = penalty


class Agent:
    def __init__(self, agentID, location, kind, pgMin, pgMax, costFn, Pd, Qd,
                 Tg, Td):
        self.agentID = agentID
        self.location = location
        self.kind = kind
        self.pgMin = pgMin
        self.pgMax = pgMax
        self.costFn = costFn
        self.Pd = Pd
        self.Qd = Qd
        self.Tg = Tg
        self.Td = Td
        # Variable to limit utility buying/selling if causing congestion
        self.utilityfixed = 0
        self.pdutility = 0
        self.pgutility = 0
        self.gpenalty = 0
        self.dpenalty = 0
