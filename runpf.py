# DLMP.py
from py2P.NetworkLoad import networkload
from py2P.DLMP import calculatedlmp
from datetime import datetime
from pandas import DataFrame, Series


def runpf(testsystem, set_conic, relaxed):

    buses, lines, generators, sBase, datamat, ppc = networkload(testsystem)

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

    # Consider reducing cyclomatic complexity by moving shadow price extraction
    B_gn = {}
    for i in buses:
        B_gn[i] = []
    for g in generators:
        B_gn[generators[g].location].append(g)

    dispatch_peerG = {}
    for g in generators:
        dispatch_peerG[g] = generators[g].Pg

    LineState = Series(data=dict)
    if relaxed == 1:
        for li in lines:
            LineState[li] = 2
    else:
        for li in lines:
            LineState[li] = 0

    SMP = generators[root].cost[1]
    start_time = datetime.now()
    if relaxed == 1:
        status, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
            calculatedlmp(
                dispatch_peerG, buses, generators, lines, sBase, SMP, gensetP, gensetU, conic=set_conic, LineInfo=LineState
                )
    else:
        status, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
            calculatedlmp(
                dispatch_peerG, buses, generators, lines, sBase, SMP, gensetP, gensetU, conic=set_conic
                )
    # status, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo = \
    #     calculatedlmp(
    #         dispatch_peerG, buses, generators, lines, sBase, SMP, gensetP, gensetU, conic=0
    #         )
    print("total time: ", datetime.now() - start_time)
    return status, dlmp, pgextra, NodeInfo, LineInfo, DLMPInfo, GenInfo
