#penaltyutil.py


def vm_p_sum(b, trades, agents, vm_p_sc):
    vsc_pos_sum = sum(abs(round(
        vm_p_sc[b, agents[trades[w].As].location]
        - vm_p_sc[b, agents[trades[w].Ab].location], 4))
        for w in trades if round(vm_p_sc[b, agents[trades[w].As].location]
                                 - vm_p_sc[b, agents[trades[w].Ab].location], 4)
        > 0)

    vsc_neg_sum = sum(abs(round(
        vm_p_sc[b, agents[trades[w].As].location]
        - vm_p_sc[b, agents[trades[w].Ab].location], 4))
        for w in trades if round(vm_p_sc[b, agents[trades[w].As].location]
                                 - vm_p_sc[b, agents[trades[w].Ab].location], 4)
        < 0)

    return vsc_pos_sum, vsc_neg_sum
