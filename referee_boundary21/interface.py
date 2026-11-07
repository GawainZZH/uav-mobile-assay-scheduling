from service_separator25_joint.interface import certify as original

def safe_fixed_cut(a,b,c,diameter):
    return a<=c<=b and max(c-a,b-c)>=diameter

def nearest_safe_cut(a,b,c,diameter):
    if b-a<diameter:return None
    intervals=[(a,b-diameter),(a+diameter,b)]
    candidates=[min(max(c,lo),hi)for lo,hi in intervals if lo<=hi]
    return min(candidates,key=lambda t:(abs(t-c),t))

def certify(inst,boundaries):
    result=original(inst,boundaries);diameter=result['maximum_shortest_road']
    for cut in result['cuts']:
        cut['previous_two_sided_valid']=cut['valid']
        cut['valid']=safe_fixed_cut(cut['left_clear'],cut['right_first'],cut['epoch'],diameter)
        cut['normalization_side']='before'if cut['valid']and cut['left_margin']>=diameter else'after'if cut['valid']else None
    result['valid']=all(c['valid']for c in result['cuts']);result['rule']='one-sided sufficient parked cut'
    return result
