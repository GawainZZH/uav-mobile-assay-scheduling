from referee_boundary21.interface import safe_fixed_cut,nearest_safe_cut
from service_separator25_joint.interface import certify as original

def safe_joint_cuts(a,b,cuts,diameter):
    points=[a]+sorted(cuts)+[b]
    return a<=b and all(a<=c<=b for c in cuts)and max((v-u for u,v in zip(points,points[1:])),default=b-a)>=diameter

def certify(inst,boundaries):
    result=original(inst,boundaries);diameter=result['maximum_shortest_road'];clusters={}
    for cut in result['cuts']:
        cut['previous_two_sided_valid']=cut['valid']
        cut['valid']=safe_fixed_cut(cut['left_clear'],cut['right_first'],cut['epoch'],diameter)
        if cut['valid']:clusters.setdefault((cut['left_clear'],cut['right_first']),[]).append(cut)
    for (a,b),cluster in clusters.items():
        overlapping=any(max(a,x)<min(b,y)for x,y in clusters if(x,y)!=(a,b))
        jointly_safe=not overlapping and safe_joint_cuts(a,b,[c['epoch']for c in cluster],diameter)
        for cut in cluster:
            cut['valid']=jointly_safe
            cut['joint_normalization']=jointly_safe
            cut['normalization_side']='before'if jointly_safe and cut['left_margin']>=diameter else'after'if jointly_safe else None
    result['valid']=all(c['valid']for c in result['cuts']);result['rule']='joint one-sided sufficient parked cuts'
    return result
