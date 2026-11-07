"""Min-plus task-set partitions; table unknowns are not infeasibilities."""
from fractions import Fraction as F


def partition(tables,counts,njobs):
    counts=tuple(counts)
    if len(tables)!=len(counts)or sum(counts)!=njobs or any(type(v)is not int or v<0 for v in counts):
        raise ValueError('consistent complete responsibility counts required')
    full=(1<<njobs)-1;states={0:(F(0),())};transitions=0
    for table,count in zip(tables,counts):
        new={}
        for mask,(cost,path)in states.items():
            available=full^mask;chosen=available
            while True:
                if chosen.bit_count()==count and table.get(chosen)is not None:
                    value=cost+F(table[chosen]);target=mask|chosen;transitions+=1
                    if target not in new or value<new[target][0]:new[target]=(value,path+(chosen,))
                if chosen==0:break
                chosen=(chosen-1)&available
        states=new
    if full not in states:return dict(status='infeasible',objective=None,subsets=None,transitions=transitions)
    value,path=states[full]
    return dict(status='optimal',objective=value,subsets=path,transitions=transitions)


def combine_intervals(tables,counts,njobs):
    lower=[];upper=[]
    for table,count in zip(tables,counts):
        low={};up={}
        for mask in range(1<<njobs):
            if mask.bit_count()!=count:continue
            r=table.get(mask,{})
            if r.get('status')=='infeasible':low[mask]=None;up[mask]=None;continue
            low[mask]=F(r.get('lower_bound',0));value=r.get('objective')
            if value is not None:
                if r.get('plan')is None:raise ValueError('unverified upper has no plan')
                if low[mask]>F(value):raise ValueError('subset lower exceeds its physical upper')
                up[mask]=F(value)
            else:up[mask]=None
        lower.append(low);upper.append(up)
    lo=partition(lower,counts,njobs);hi=partition(upper,counts,njobs)
    if lo['objective']is None and hi['objective']is not None:raise ValueError('false lower infeasibility')
    if lo['objective']is not None and hi['objective']is not None and lo['objective']>hi['objective']:
        raise ValueError('partition interval contradiction')
    closed=lo['objective']is not None and lo['objective']==hi['objective']
    return dict(status='infeasible'if lo['objective']is None else'optimal'if closed else'feasible'if hi['objective']is not None else'unresolved',
        optimal=closed,lower_bound=lo['objective'],objective=hi['objective'],
        lower_partition=lo['subsets'],upper_partition=hi['subsets'],
        lower_transitions=lo['transitions'],upper_transitions=hi['transitions'])
