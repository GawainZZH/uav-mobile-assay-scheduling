"""Atomic lower-band updates; no oracle call or new feasibility assertion."""
from dataclasses import replace
from fractions import Fraction as F
from service_separator25_joint.solve import signature
from block_flow_solver.bands import exact
from block_flow_solver.matrix import row


def apply_prices(m,prices):
    d=m.domain
    if d.inst is None or d.stages is None:raise ValueError('audited physical stage domain required')
    identities=[signature(replace(local,groups=1))for local,ids in d.stages]
    compatible={k:[]for k in range(len(d.stages))}
    for p in prices:
        for k,key in enumerate(identities):
            if key==p.instance_signature and p.crew_capacity==d.inst.drones:compatible[k].append(p)
    changes=[]
    for i,b in enumerate(m.bands):
        if b['infeasible']:continue
        k=b['stage'];ids=d.stages[k][1];S=tuple(ids.index(j)for j in b['subset'])
        if 'local_subset'in b and S!=tuple(b['local_subset']):raise ValueError('global/local responsibility mapping differs')
        best=F(str(b['lower']));witness=None
        for p in compatible[k]:
            bound=p.lower(S,b['start'],b['end'])
            if bound>best:best=bound;witness=p
        if best>F(str(b['upper'])):raise ValueError(('shared price contradicts a retained upper',i,best,b['upper']))
        if not b['subset']:
            if best>F(str(b['lower'])):raise ValueError('price contradicts an exact empty road family')
            continue
        if witness is not None:
            if i not in m.z:raise ValueError('nonempty family cost coordinate required')
            changes.append(dict(family=i,stage=k,subset=tuple(b['subset']),local_subset=S,start=b['start'],end=b['end'],
                old_lower=F(str(b['lower'])),lower=best,certificate_id=witness.certificate_id))
    # All contradictions are checked before changing any row or bound.
    for change in changes:
        i=change['family'];b=m.bands[i];b['lower']=change['lower'];b['exact']=exact(b)
        if b['exact']:b['proof_scope']='transferred_rational_LP_price_and_physical_plan'
        b['lower_price_certificate']=change['certificate_id']
        e={j:change['lower']for j in m.family_columns[i]};e[m.z[i]]=F(-1)
        m.rows.append(row(e,0,False,'shared_G1_lower'))
    return changes
