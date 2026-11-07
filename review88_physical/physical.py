"""Independent interval replay, without graph transitions or optimization."""
from collections import Counter
from fractions import Fraction
from .types import Route


def validate(inst, plan, require_all=True):
    """Check one group's concrete route/truck/assay calendar and exact cost."""
    errors = []
    road_cost = Fraction(0)
    flight_cost = Fraction(0)
    age_cost = Fraction(0)
    stationary = [None]*inst.H
    drone_use = [0]*inst.H
    station_use = [0]*inst.H
    nodes = {(0, 0)}
    last = (0, 0)
    routes = []
    try:
        source_routes = plan['routes']
        trail = plan['trail']
        assays = plan['assays']
    except (KeyError, TypeError):
        return dict(valid=False, errors=['plan requires routes, trail and assays'])
    for raw in source_routes:
        try:
            routes.append(Route(**raw) if isinstance(raw, dict) else raw)
        except (TypeError, ValueError):
            errors.append('malformed sortie')
    for arc in trail:
        if len(arc) != 4 or any(type(x) is not int for x in arc):
            errors.append('malformed truck arc')
            continue
        p, t, q, u = arc
        if not (0 <= p < len(inst.pads) and 0 <= q < len(inst.pads) and 0 <= t < u <= inst.H):
            errors.append('truck arc domain')
            continue
        if (p, t) != last:
            errors.append('truck path continuity')
        nodes.add((p, t))
        nodes.add((q, u))
        if p == q:
            if u != t+1:
                errors.append('waiting arcs must last one slot')
            for slot in range(t, u):
                stationary[slot] = p
        else:
            if u-t != inst.road(p, q):
                errors.append('road duration')
            road_cost += (u-t)*inst.road_cost
        last = (q, u)
    if last != (0, inst.H):
        errors.append('truck must end at depot at H')
    covered = []
    expected = Counter()
    good_routes = []
    for route in routes:
        if not isinstance(route, Route):
            errors.append('malformed sortie')
            continue
        p, launch, q, recovery, j, sample = (route.p, route.launch, route.q,
                                           route.recovery, route.job, route.sample)
        if any(type(x) is not int for x in (p, launch, q, recovery, j, sample)):
            errors.append('sortie clocks and indices must be integers')
            continue
        if not (0 <= j < len(inst.jobs) and 0 <= p < len(inst.pads) and 0 <= q < len(inst.pads)
                and 0 <= launch < recovery <= inst.H):
            errors.append('sortie domain')
            continue
        job = inst.jobs[j]
        covered.append(j)
        expected[(j, sample, recovery)] += 1
        good_routes.append(route)
        if not job.lo <= sample <= job.hi:
            errors.append('sampling window')
        if launch+inst.fly(p, len(inst.pads)+j) > sample:
            errors.append('outward flight reservation')
        if sample+inst.sampling_duration+inst.fly(len(inst.pads)+j, q) > recovery:
            errors.append('sampling and return reservation')
        if recovery-launch > inst.endurance:
            errors.append('flight endurance')
        if (p, launch) not in nodes or (q, recovery) not in nodes:
            errors.append('sortie truck rendezvous endpoint')
        end_ground = recovery+inst.ground_duration
        if end_ground > inst.H:
            errors.append('final ground work exceeds horizon')
        for slot in range(launch, min(end_ground, inst.H)):
            drone_use[slot] += 1
        for slot in range(recovery, min(end_ground, inst.H)):
            if stationary[slot] != q:
                errors.append('truck must remain stationary for recovery ground work')
                break
        flight_cost += inst.flight_fixed+inst.air_cost*(recovery-launch)
    if len(covered) != len(set(covered)):
        errors.append('duplicate sampled task')
    if require_all and set(covered) != set(range(len(inst.jobs))):
        errors.append('exact task coverage')
    if len(routes) > inst.batteries:
        errors.append('finite precharged pack inventory exceeded')
    if any(value > inst.drones for value in drone_use):
        errors.append('UAV occupancy includes recovery ground work')
    actual = Counter()
    for assignment in assays:
        if len(assignment) != 4 or any(type(x) is not int for x in assignment):
            errors.append('malformed assay start')
            continue
        j, sample, recovery, start = assignment
        if not 0 <= j < len(inst.jobs):
            errors.append('assay task domain')
            continue
        actual[(j, sample, recovery)] += 1
        job = inst.jobs[j]
        finish = start+inst.assay_duration
        if not (recovery+inst.handoff_delay <= start <= sample+job.shelf and finish <= min(job.deadline, inst.H)
                and 0 <= start < finish):
            errors.append('assay release, shelf or completion deadline')
        pads = set()
        for slot in range(max(0, start), min(finish, inst.H)):
            station_use[slot] += 1
            pads.add(stationary[slot])
        if None in pads or len(pads) != 1:
            errors.append('continuous stationary assay reservation')
        age_cost += job.weight*(finish-sample)
    if actual != expected:
        errors.append('one assay start for every returned sample')
    if any(value > inst.stations for value in station_use):
        errors.append('assay station capacity')
    return dict(valid=not errors, errors=list(dict.fromkeys(errors)),
                objective=road_cost+flight_cost+age_cost, road_cost=road_cost,
                flight_cost=flight_cost, age_cost=age_cost, covered=tuple(sorted(covered)),
                max_drone_use=max(drone_use, default=0), max_station_use=max(station_use, default=0),
                packs_used=len(routes))


def decode(inst, routes, actions, required=(), forbidden=()):
    """Match anonymous starts to concrete returned samples, then replay."""
    selected = [routes[rid] for action in actions for rid in action['launches']]
    returns = {}
    for route in selected:
        returns.setdefault(route.recovery, []).append(route)
    pending = []
    trail = []
    assays = []
    for action in actions:
        t = action['t']
        for route in returns.get(t, ()):
            j, s, e = route.job, route.sample, route.recovery
            job = inst.jobs[j]
            latest = min(inst.H-inst.assay_duration, s+job.shelf, job.deadline-inst.assay_duration)
            pending.append((job.weight, latest, j, s, e))
        for token in action['starts']:
            matching = [item for item in pending if item[:2] == token and item[4]+inst.handoff_delay <= t]
            if not matching:
                raise ValueError('anonymous start has no concrete returned sample')
            chosen = min(matching)
            pending.remove(chosen)
            w, latest, j, s, e = chosen
            assays.append((j, s, e, t))
        trail.append((action['p'], t, action['q'], action['u']))
    plan = dict(routes=selected, trail=trail, assays=assays)
    check = validate(inst, plan, require_all=False)
    covered = set(check.get('covered', ()))
    if pending or not set(required) <= covered or covered & set(forbidden) or not check['valid']:
        raise ValueError(f'invalid decoded physical plan: {check["errors"]}')
    plan['check'] = check
    return plan


def validate_fleet(inst, group_plans):
    """Independent per-group resources and global exactly-once task coverage."""
    errors=[];checks=[];covered=[]
    if len(group_plans)!=inst.groups:errors.append('fleet group cardinality')
    for plan in group_plans:
        c=validate(inst,plan,require_all=False);checks.append(c)
        errors.extend(c['errors']);covered.extend(c.get('covered',()))
    if Counter(covered)!=Counter(range(len(inst.jobs))):errors.append('global exact task coverage')
    return dict(valid=not errors,errors=list(dict.fromkeys(errors)),group_checks=checks,
                objective=sum((c['objective'] for c in checks),Fraction(0)),
                road_cost=sum((c['road_cost'] for c in checks),Fraction(0)),
                flight_cost=sum((c['flight_cost'] for c in checks),Fraction(0)),
                age_cost=sum((c['age_cost'] for c in checks),Fraction(0)),
                covered=tuple(sorted(covered)),packs_used=[c['packs_used'] for c in checks])
