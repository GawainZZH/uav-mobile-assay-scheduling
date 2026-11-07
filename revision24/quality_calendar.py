"""Prescribed successful analytical-QC blocks and verified result release.

Every block closes ALL stations of EVERY group and requires parked trucks.
The blocks represent reserved IPC/blank/batch-preparation work, not measured
method durations. Calibration failure and repeat analyses are outside this
deterministic model. Original physical validation remains mandatory separately.
"""
from collections import Counter
from fractions import Fraction


def _blocks(blocks):
    result = tuple(tuple(block) for block in blocks)
    if not result or any(len(b) != 2 or any(type(v) is not int for v in b)
                         or b[0] < 0 or b[1] <= b[0] for b in result):
        raise ValueError('ordered positive integer QC blocks required')
    if any(a[1] > b[0] for a, b in zip(result, result[1:])):
        raise ValueError('overlapping or unordered QC blocks')
    return result


def _rounds(inst, blocks):
    blocks = _blocks(blocks)
    if len(blocks) % 3 or blocks[0][0] != 0 or blocks[-1][1] > inst.H:
        raise ValueError('three QC blocks per complete round required')
    starts = [blocks[i][0] for i in range(0, len(blocks), 3)]
    ends = starts[1:] + [inst.H]
    rounds = tuple((a, end, blocks[3*i+2][0]) for i, (a, end) in enumerate(zip(starts, ends)))
    if any(blocks[3*i+2][1] > end for i, (_, end, _) in enumerate(rounds)):
        raise ValueError('QC work crosses round boundary')
    return rounds


def calendar(inst, duration, period=60):
    """Return synchronized blocks; each whole round has at most eight tasks."""
    if type(duration) is not int or duration not in (2, 4):
        raise ValueError('QC duration must be 2 or 4')
    if type(period) is not int or period < 52+duration or inst.H % period:
        raise ValueError('horizon must contain whole QC rounds')
    totals = Counter(j.lo // period for j in inst.jobs)
    if max(totals.values(), default=0) > 8:
        raise ValueError('at most eight tasks per round')
    if any(j.hi >= (j.lo // period+1)*period for j in inst.jobs):
        raise ValueError('collection window crosses round boundary')
    return tuple((r+offset, r+offset+duration)
                 for r in range(0, inst.H, period) for offset in (0, 30, 52))


def release_epoch(completion, blocks):
    """End of the first verification block starting at/after completion."""
    return next((end for start, end in blocks if start >= completion), None)


def _round_for(job, rounds):
    return next(((start, end, last) for start, end, last in rounds
                 if start <= job.lo < end), None)


def eligibility(inst, j, s, e, k, S, blocks):
    """Assay feasibility, same-round closure, and window-opening release cap."""
    if any(type(v) is not int for v in (j, s, e, k)) or not 0 <= j < len(inst.jobs):
        return False
    job = inst.jobs[j]; finish = k+inst.assay_duration
    if not (job.lo <= s <= job.hi and e >= s+inst.sampling_duration
            and 0 <= k and e+inst.handoff_delay <= k <= s+job.shelf
            and finish <= min(job.deadline, inst.H)):
        return False
    if any(k < end and start < finish for start, end in blocks):
        return False
    bounds = _round_for(job, _rounds(inst, blocks))
    if bounds is None or job.hi >= bounds[1] or finish > bounds[2]:
        return False
    release = release_epoch(finish, blocks)
    return release is not None and release <= min(inst.H, job.lo+S, bounds[1])


def apply_model(m, blocks, S):
    """Mutate common matrix bounds/RHS; caller must rebuild its sparse matrix."""
    blocks = _blocks(blocks); rounds = _rounds(m.inst, blocks)
    road_disabled = assay_disabled = route_disabled = 0
    for (g, i), col in m.zids.items():
        p, t, q, u = m.arcs[i]
        if p != q and any(t < end and start < u for start, end in blocks):
            m.M.upper[col] = 0; road_disabled += 1
    for g in range(m.inst.groups):
        for start, end in blocks:
            for t in range(start, end):
                row = m.M.index['slot', g, t]
                m.M.rhs[row] = min(m.M.rhs[row], -m.inst.stations)
    for (g, j, s, e, k), col in m.assignments.items():
        if not eligibility(m.inst, j, s, e, k, S, blocks):
            m.M.upper[col] = 0; assay_disabled += 1
    for (g, i), col in getattr(m, 'rids', {}).items():
        route = m.routes[i]; bounds = _round_for(m.inst.jobs[route.job], rounds)
        if (bounds is None or route.launch < bounds[0]
                or route.recovery+m.inst.ground_duration > bounds[1]):
            m.M.upper[col] = 0; route_disabled += 1
    return dict(road_columns_disabled=road_disabled, assay_columns_disabled=assay_disabled,
                route_columns_disabled=route_disabled, quality_blocks_per_group=len(blocks))


def release_cost(m, blocks):
    """Map assignment column to rational weighted verified-completion cost."""
    result = {}
    for (g, j, s, e, k), col in m.assignments.items():
        release = release_epoch(k+m.inst.assay_duration, blocks)
        result[col] = m.inst.jobs[j].weight*release if release is not None else Fraction(0)
    return result


def validate_quality(inst, fleet, blocks, S):
    """Independently inspect QC intervals and release; no physical-cost replay."""
    errors = []; releases = [None]*len(inst.jobs); covered = Counter(); batch_sizes = Counter()
    try:
        blocks = _blocks(blocks); rounds = _rounds(inst, blocks)
        plans = fleet['group_plans'] if isinstance(fleet, dict) else fleet
        if len(plans) != inst.groups:
            errors.append('QC group total mismatch')
        by_round = Counter()
        for job in inst.jobs:
            bounds = _round_for(job, rounds)
            if bounds is None or job.hi >= bounds[1]:
                errors.append('collection window crosses round boundary')
            else:
                by_round[bounds[0]] += 1
        if max(by_round.values(), default=0) > 8:
            errors.append('more than eight tasks in QC round')
        for g, plan in enumerate(plans):
            parks = [0]*inst.H; occupation = [0]*inst.H
            for p, t, q, u in plan['trail']:
                if not 0 <= t < u <= inst.H:
                    errors.append('QC truck interval outside horizon'); continue
                if p == q:
                    for slot in range(t, u): parks[slot] += 1
                elif any(t < end and start < u for start, end in blocks):
                    errors.append('road overlaps quality block')
            for start, end in blocks:
                if any(parks[t] != 1 for t in range(start, end)):
                    errors.append('quality block requires continuous parking')
            for j, s, e, k in plan['assays']:
                if not 0 <= j < len(inst.jobs):
                    errors.append('assay task outside domain'); continue
                covered[j] += 1
                if not eligibility(inst, j, s, e, k, S, blocks):
                    errors.append('assay violates QC, original limits, or release promise')
                finish = k+inst.assay_duration
                for t in range(max(0, k), min(inst.H, finish)): occupation[t] += 1
                release = release_epoch(finish, blocks); releases[j] = release
                if release is not None: batch_sizes[g, release] += 1
            for t, number in enumerate(occupation):
                busy = any(start <= t < end for start, end in blocks)
                if number > (0 if busy else inst.stations):
                    errors.append('assay station capacity including QC exceeded')
            for route in plan['routes']:
                def value(name):
                    return route[name] if isinstance(route, dict) else getattr(route, name)
                j = value('job'); bounds = _round_for(inst.jobs[j], rounds)
                if (bounds is None or value('launch') < bounds[0]
                        or value('recovery')+inst.ground_duration > bounds[1]):
                    errors.append('flight or ground obligation crosses round boundary')
        if any(size > 10 for size in batch_sizes.values()):
            errors.append('over ten samples between consecutive quality blocks')
        if covered != Counter(range(len(inst.jobs))):
            errors.append('one quality-released assay required per task')
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        errors.append('malformed quality calendar: '+str(exc))
    weighted = sum((job.weight*r for job, r in zip(inst.jobs, releases) if r is not None), Fraction(0))
    return dict(valid=not errors, errors=list(dict.fromkeys(errors)), release_epochs=releases,
                weighted_release=weighted, total_quality_blocks=len(blocks)*inst.groups,
                quality_blocks_per_group=len(blocks),
                batch_sizes=[dict(group=g, release=release, samples=size)
                             for (g, release), size in sorted(batch_sizes.items())],
                physical_validation_included=False)
