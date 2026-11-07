"""Primal-only public prefixes and total charged workflow clocks."""
from copy import deepcopy
import hashlib
import json
import math
import time


def load_prefix(path, expected_sha, input_sha):
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=expected_sha:
        raise ValueError('public prefix hash mismatch')
    row=json.loads(data)
    if row['input_sha256']!=input_sha:
        raise ValueError('public prefix input mismatch')
    charge=row['charge_seconds'];source=row['source_elapsed_seconds'];cutoff=row['cutoff_seconds']
    if not all(math.isfinite(v) for v in (charge,source,cutoff)) or not 0<=cutoff<=source<=charge:
        raise ValueError('invalid public prefix charge or cutoff')
    fleets=[]
    for event in row['candidates']:
        available=event['available_seconds']
        if not math.isfinite(available) or not 0<=available<=cutoff:
            raise ValueError('candidate unavailable at prefix cutoff')
        # Original objectives, duals, trees and cached contexts are never accepted.
        fleets.append({'group_plans':deepcopy(event['fleet']['group_plans'])})
    return dict(charge_seconds=charge,fleets=fleets,prefix_sha256=expected_sha)


def remaining(total, origin, charge, *, now=None):
    now=time.perf_counter() if now is None else now
    if not all(math.isfinite(v) for v in (total,origin,charge,now)) or total<=0 or charge<0 or now<origin:
        raise ValueError('invalid charged clock')
    left=total-charge-(now-origin)
    if left<=0:raise TimeoutError('charged allowance exhausted')
    return left
