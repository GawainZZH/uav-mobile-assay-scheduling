# Staggered verification calendars

This exact, four-task example studies group-dependent verification calendars.
It is separate from the 24-input SFBC performance comparison and requires no
numerical optimizer.

```console
python -m unittest tests.test_calendar41 -v
python -B -m calendar41.mechanism --output new-runs/staggered.json
```

Two pads are at (0,0) and (5,0). Two earlier tasks are at (0,1), with sampling
epoch 1 and assay deadline 5; two subsequent tasks are at (5,1), with epoch 16
and deadline 20. Each group has one UAV, two precharged packs and one station.
Endurance is three slots, assays take two slots, ground service takes three
slots and handoff delay is zero. Travel reservations follow unit Manhattan
duration. The horizon is 30. All values are explicit in `instance()`.

Both groups reserve verification [5,6) and [24,25). Additional reservations
are [13,14) for the first group and [8,10) for the second. Assays cannot overlap
verification; the truck must be parked during each reservation. A completed
assay is released at the end of the first subsequent verification interval.

Geometry and endurance force the earlier sorties A0--A3 and subsequent sorties
B15--B18. Assay starts are forced to 3 and 18. Each group must serve one task
from each phase. The four responsibility assignments and five possible outward
road starts per group produce exactly 100 candidates, covering the whole
instance rather than a prescribed favorable branch. The final road is forced
to [25,30).

Exactly 12 fleets are feasible. They normalize to four fleets at cut 8, with
roads [8,13) and [10,15). Every fleet has transportation cost 40, age cost 16
and weighted verified-release sum 62. The complete physical group-column DW
value is also 56: columns cover at most two tasks, so coverage of four tasks
forces every positive column to cover two; each such column costs 28.

The largest common free road gap is only three slots, although each individual
calendar has a five-slot gap. Normalization therefore can use different gaps
for different groups. At cut 8, both groups have state (A, one remaining pack).
Exchanging suffixes preserves basic physical feasibility and task coverage but
violates both verification calendars. Group labels preserve the boundary
representation; anonymous pooling requires identical continuation domains.

`data/staggered-calendar/results.json` preserves the complete feasible fleets,
normalizations and exchange counterexample. The code replays the original
physical validator and independently checks calendar occupation and release.
This example is not an implementation or performance evaluation of SFBC with
heterogeneous calendars.
