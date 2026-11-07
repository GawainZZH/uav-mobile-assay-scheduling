# UAV sampling and mobile assay scheduling

This research code implements homogeneous-group stock flow branch-and-check
(SFBC) and a jointly strengthened monolithic MILP. The algorithms coordinate
UAV collections, parked assays, vehicle rendezvous and finite precharged stocks.
Qualifying temporal boundaries allow block concatenation. Large merged blocks
can remain difficult for both methods.

This release contains the **24-input stock flow comparison** from the revised
manuscript. It is the new algorithm's reproducibility component. Older unit-service,
verified-release and operational comparisons are not bundled here. No manuscript,
internal referee report, solver license or historical approval chain is included.

## Environment

Recorded environment: Python 3.13.5, PySCIPOpt 6.2.1, NumPy 2.5.0,
SciPy 1.18.0, highspy 1.15.1 and psutil 7.1.3, on Windows 11 with an
Intel Core i9-14900KF. SCIP is the comparison solver; these results are not a
Gurobi comparison. The pinned Python packages supply the numerical dependencies.
No claim of runtime portability across operating systems follows from the archive.

```console
python -m venv .venv
# Activate .venv using your operating system's command, then:
python -m pip install -r requirements.txt
python -B reproduce.py verify
python -B reproduce.py regenerate
python -B reproduce.py summarize
python -B reproduce.py replay
```

`verify` validates the public file inventory. `regenerate` reconstructs all
24 raw inputs and common candidates without numerical optimization, comparing
them with the archived files. `summarize` recomputes the complete 48-arm results
from the original receipts, with fixed arm identities and result-to-receipt bindings.
`replay` verifies all 48 complete result hashes and physically validates the
corresponding fleet intervals. It does not re-solve an optimization problem.
Each operation validates file hashes before loading
scientific modules. Checksums detect drift, not malicious changes to a manifest
and all of its files. The separately published archive digest identifies a release.

## Run one instance

```console
python -B reproduce.py run --case DEMO --method F --seconds 20 --output new-runs/demo-F.json
python -B reproduce.py run --case DEMO --method S --seconds 20 --output new-runs/demo-S.json
python -B reproduce.py run --case RSC_A_s11 --method F --seconds 180 --output new-runs/A11-F.json
```

`F` denotes SFBC; `S` denotes the strengthened MILP. `DEMO` is a four-task
functional example, excluded from every performance statistic. Published case
names are listed in `data/fixed48-results.csv`. Every call cold-generates and
physically validates the common candidate. It does not import old prices,
incumbent vectors or search trees. An existing output file is never overwritten.

The convenience command reserves up to fifteen seconds of the requested internal
allowance for completion. It does not impose an operating-system process deadline
or monitor process-tree memory. Use an external process manager when enforcing
such caps. The original comparison used 180 seconds internally, a 165-second
native deadline, 210 seconds externally, one thread and 2 GiB. The public command
calls the same frozen scientific entry but replaces private admission and
publication machinery. New wall times must not replace the frozen observations.

## Results and interpretation

All 24 inputs and both methods remain in the denominator. Unproved runs have
PAR2 420 seconds. SFBC/MILP prove 18/4 optima, with PAR2 118.474088/365.739750
seconds, a 67.606997% decrease. Four-pad and merged-block failures remain visible.
Median sampled peak memory is 370.53/762.75 MiB; the maximum reverses this comparison
at 1663.29/1135.91 MiB. No uniform memory advantage is asserted.

`data/summaries/` contains the physical plans and intervals. `data/receipts/`
contains the original timing and process-tree memory observations. Complete
matrix, pricing and native result JSON files are preserved in `data/raw-results.zip`.
`data/analysis.json` and `data/fixed48-results.csv` preserve the reported statistics.
Rational price calculations and independent physical plan validation are distinct
from SCIP's numerical branch-and-bound conclusions. A finite upper cap for an
unverified local family is not a feasible fleet. Resource exhaustion remains an
unresolved outcome.

The original scientific Python files are copied byte for byte. `PROVENANCE.json`
maps those files to their frozen hashes. The small `reproduction40` adapter and
the root command are new release tooling. Some transitive modules contain legacy
helper functions; no historical benchmark command is invoked by this interface.
