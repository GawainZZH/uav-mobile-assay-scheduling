# Archived scientific studies

The [r41-legacy-evidence release](https://github.com/GawainZZH/uav-mobile-assay-scheduling/releases/tag/r41-legacy-evidence)
holds 28 content-addressed ZIP parts (572,621,768 compressed bytes in total).
The checked catalog describes 74,514 logical scientific files and 71,233 distinct
objects. Repeated files share an object. The data include original inputs,
candidate pools, complete saved native outputs, receipts, phase traces, scientific
summaries, source snapshots and the 17 assets used by the current 11 figure groups.
The current four tables are also mapped in `paper-map.json`.

The 28 cohorts are indexed in `collections.json`; per-cohort frozen Python source
maps are embedded in `catalog.json.gz`. A cohort may include several methods and
budgets. Retain the original configuration and arm identifiers when recomputing
statistics. Interrupted, infeasible and unresolved outcomes remain in the data.
The two later 72-arm studies are included in full, rather than only as plotted
summaries. The new stock-flow 48-arm comparison remains in the main repository's
`data/` directory and is not counted again in this archive.

## Read or restore

From the repository root, using Python's standard library:

```console
python -B legacy_evidence.py list
python -B legacy_evidence.py fetch --collection revision39_coordination_20261106
python -B legacy_evidence.py verify --collection revision39_coordination_20261106
python -B legacy_evidence.py materialize --collection revision39_coordination_20261106 --output restored/coordination
```

Use an exact collection identifier from `collections.json`. `valuation_reference`
is an alias for `fleet_sample_assay_second_review_20261016`. Omitting the collection
from `fetch` or `verify` selects the full archive. Parts contain multiple cohorts;
fetching one cohort can download shared parts. Default cache: `legacy-cache/`.
The commands verify the repository inventory, part hashes, exact archive member
population, object hashes and byte lengths. Restoration requires a new destination;
it adds the shared scientific Python files, then applies the selected cohort's
frozen versions. No numerical optimizer is called by these commands.

## Scientific use

`paper-map.json` identifies the saved evidence for each current figure and table.
Scientific plotting functions and their dependencies are retained with the source
files. Original plotting controllers can require local fonts, TeX and private
publication metadata; importing the relevant scientific function or reading the
saved numerical source avoids that metadata requirement. The saved figure PDFs
are also provided. A plot's visual values do not replace its raw study outputs.

The restored `reproducibility/reproduce_sample_assay.py` and
`reproducibility/reproduce_second_review.py` provide portable one-instance entries
for the corresponding archived model families. Run `--help` inside the restored
directory for the required arguments. Native solving requires the dependencies
reported with that cohort. These entries do not reproduce the process manager's
historical timing or replace the original measurements. No legacy solve was rerun
to create this release. The stock-flow algorithm has a separate executable entry,
`reproduce.py`, with its own input regeneration and result replay.

Some original batch controllers validate a private admission chain before running.
Those controllers are retained as scientific source history, but the admission
chain is not a requirement for the public restoration or one-instance entries.
Model source versions and numerical inputs are identified by SHA256 in the catalog.
When a machine-specific path was anonymized, `original_sha256` identifies the
original bytes and `sha256` identifies the delivered bytes. Numerical fields were
preserved. Internal reviewer correspondence, manuscript text, credentials and
solver licenses are excluded.

The main paper distinguishes the unit-service reference, non-unit parked assays
and verified-release policies. A restored source version belongs to its named
study; do not pool different objective definitions or substitute an older method
for the stock-flow comparison. The archive supports inspection and recomputation
from saved evidence. Byte verification does not independently prove a historical
solver's global optimality conclusion.
