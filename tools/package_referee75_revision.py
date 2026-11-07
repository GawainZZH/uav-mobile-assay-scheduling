"""Prepare or build additive Referee75 packages; import and --prepare never write.

No optimizer, experiment runner, manuscript builder, or verifier is executed.
Creation requires final, hash-bound verification supplied by the revision owner.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import zipfile


DATE = '20261018'
PREFIX = 'TS_Referee75_'
PAPER = 'manuscript/sample-assay-ts'
BASELINE = 'submission_packages/TS_Round3_Major78_Research_Reproduction_20261017.zip'
BASELINE_SHA256 = '424a14892f230ce3aa5f83a915ecbbd290a6d3e67b75f5bb0134cfdafbbc5f60'
FINAL_VERIFICATION = 'docs/2026-10-18-referee75-final-verification.json'
RESPONSE = 'docs/2026-10-18-referee75-response.md'
MECHANISM_STATUS = 'docs/data/referee75_20261018/mechanism-final-status.json'
FROZEN_MANIFEST = 'results/algorithm_recovery_closure_validation_20261017/manifest.json'
ARCHIVE_MANIFEST = 'REFEREE75_PACKAGE_MANIFEST.json'
TEMPLATE_NAMES = {'informs4.cls', 'eqndefns-left.sty', 'informs2014trsc.bst', 'informs_Logo.pdf'}
COHORT_NAMES = {
    'algorithm_recovery_closure_pilot_20261017',
    'algorithm_recovery_closure_validation_20261017',
    'algorithm_recovery_pooled_pilot_20261017',
    'algorithm_recovery_root_pilot_20261017',
    'algorithm_recovery_validation_20261017',
}
MECHANISM_EVIDENCE_DIRS = ('referee75_mechanism/freeze-preview',
                           'referee75_mechanism/formal-run')


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def encoded(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode('utf8')


def read_json(path):
    return json.loads(path.read_text(encoding='utf8'))


def safe_path(root, name):
    """Confine every newly selected file to the explicitly supplied project."""
    root = Path(root).resolve()
    if '\\' in str(name) or PurePosixPath(str(name)).is_absolute():
        raise ValueError(f'unsafe path: {name}')
    target = (root / name).resolve()
    if not target.is_relative_to(root):
        raise ValueError(f'path escapes outside project: {name}')
    return target


def tex_dependencies(paper, entry='main.tex'):
    """Return recursive literal TeX inputs and graphics, without build products.

    Resolution supports project-relative paths, paths relative to the including
    source, and literal graphicspath directories. Dynamic TeX paths fail closed.
    """
    paper = Path(paper).resolve()
    selected = set()
    pending = [safe_path(paper, entry)]
    graphic_dirs = []

    def resolve(item, parent, extensions, graphic=False):
        if re.search(r'[\\{}#$]', item):
            raise ValueError(f'dynamic TeX dependency needs explicit resolution: {item}')
        safe_path(paper, item)
        candidates = [parent / item, paper / item]
        if graphic:
            candidates += [folder / item for folder in graphic_dirs]
            candidates.append(paper / 'figures' / item)
        for candidate in candidates:
            for suffix in ('',) if candidate.suffix else extensions:
                path = Path(str(candidate) + suffix).resolve()
                if not path.is_relative_to(paper):
                    raise ValueError(f'TeX dependency escapes outside project: {item}')
                if path.is_file():
                    return path
        raise FileNotFoundError(f'TeX dependency missing: {item}')

    while pending:
        path = pending.pop()
        name = path.relative_to(paper).as_posix()
        if name in selected:
            continue
        if not path.is_file():
            raise FileNotFoundError(path)
        selected.add(name)
        text = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text(encoding='utf8'))
        for block in re.findall(r'\\graphicspath\s*\{((?:\s*\{[^{}]+\})+)\s*\}', text):
            for item in re.findall(r'\{([^{}]+)\}', block):
                graphic_dirs.append(safe_path(paper, item))
        for item in re.findall(r'\\(?:input|include)\s*\{([^{}]+)\}', text):
            pending.append(resolve(item, path.parent, ('.tex',)))
        for item in re.findall(r'\\includegraphics\*?(?:\s*\[[^\]]*\])?\s*\{([^{}]+)\}', text):
            found = resolve(item, path.parent, ('.pdf', '.png', '.jpg', '.jpeg', '.eps'), True)
            selected.add(found.relative_to(paper).as_posix())
        for items in re.findall(r'\\(?:bibliography|addbibresource)\s*\{([^{}]+)\}', text):
            for item in items.split(','):
                found = resolve(item.strip(), path.parent, ('.bib',))
                selected.add(found.relative_to(paper).as_posix())
    return selected


def verified_template_assets(paper):
    paper = Path(paper)
    expected = read_json(paper / 'template-assets.json')
    if set(expected) != TEMPLATE_NAMES:
        raise ValueError('official template asset dictionary has unexpected names')
    for name, expected_hash in expected.items():
        if digest(safe_path(paper, name).read_bytes()) != expected_hash:
            raise ValueError(f'official template hash mismatch: {name}')
    return expected


def write_preserved_archive(baseline, destination, overlay, metadata=None):
    """Preserve every prior file byte, relocating only changed-name conflicts.

    The returned root manifest hashes all payloads, excluding itself. Historical
    cache/binary files remain intact. ZIP compression bytes are not preserved;
    each original member's uncompressed bytes and mapping are preserved.
    """
    destination = Path(destination)
    if destination.exists():
        raise FileExistsError(f'refuse replacement: {destination}')
    overlay = dict(overlay)
    for name, raw in overlay.items():
        safe_path(destination.parent, name)
        if not isinstance(raw, bytes):
            raise TypeError(f'overlay content must be bytes: {name}')
    inventory = {}
    preserved = {}
    with zipfile.ZipFile(baseline) as prior:
        items = prior.infolist()
        names = [item.filename for item in items]
        if len(names) != len(set(names)):
            raise ValueError('baseline contains duplicate member names')
        occupied = set(names) | set(overlay) | {ARCHIVE_MANIFEST}

        def historical(name):
            base = 'history/round3_archive_20261017/' + name
            result = base
            serial = 1
            while result in occupied:
                result = f'{base}.preserved-{serial}'
                serial += 1
            occupied.add(result)
            return result

        with zipfile.ZipFile(destination, 'x', zipfile.ZIP_DEFLATED, compresslevel=7) as out:
            def add(name, raw):
                if name in inventory:
                    raise ValueError(f'duplicate output member: {name}')
                out.writestr(name, raw)
                inventory[name] = dict(bytes=len(raw), sha256=digest(raw))

            for item in items:
                name = item.filename
                raw = prior.read(item)
                conflict = name == ARCHIVE_MANIFEST or (name in overlay and overlay[name] != raw)
                location = historical(name) if conflict else name
                add(location, raw)
                preserved[name] = dict(location=location, bytes=len(raw), sha256=digest(raw))
                if not conflict:
                    overlay.pop(name, None)
            for name, raw in sorted(overlay.items()):
                if name == ARCHIVE_MANIFEST:
                    raise ValueError('root package manifest is reserved')
                add(name, raw)
            manifest = dict(metadata or {})
            manifest.update(schema='referee75-package-v1', date='2026-10-18', files=inventory,
                            prior_archive=dict(sha256=digest(Path(baseline).read_bytes()),
                                               files=len(preserved), original_byte_mapping=preserved),
                            file_count_including_this_manifest=len(inventory) + 1)
            out.writestr(ARCHIVE_MANIFEST, encoded(manifest))
    return manifest


def _local_import_closure(root, selected):
    """Include available local Python dependencies without importing/executing."""
    pending = [name for name in selected if name.endswith('.py')]
    visited = set()
    while pending:
        name = pending.pop()
        if name in visited:
            continue
        visited.add(name)
        tree = ast.parse(safe_path(root, name).read_text(encoding='utf-8-sig'), filename=name)
        for node in ast.walk(tree):
            modules = [a.name for a in node.names] if isinstance(node, ast.Import) else []
            if isinstance(node, ast.ImportFrom) and node.module and not node.level:
                modules.append(node.module)
            for module in modules:
                fragment = module.replace('.', '/')
                for candidate in (fragment + '.py', fragment + '/__init__.py',
                                  'tools/' + fragment + '.py', 'docs/' + fragment + '.py'):
                    if safe_path(root, candidate).is_file() and candidate not in selected:
                        selected.add(candidate)
                        pending.append(candidate)
                # Package initializers are also required for dotted imports.
                parts = fragment.split('/')
                for count in range(1, len(parts)):
                    init = '/'.join(parts[:count]) + '/__init__.py'
                    if safe_path(root, init).is_file() and init not in selected:
                        selected.add(init)
                        pending.append(init)


def prepare_manifest(root, verification=FINAL_VERIFICATION, mechanism_status=MECHANISM_STATUS):
    """Read-only inventory and readiness report. Missing evidence is a blocker."""
    root = Path(root).resolve()
    blockers = []
    selected = set()
    source = set()
    final = {}
    frozen_sources = {}

    def require(name):
        if not safe_path(root, name).is_file():
            blockers.append(f'missing required file: {name}')
            return False
        selected.add(name)
        return True

    for name in (BASELINE, verification, mechanism_status, RESPONSE, FROZEN_MANIFEST,
                 PAPER + '/main.pdf', 'reproducibility/requirements.txt'):
        require(name)
    # The baseline is preserved through the archive merge, never nested/copied.
    selected.discard(BASELINE)
    if (root / BASELINE).is_file() and digest((root / BASELINE).read_bytes()) != BASELINE_SHA256:
        blockers.append('baseline archive hash differs from the inspected Round3 archive')
    try:
        source = tex_dependencies(root / PAPER)
        source.update(verified_template_assets(root / PAPER))
        source.update(('template-assets.json', 'BUILD.ps1'))
        for name in source:
            require(PAPER + '/' + name)
    except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as exc:
        blockers.append(f'LaTeX source closure: {exc}')
    if (root / FROZEN_MANIFEST).is_file():
        frozen_sources = read_json(root / FROZEN_MANIFEST)['sources']
        if len(frozen_sources) != 209:
            blockers.append('frozen algorithm source inventory must contain 209 files')
        for name, expected in frozen_sources.items():
            if require(name) and digest(safe_path(root, name).read_bytes()) != expected:
                blockers.append(f'frozen algorithm source hash mismatch: {name}')
    found_cohorts = {p.name for p in (root / 'results').glob('algorithm_recovery*') if p.is_dir()}
    if found_cohorts != COHORT_NAMES:
        blockers.append('fixed cohort directory set differs: missing=' +
                        repr(sorted(COHORT_NAMES - found_cohorts)) + ', unexpected=' +
                        repr(sorted(found_cohorts - COHORT_NAMES)))
    cohorts = [root / 'results' / name for name in sorted(COHORT_NAMES)]
    cohort_files = {p.relative_to(root).as_posix() for folder in cohorts
                    for p in folder.rglob('*') if p.is_file()}
    evidence_files = cohort_files | {RESPONSE, mechanism_status}
    for name in MECHANISM_EVIDENCE_DIRS:
        folder = root / name
        if not folder.is_dir():
            blockers.append(f'missing mechanism evidence directory: {name}')
        records = {p.relative_to(root).as_posix() for p in folder.rglob('*') if p.is_file()}
        selected.update(records)
        evidence_files.update(records)
    for folder in sorted((root / 'results').glob('referee75_mechanism*')) + [root / 'docs/data/referee75_20261018']:
        if folder.is_dir():
            evidence_files.update(p.relative_to(root).as_posix() for p in folder.rglob('*') if p.is_file())
    # Preserve every raw/process/stdout/stderr/canonical file, without suffix filters.
    for folder in cohorts + sorted((root / 'results').glob('referee75_mechanism*')):
        if folder.is_dir():
            selected.update(p.relative_to(root).as_posix() for p in folder.rglob('*') if p.is_file())
    patterns = (
        'algorithm_recovery*/**/*.py', 'referee75_mechanism/**/*.py',
        'docs/data/algorithm-recovery/**/*', 'docs/data/referee75_20261018/**/*',
        'docs/2026-10-17-*recovery*.*', 'docs/2026-10-17-root-closure*.*',
        'docs/2026-10-18-referee75*.*', 'docs/plans/*referee75*.md',
        'docs/plans/*algorithm-recovery*.md', 'tools/*algorithm_recovery*.py',
        'tools/*recovery_closure*.py', 'tools/referee75*.py',
        'tools/package_referee75_revision.py', 'tools/verify_referee75_revision.py',
        'tools/verify_referee75_final.py', 'docs/referee75_failure_archive_audit.py',
        'tools/verify_pooled_manuscript.py', 'tools/reproduce_referee75_example.py',
        'tests/*algorithm_recovery*.py', 'tests/*recovery_closure*.py',
        'tests/*recovery_parent*.py', 'tests/*referee75*.py',
    )
    for pattern in patterns:
        selected.update(p.relative_to(root).as_posix() for p in root.glob(pattern) if p.is_file())
    # Only bounded, nonsensitive manuscript snapshots are new history overlays.
    for snapshot in ('before-referee75_20261018', 'before-pooled-rewrite_20261018'):
        for name in ('main.tex', 'main.pdf', 'pooled-study.tex', 'pooled-method.tex',
                     'appendix-algorithm.tex', 'appendix-model.tex', 'appendix-proofs.tex',
                     'appendix-results.tex'):
            path = root / PAPER / 'history' / snapshot / name
            if path.is_file():
                selected.add(path.relative_to(root).as_posix())
    _local_import_closure(root, selected)
    if (root / verification).is_file():
        final = read_json(root / verification)
        if final.get('status') != 'PASS':
            blockers.append('final-verification status must be PASS')
        hashes = final.get('hashes', {})
        for name in source | {'main.pdf'}:
            relative = PAPER + '/' + name
            expected = hashes.get(relative, hashes.get(name))
            if expected != digest(safe_path(root, relative).read_bytes()):
                blockers.append(f'final-verification missing or stale hash: {relative}')
        research_inputs = final.get('research_inputs_inventory', {})
        if not isinstance(research_inputs, dict) or set(research_inputs) != cohort_files:
            blockers.append('research input inventory differs from complete fixed-cohort file set')
        extra_hashes = final.get('extra_evidence_hashes', {})
        if not isinstance(extra_hashes, dict):
            blockers.append('extra_evidence_hashes must be a root-relative hash dictionary')
            extra_hashes = {}
        for name in sorted(evidence_files - set(extra_hashes)):
            blockers.append(f'final-verification missing evidence hash: {name}')
        # Validate every recorded entry, including a deleted file that scanning
        # alone could no longer discover. Neither a failed arm nor its raw log
        # can disappear between the final audit and packaging.
        for label, hashes in (('research input inventory', research_inputs),
                              ('evidence', extra_hashes)):
            if not isinstance(hashes, dict):
                continue
            for name, record in hashes.items():
                expected = record.get('sha256') if isinstance(record, dict) else record
                path = safe_path(root, name)
                if not path.is_file():
                    blockers.append(f'missing evidence from {label}: {name}')
                elif expected != digest(path.read_bytes()):
                    blockers.append(f'stale {label} hash: {name}')
    if (root / mechanism_status).is_file():
        mechanism = read_json(root / mechanism_status)
        if mechanism.get('status') != 'PASS':
            blockers.append('mechanism final status must be PASS (including censored/negative results)')
    if (root / RESPONSE).is_file() and not (root / RESPONSE).read_text(encoding='utf8').strip():
        blockers.append('English referee response is empty')
    inventory = {}
    for name in sorted(selected):
        raw = safe_path(root, name).read_bytes()
        inventory[name] = dict(bytes=len(raw), sha256=digest(raw))
    return dict(status='NOT_READY' if blockers else 'READY', blockers=blockers,
                schema='referee75-package-plan-v1', source_files=sorted(source),
                selected_files=inventory, frozen_209_sources=frozen_sources,
                verification=verification, mechanism_status=mechanism_status,
                pages=final.get('pages'), output_prefix=PREFIX)


SOURCE_README = '''# Referee75 Transportation Science source

Run BUILD.ps1 here (pdfLaTeX, three passes), or run pdflatex with
-interaction=nonstopmode -halt-on-error -file-line-error main.tex three times.
The class, equation style, bibliography style and logo are verified against
template-assets.json. This ZIP contains the actual recursive input/graphics
closure and build instructions; it contains no completed manuscript PDF.
The separate manuscript, English response and full reproduction package match
this revision. Source hashes are in source-package-manifest.json.
'''

RESEARCH_README = '''# Referee75 research reproduction

The active manuscript is manuscript/sample-assay-ts/main.tex and main.pdf.
The matching English response is docs/2026-10-18-referee75-response.md.
REFEREE75_PACKAGE_MANIFEST.json maps EVERY original Round3 archive member to
its exact preserved bytes; changed-name conflicts live under
history/round3_archive_20261017/. Historical listed caches are retained.
The root manifest hashes all payload files except itself and records the frozen
209-file algorithm inventory, final verification, and companion ZIP hashes.

Install reproducibility/requirements.txt in a separate Python environment.
Commands below are offline audits of saved evidence. They do not repeat the
timed optimization experiments. Run from this extracted archive root:

    python tools/analyze_algorithm_recovery.py --root . --folder results/algorithm_recovery_closure_validation_20261017 --out closure-reanalysis.json
    python tools/replay_algorithm_recovery_bounds.py --root . --folder results/algorithm_recovery_closure_validation_20261017 --out closure-dag-replay.json
    python tools/replay_recovery_closure.py --folder results/algorithm_recovery_closure_validation_20261017 --out closure-box-replay.json

The last two commands reconstruct the frozen graphs and may require memory and
time; they make zero optimizer calls. The complete five algorithm-recovery
cohorts include failures, raw results, process measurements and stdout/stderr.
The Referee75 mechanism cohort, including censored/negative outcomes, is retained
separately and does not replace the algorithm performance cohort. Its manifest
and saved final audit define its scope. A PASS archive audit certifies retained
records, not experiment success; inspect experiment_status and coverage counts.
No new experiment is launched by this
packager. Historical writer/build/package scripts are archival evidence and
MUST NOT be used to regenerate the active revision. Build only the current
LaTeX closure using its BUILD.ps1 or the separate source ZIP.

tools/referee75_statistics_audit.py and tools/referee75_calendar_audit.py replay
saved evidence; they write their documented report outputs. The worked-example
script rebuilds the explicitly specified small example and writes
docs/data/referee75_20261018/worked-example.json. It is not an experiment runner.
Run such report writers only in a disposable extracted copy. The document
verifier requires pdflatex-generated aux/log plus pdfinfo; it is separate from
data-only replays. Author-specific unresolved facts remain author placeholders.
'''


def create_packages(root, verification=FINAL_VERIFICATION, mechanism_status=MECHANISM_STATUS, suffix=''):
    """Build new filenames only after all final gates pass; no implicit solving."""
    root = Path(root).resolve()
    plan = prepare_manifest(root, verification, mechanism_status)
    if plan['status'] != 'READY':
        raise ValueError('Referee75 packaging not ready: ' + '; '.join(plan['blockers']))
    submission = root / 'submission_packages'
    if suffix and not re.fullmatch(r'_[a-zA-Z0-9]+', suffix):
        raise ValueError('suffix must be an underscore followed by letters/digits')
    outputs = {kind: submission / f'{PREFIX}{label}_{DATE}{suffix}.{extension}' for kind, label, extension in (
        ('source', 'LaTeX', 'zip'), ('response', 'Response', 'zip'),
        ('research', 'Research_Reproduction', 'zip'), ('pdf', 'Manuscript', 'pdf'),
        ('delivery', 'Delivery', 'json'))}
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(f'refuse replacement: {path}')
    source = {name: safe_path(root / PAPER, name).read_bytes() for name in plan['source_files']}
    source['README-REFEREE75.md'] = SOURCE_README.encode()
    source['source-package-manifest.json'] = encoded({name: digest(raw) for name, raw in sorted(source.items())})
    with zipfile.ZipFile(outputs['source'], 'x', zipfile.ZIP_DEFLATED, compresslevel=7) as archive:
        for name, raw in sorted(source.items()):
            archive.writestr(name, raw)
    with zipfile.ZipFile(outputs['response'], 'x', zipfile.ZIP_DEFLATED, compresslevel=7) as archive:
        archive.writestr('Response_to_Reviewers.md', (root / RESPONSE).read_bytes())
    overlay = {name: safe_path(root, name).read_bytes() for name in plan['selected_files']}
    # Recheck the selected inventory immediately before writing its bytes.
    for name, raw in overlay.items():
        if digest(raw) != plan['selected_files'][name]['sha256']:
            raise ValueError(f'file changed during packaging: {name}')
    overlay['README-REFEREE75.md'] = RESEARCH_README.encode()
    review = dict(version='Referee75-20261018', status='final_verification_passed',
                  target='Transportation Science', verification=verification,
                  mechanism_final_status=mechanism_status, pages=plan['pages'],
                  mechanism_archive_audit=read_json(root / mechanism_status),
                  manuscript_sha256=digest((root / PAPER / 'main.tex').read_bytes()),
                  manuscript_pdf_sha256=digest((root / PAPER / 'main.pdf').read_bytes()))
    overlay[PAPER + '/review-referee75.json'] = encoded(review)
    manifest = write_preserved_archive(root / BASELINE, outputs['research'], overlay, dict(
        status='created_pending_extracted_verification', frozen_209_sources=plan['frozen_209_sources'],
        selected_inventory=plan['selected_files'], review=review,
        source_zip_sha256=digest(outputs['source'].read_bytes()),
        response_zip_sha256=digest(outputs['response'].read_bytes())))
    with outputs['pdf'].open('xb') as out:
        out.write((root / PAPER / 'main.pdf').read_bytes())
    delivery = dict(status='created_pending_extracted_verification',
                    prior_files_preserved=manifest['prior_archive']['files'],
                    research_files=manifest['file_count_including_this_manifest'],
                    source_files=len(source), pages=plan['pages'], artifacts={
                        kind: dict(path=str(path), bytes=path.stat().st_size,
                                   sha256=digest(path.read_bytes()))
                        for kind, path in outputs.items() if kind != 'delivery'})
    with outputs['delivery'].open('xb') as out:
        out.write(encoded(delivery))
    return delivery


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--verification', default=FINAL_VERIFICATION)
    parser.add_argument('--mechanism-status', default=MECHANISM_STATUS)
    parser.add_argument('--prepare', action='store_true', help='read-only readiness/inventory; no ZIPs')
    parser.add_argument('--suffix', default='', help='new version suffix, e.g. _v2; existing outputs are never replaced')
    args = parser.parse_args(argv)
    result = prepare_manifest(args.root, args.verification, args.mechanism_status) if args.prepare else create_packages(args.root, args.verification, args.mechanism_status, args.suffix)
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
    return 0 if result['status'] != 'NOT_READY' else 2


if __name__ == '__main__':
    raise SystemExit(main())
