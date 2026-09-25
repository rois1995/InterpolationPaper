#!/usr/bin/env python3
"""Build a flat Elsevier submission export without changing repository inputs.

Run after compiling paper.tex. Defaults to a strict author-information check.
Use --allow-incomplete only to prepare a review draft, never a final submission.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def pending_author_items() -> list[str]:
    pending = []
    patterns = ('to be supplied', '[Authors:', '[Affiliation', '[Author names', '[City]', '[Country]')
    for rel in ('submission/author_details.tex', 'submission/author_declarations.tex'):
        text = (ROOT / rel).read_text(encoding='utf-8')
        visible = '\n'.join(line.split('%', 1)[0] for line in text.splitlines())
        if any(token in visible for token in patterns):
            pending.append(rel)
    return pending

def get_dependencies(main: str) -> dict[str, str]:
    """Discover literal CSV inputs and the two author-supplied TeX files."""
    paths = set()
    for match in re.finditer(r'\{((?:\\DataDir/|\\GeneratedDir/|data/|generated/|RANS_PostProcessing/)[^{}]+\.csv)\}', main):
        value = match.group(1).replace('\\DataDir/', 'data/').replace('\\GeneratedDir/', 'generated/')
        paths.add(value)
    paths.update(('References.bib', 'paper.bbl', 'submission/author_details.tex',
                  'submission/author_declarations.tex', 'submission/figures/Figure_1.pdf',
                  'submission/figures/Figure_2.pdf', 'submission/figures/figure_captions.txt'))
    mapping: dict[str, str] = {}
    reverse: dict[str, str] = {}
    for rel in sorted(paths):
        src = (ROOT / rel).resolve()
        if not src.is_relative_to(ROOT) or not src.is_file():
            raise FileNotFoundError(f'Missing or unsafe manuscript input: {rel}')
        name = src.name
        if name in reverse and reverse[name] != rel:
            raise ValueError(f'Flattening collision: {reverse[name]} and {rel}')
        reverse[name] = rel
        mapping[rel] = name
    return mapping

def run_build(directory: Path) -> None:
    exe = shutil.which('latexmk')
    if exe is None:
        raise RuntimeError('latexmk is required; install a current TeX distribution.')
    command = [exe, '-g', '-pdf', '-interaction=nonstopmode', '-halt-on-error']
    if shutil.which('bibtex') is None:
        fallback = shutil.which('bibtex.original') or shutil.which('bibtex8')
        if fallback is None:
            raise RuntimeError('BibTeX is required.')
        command += ['-e', '$bibtex = "' + fallback.replace('\\', '/') + ' %O %B";']
    command.append('paper.tex')
    result = subprocess.run(command, cwd=directory, text=True, capture_output=True)
    (directory / 'export_build.log').write_text(result.stdout + result.stderr, encoding='utf-8')
    if result.returncode:
        raise RuntimeError('Submission build failed; inspect export_build.log in the temporary export.')
    log = (directory / 'paper.log').read_text(errors='replace')
    if re.search(r'undefined (?:references|citations)|Overfull \\[hv]box|! LaTeX Error', log):
        raise RuntimeError('Submission build has unresolved references/citations or overflowing content.')

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'dist'/'computers_fluids')
    parser.add_argument('--zip', dest='zip_path', type=Path, default=ROOT/'dist'/'ComputersFluids_submission_sources.zip')
    parser.add_argument('--allow-incomplete', action='store_true')
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--no-build', action='store_true', help='Copy source only; do not claim a tested submission build.')
    args = parser.parse_args()
    pending = pending_author_items()
    if pending:
        print('Author confirmation remains in: ' + ', '.join(pending), file=sys.stderr)
        if not args.allow_incomplete:
            print('Complete these fields before final export; use --allow-incomplete for a review draft.', file=sys.stderr)
            return 2
    if args.check_only:
        print('Author-information preflight passed.' if not pending else 'Review draft: author fields pending.')
        return 0
    main_tex = (ROOT/'paper.tex').read_text(encoding='utf-8')
    mapping = get_dependencies(main_tex)
    hashes = {rel: sha256(ROOT/rel) for rel in mapping}
    hashes['paper.tex'] = sha256(ROOT/'paper.tex')
    out = args.output_dir.resolve()
    if out == ROOT or out in (ROOT/'data', ROOT/'generated', ROOT/'RANS_PostProcessing', ROOT/'submission'):
        raise ValueError('Output must be a separate export directory, not a source directory.')
    if out.exists() and any(out.iterdir()):
        marker = out / '.cf_generated_export'
        if not marker.is_file() or marker.read_text().strip() != 'generated by prepare_cf_submission.py':
            raise FileExistsError(f'Refusing to replace an unrecognized nonempty directory: {out}')
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    (out / '.cf_generated_export').write_text('generated by prepare_cf_submission.py\n')
    for rel, name in mapping.items():
        shutil.copyfile(ROOT/rel, out/name)
    flat = main_tex.replace(r'\newcommand{\DataDir}{data}', r'\newcommand{\DataDir}{.}')
    flat = flat.replace(r'\newcommand{\GeneratedDir}{generated}', r'\newcommand{\GeneratedDir}{.}')
    for rel, name in mapping.items():
        flat = flat.replace('{'+rel+'}', '{'+name+'}')
    (out/'paper.tex').write_text(flat, encoding='utf-8')
    if not args.no_build:
        run_build(out)
    for rel, expected in hashes.items():
        if sha256(ROOT/rel) != expected:
            raise RuntimeError(f'Source changed during export: {rel}')
    report = {'draft_author_fields_pending': pending, 'compiled': not args.no_build,
              'source_files_sha256': hashes, 'source_to_flat': mapping,
              'source_inputs_unchanged': True}
    report_path = out.parent / (out.name+'_manifest.json')
    report_path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    zp = args.zip_path.resolve(); zp.parent.mkdir(parents=True, exist_ok=True)
    names = {'paper.tex', *mapping.values()}
    if not args.no_build:
        names.add('paper.pdf')
    with zipfile.ZipFile(zp, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(names):
            archive.write(out/name, arcname=name)
    print(f'Flat source archive: {zp}')
    print(f'Export manifest: {report_path}')
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
