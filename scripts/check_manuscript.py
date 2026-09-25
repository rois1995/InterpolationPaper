#!/usr/bin/env python3
"""Mechanical preflight for this manuscript; not a substitute for scientific review.

Checks cross-references, citation keys, acronym introductions, supplied metadata,
and the latest compilation log. Author confirmations are reported separately;
--submission makes unresolved author fields fail the check.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys
ROOT=Path(__file__).resolve().parents[1]

def uncomment(s: str) -> str:
    return re.sub(r'(?<!\\)%[^\n]*','',s)

def macro_body(s: str, name: str) -> str:
    m=re.search(r'\\newcommand\{\\'+re.escape(name)+r'\}\{',s)
    if not m: raise ValueError(f'Missing author macro: {name}')
    start=m.end(); end=start; depth=1
    while depth:
        if end>=len(s): raise ValueError('Unbalanced author declaration')
        c=s[end]
        if c=='{' and s[end-1]!='\\': depth+=1
        if c=='}' and s[end-1]!='\\': depth-=1
        end+=1
    return s[start:end-1]

def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--submission',action='store_true')
    ap.add_argument('--output',type=Path,default=ROOT/'submission'/'MANUSCRIPT_PREFLIGHT.json')
    args=ap.parse_args()
    full=uncomment((ROOT/'paper.tex').read_text())
    s=full.split(r'\begin{document}',1)[1]
    author=uncomment((ROOT/'submission/author_declarations.tex').read_text())
    for name in ['FundingStatement','AuthorContributionStatement','CompetingInterestStatement','ResearchAssistanceStatement','ManuscriptAssistanceStatement']:
        s=s.replace('\\'+name,macro_body(author,name))
    bib=(ROOT/'References.bib').read_text()
    keys=re.findall(r'@\w+\s*\{\s*([^,\s]+)',bib)
    citekeys=set()
    for m in re.finditer(r'\\cite\w*\*?(?:\[[^]]*\])*\{([^}]+)\}',s): citekeys.update(m.group(1).split(','))
    labels=re.findall(r'\\label\{([^}]+)\}',s)
    # named pgfplots legend is written to the auxiliary file rather than \label.
    named=re.findall(r'legend to name=([^,\s]+)',s)
    refs=set(re.findall(r'\\(?:eqref|ref|pageref|autoref)\*?\{([^}]+)\}',s))
    errors=[]
    for x in refs-set(labels)-set(named): errors.append('Undefined source reference: '+x)
    for x in citekeys-set(keys): errors.append('Missing bibliography key: '+x)
    duplicates=[x for x,n in Counter(labels).items() if n>1]
    errors.extend('Duplicate label: '+x for x in duplicates)
    errors.extend('Duplicate bibliography key: '+x for x,n in Counter(keys).items() if n>1)
    abstract=re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',s,re.S).group(1)
    abstract_words=len(abstract.split())
    if re.search(r'\\cite',abstract): errors.append('Review citation in abstract')
    definitions={
      'NN':r'NN for nearest-neighbor Taylor interpolation',
      'BAR':r'BAR for barycentric interpolation',
      'GP':r'GP for Galerkin projection',
      'RANS':r'Reynolds-averaged Navier--Stokes \(RANS\)',
      'BDF2':r'second-order backward differentiation formula \(BDF2\)',
      'BJ':r'Barth--Jespersen \(BJ\)',
      'GC':r'global correction, denoted GC',
      'SF':r'safety factor \(SF\)',
      'NACA':r'National Advisory Committee for Aeronautics \(NACA\)',
      'SA':r'Spalart--Allmaras \(SA\)',
      'SST':r'shear-stress-transport \(SST\)',
      'SU2':r'Stanford University Unstructured \(SU2\)',
      'HLLC':r'Harten--Lax--van Leer contact \(HLLC\)',
      'MUSCL':r'monotonic upstream-centered scheme for conservation laws \(MUSCL\)',
      'CFL':r'Courant--Friedrichs--Lewy \(CFL\)',
      'FGMRES':r'flexible generalized minimal residual \(FGMRES\)',
      'ILU':r'incomplete lower--upper \(ILU\)',
      'RMS':r'root-mean-square \(RMS\)',
      'AI':r'artificial intelligence \(AI\)',
      'CRediT':r'Contributor Roles Taxonomy \(CRediT\)',
    }
    acronym_checks={}
    for key,definition in definitions.items():
        first=re.search(r'(?<![A-Za-z0-9])'+key+r'(?:[1234])?(?![A-Za-z0-9])',s)
        defined=re.search(definition,s)
        ok=bool(first and defined and defined.start()<=first.start()<defined.end())
        acronym_checks[key]=ok
        if not ok: errors.append('Review first-use definition: '+key)
    kw=re.search(r'\\begin\{keyword\}(.*?)\\end\{keyword\}',s,re.S).group(1)
    keywords=[x.strip() for x in kw.split(r'\sep')]
    hl=[x[2:].strip() for x in (ROOT/'submission/highlights.txt').read_text().splitlines() if x.startswith('- ')]
    if len(hl) not in range(3,6): errors.append('Highlights count outside publisher 3–5 guidance')
    if any(len(x)>85 for x in hl): errors.append('Highlight exceeds 85 characters')
    ai=s.index(r'\section*{Declaration of generative AI')
    if ai<s.index(r'\appendix'): errors.append('Writing assistance declaration precedes appendices')
    if r'\section' in s[ai+len(r'\section*'):s.index(r'\bibliographystyle')]: errors.append('Another section follows writing assistance declaration')
    author_files=['submission/author_details.tex','submission/author_declarations.tex']
    pending=[x for x in author_files if re.search(r'to be supplied|\[Authors:',uncomment((ROOT/x).read_text()))]
    # Do not mistake the absence of author facts for a completed submission.
    if args.submission and pending: errors.append('Author confirmation fields remain: '+', '.join(pending))
    log=(ROOT/'paper.log').read_text(errors='replace') if (ROOT/'paper.log').exists() else ''
    layout_warnings=re.findall(r'Overfull \\[hv]box[^\n]*|undefined references[^\n]*|undefined citations[^\n]*|! LaTeX Error[^\n]*|destination with the same identifier[^\n]*',log)
    errors.extend(layout_warnings)
    csvs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*.csv')) if 'dist' not in p.parts}
    blg=(ROOT/'paper.blg').read_text(errors='replace') if (ROOT/'paper.blg').exists() else ''
    report={'mechanical_errors':errors,'author_confirmations_pending':pending,
      'journal_specific_guide_status':'Not verified: live Computers & Fluids guide returned HTTP 403; confirm portal requirements.',
      'abstract_words':abstract_words,'keyword_count':len(keywords),'highlight_characters':[len(x) for x in hl],
      'defined_acronyms_at_first_source_use':acronym_checks,
      'source_labels':len(labels),'bibliography_entries':len(keys),'cited_keys':len(citekeys),
      'bibtex_notices':[x for x in blg.splitlines() if x.startswith('Warning--')],
      'csv_sha256':csvs,
      'scope':'Mechanical source checks; does not establish scientific validity, detect every undefined symbol, or verify every cited paper.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='csv_sha256'},indent=2,ensure_ascii=False))
    return 1 if errors else 0
if __name__=='__main__': sys.exit(main())
