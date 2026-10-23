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
    # Acronyms and a word stem of their expansion. The wording is free: the first use of
    # each acronym in a scope must itself be the definition "full name (ACR)", with an
    # optional order suffix such as NN1 or CV2, and the stem within the preceding words.
    expansions={
      'NN':'nearest-neighbor','BAR':'barycentric','GP':'Galerkin projection',
      'CV':'control-volume transfer','GC':'global correction',
      'RANS':'Reynolds-averaged Navier--Stokes',
      'BDF2':'backward differentiation formula','BJ':'Barth--Jespersen',
      'SF':'safety factor',
      'NACA':'National Advisory Committee for Aeronautics','SA':'Spalart--Allmaras',
      'SST':'shear-stress-transport','SU2':'Stanford University Unstructured',
      'HLLC':'Harten--Lax--van Leer contact',
      'MUSCL':'monotonic upstream-centered scheme for conservation laws',
      'CFL':'Courant--Friedrichs--Lewy','FGMRES':'flexible generalized minimal residual',
      'ILU':'incomplete lower--upper','BiCGSTAB':'biconjugate-gradient',
      'RMS':'root-mean-square','AI':'artificial intelligence',
      'CRediT':'Contributor Roles Taxonomy',
    }
    suffixed={'NN','BAR','GP','CV'}  # method families carry an order suffix
    # The abstract must stand alone, so it is checked as its own scope; the main text
    # (after \maketitle) is checked again from its own start.
    scopes={'abstract':abstract,'main text':s.split(r'\maketitle',1)[1]}
    def norm(x: str) -> str:
        # Compare expansions without case, line breaks or the en-dash spelling.
        return re.sub(r'\s+',' ',x.replace('--','-')).lower()
    acronym_checks={}
    for scope,text in scopes.items():
        for key,stem in expansions.items():
            pat=r'(?<![A-Za-z0-9\\])'+re.escape(key)+(r'[1-3]?' if key in suffixed else '')+r'(?![A-Za-z0-9])'
            first=re.search(pat,text)
            if not first:
                acronym_checks[f'{scope}:{key}']='unused'; continue
            before=text[max(0,first.start()-160):first.start()]
            after=text[first.end():first.end()+1]
            # Accepted forms: "full name (ACR)" and "full name, denoted ACR".
            in_parens=before.endswith('(') and after==')'
            denoted=bool(re.search(r'denoted\s*$',before))
            ok=(in_parens or denoted) and norm(stem) in norm(before)
            acronym_checks[f'{scope}:{key}']='defined at first use' if ok else 'used before definition'
            if not ok:
                ln=text[:first.start()].count('\n')+1
                errors.append(f'Review first-use definition: {key} ({scope}, line {ln} of scope: '
                              +re.sub(r'\s+',' ',text[max(0,first.start()-50):first.end()+10])+')')
    kw=re.search(r'\\begin\{keywords\}(.*?)\\end\{keywords\}',s,re.S).group(1)
    keywords=[x.strip() for x in kw.split(r'\sep')]
    if len(keywords) not in range(1,8): errors.append('Keyword count outside journal 1–7 requirement')
    if abstract_words>250: errors.append('Abstract exceeds 250 words')
    hl=[x[2:].strip() for x in (ROOT/'submission/highlights.txt').read_text().splitlines() if x.startswith('- ')]
    if len(hl) not in range(3,6): errors.append('Highlights count outside publisher 3–5 guidance')
    if any(len(x)>85 for x in hl): errors.append('Highlight exceeds 85 characters')
    # The highlights page of the CAS manuscript must repeat the separate highlights file.
    hl_env=re.search(r'\\begin\{highlights\}(.*?)\\end\{highlights\}',s,re.S)
    hl_tex=[x.strip() for x in hl_env.group(1).split(r'\item')[1:]] if hl_env else []
    if hl_tex!=hl: errors.append('Highlights in paper.tex differ from submission/highlights.txt')
    ai=s.index(r'\section*{Declaration of generative AI')
    if ai<s.index(r'\appendix'): errors.append('Writing assistance declaration precedes appendices')
    if r'\section' in s[ai+len(r'\section*'):s.index(r'\bibliographystyle')]: errors.append('Another section follows writing assistance declaration')
    author_files=['submission/author_details.tex','submission/author_declarations.tex']
    pending=[x for x in author_files if re.search(r'to be supplied|\[Authors:',uncomment((ROOT/x).read_text()))]
    # Do not mistake the absence of author facts for a completed submission.
    if args.submission and pending: errors.append('Author confirmation fields remain: '+', '.join(pending))
    supplement=uncomment((ROOT/'supplementary.tex').read_text())
    supplement_labels=re.findall(r'\\label\{([^}]+)\}',supplement)
    supplement_refs=set(re.findall(r'\\(?:eqref|ref)\*?\{([^}]+)\}',supplement))
    errors.extend('Undefined supplementary reference: '+x for x in supplement_refs-set(supplement_labels))
    errors.extend('Duplicate supplementary label: '+x for x,n in Counter(supplement_labels).items() if n>1)
    for document in ('paper','supplementary'):
        log_path=ROOT/(document+'.log')
        if not log_path.exists():
            errors.append('Build required: '+document+'.tex')
            continue
        log=log_path.read_text(errors='replace')
        layout_warnings=re.findall(r'Overfull \\[hv]box[^\n]*|undefined references[^\n]*|undefined citations[^\n]*|! LaTeX Error[^\n]*|destination with the same identifier[^\n]*',log)
        errors.extend(document+': '+warning for warning in layout_warnings)
    csvs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.rglob('*.csv')) if 'dist' not in p.parts}
    blg=(ROOT/'paper.blg').read_text(errors='replace') if (ROOT/'paper.blg').exists() else ''
    report={'mechanical_errors':errors,'author_confirmations_pending':pending,
      'journal_specific_guide_status':'Journal of Computational and Applied Mathematics Guide for Authors (ScienceDirect) checked 2026-10-23: abstract <=250 words, 1-7 keywords, 3-5 highlights of <=85 characters, numbered references; confirm portal requirements at submission.',
      'abstract_words':abstract_words,'keyword_count':len(keywords),'highlight_characters':[len(x) for x in hl],
      'defined_acronyms_at_first_source_use':acronym_checks,
      'source_labels':len(labels),'bibliography_entries':len(keys),'cited_keys':len(citekeys),
      'supplementary_source_labels':len(supplement_labels),
      'bibtex_notices':[x for x in blg.splitlines() if x.startswith('Warning--')],
      'csv_sha256':csvs,
      'scope':'Mechanical source checks; does not establish scientific validity, detect every undefined symbol, or verify every cited paper.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='csv_sha256'},indent=2,ensure_ascii=False))
    return 1 if errors else 0
if __name__=='__main__': sys.exit(main())
