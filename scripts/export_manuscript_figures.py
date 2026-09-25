#!/usr/bin/env python3
"""Export the manuscript's two data-driven figures as separate vector PDFs.

Run from any directory with a current TeX distribution. Figures are extracted from
paper.tex so their plotted values, labels, and styles cannot diverge from the paper.
The scientific CSV files are read in place and are never written by this utility.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]

def strip_command(text: str, command: str) -> tuple[str, str]:
    m = re.search(r'\\'+command+r'\s*\{', text)
    if not m:
        raise ValueError(f'Missing \\{command} in figure')
    start = m.end(); level = 1; end = start
    while level:
        if end >= len(text):
            raise ValueError('Unbalanced figure command')
        c=text[end]
        if c=='{' and (end==0 or text[end-1]!='\\'): level+=1
        if c=='}' and (end==0 or text[end-1]!='\\'): level-=1
        end+=1
    return text[:m.start()]+text[end:], text[start:end-1]

def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=ROOT/'submission'/'figures')
    args=ap.parse_args(); out=args.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
    source=(ROOT/'paper.tex').read_text(encoding='utf-8')
    preamble=source.split(r'\begin{document}',1)[0]
    preamble=re.sub(r'\\documentclass(?:\[[^]]*\])?\{[^}]*\}',r'\\documentclass[12pt]{article}',preamble,count=1)
    preamble=re.sub(r'^\\(?:journal|biboptions)\{[^\n]*\}\s*$', '',preamble,flags=re.M)
    preamble+='\n'+r'\usepackage[active,tightpage]{preview}'+'\n'+r'\PreviewEnvironment{preview}'+'\n'+r'\setlength\PreviewBorder{3pt}'+'\n'
    figs=re.findall(r'\\begin\{figure\}(?:\[[^]]*\])?.*?\\end\{figure\}',source,re.S)
    if len(figs)!=2:
        raise ValueError(f'Expected two manuscript figures; found {len(figs)}. Review export before proceeding.')
    initial={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.rglob('*.csv')}
    captions=[]; manifest=[]
    with tempfile.TemporaryDirectory(prefix='cf-figure-') as tmp:
        build=Path(tmp)
        for n,fig in enumerate(figs,1):
            content,caption=strip_command(fig,'caption')
            content=re.sub(r'\\label\{[^}]*\}', '',content)
            content=re.sub(r'\\begin\{figure\}(?:\[[^]]*\])?', r'\\begin{preview}\\begin{minipage}{\\textwidth}', content)
            content=content.replace(r'\end{figure}',r'\end{minipage}\end{preview}')
            tex=build/f'Figure_{n}.tex'
            tex.write_text(preamble+'\n\\begin{document}\n'+content+'\n\\end{document}\n')
            for _ in range(2):
                result=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',f'-output-directory={build}',str(tex)],cwd=ROOT,text=True,capture_output=True)
                if result.returncode:
                    raise RuntimeError(result.stdout[-6000:]+'\n'+result.stderr)
            log=(build/f'Figure_{n}.log').read_text(errors='replace')
            if re.search(r'Overfull \\[hv]box|undefined references|! LaTeX Error',log):
                raise RuntimeError(f'Figure {n} has compilation warnings requiring review.')
            dst=out/f'Figure_{n}.pdf'; shutil.copyfile(build/f'Figure_{n}.pdf',dst)
            captions.append((n,caption))
            manifest.append({'figure':n,'file':dst.name,'sha256':hashlib.sha256(dst.read_bytes()).hexdigest()})
    final={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.rglob('*.csv')}
    if initial!=final: raise RuntimeError('Source data changed during figure export.')
    caption_text='\n\n'.join(f'Figure {n}. {text}' for n,text in captions)+'\n'
    (out/'figure_captions.tex').write_text(caption_text,encoding='utf-8')
    plain=caption_text.replace(r'\Delta C_L', 'ΔC_L').replace(r'L_2','L₂').replace('$','').replace('--','–')
    (out/'figure_captions.txt').write_text(plain,encoding='utf-8')
    (out/'figure_export_manifest.json').write_text(json.dumps({'source_tex_sha256':hashlib.sha256(source.encode()).hexdigest(),'csv_inputs_unchanged':True,'figures':manifest},indent=2)+'\n')
    print(f'Exported {len(figs)} vector PDFs and their exact manuscript captions to {out}')
if __name__=='__main__': main()
