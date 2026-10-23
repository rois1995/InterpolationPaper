# Mesh-to-mesh transfer paper: static, Euler, BDF2-history and fixed-body RANS evidence

`paper.tex` and `paper.pdf` are the manuscript. It compares four transfer families on fixed physical
domains: nearest-neighbor Taylor interpolation (NN1-NN3), barycentric interpolation (BAR2, BAR3),
continuous Galerkin projection (GP2, GP3) and the conservative control-volume transfer CV2 (median-dual
control volumes, P1 hat-function field, MeshAdaptation method `Conservative`). The evidence covers static
smooth and discontinuous fields, repeated Euler remeshing of a Riemann problem, controlled BDF2 restarts
and remeshing transients on the NACA 0012 airfoil, and fixed-body SA/SST RANS calculations with their
numerical-sensitivity controls. No moving-body/ALE calculation is part of the study.

All transfer calculations ran with MeshAdaptation `Scripts/` at commit b3054414 and SU2 v8 commit
dcee051d (see `REPRODUCE.md` and `provenance/`).

## Build

Use Python 3.8+ with `numpy` and `pandas` (see `requirements.txt`) and a TeX distribution with `latexmk`,
pdfLaTeX, BibTeX (`elsarticle-num.bst`), the STIX fonts, PGFPlots/PGFPlotstable, siunitx and the standard packages in the
preamble. The Elsevier CAS single-column class files `cas-sc.cls` and `cas-common.sty` (template v2.4) are included in the
project root.

```bash
python -m pip install -r requirements.txt
make
```

No shell escape is needed. Run from the project root. `make clean && make` rebuilds LaTeX;
`make distclean && make` also recreates the generated data.

## Contents

- `data/`: 34 summary CSVs of the static, Euler, BDF2 and Riemann analyses; `data/bdf2_raw/` holds the
  six per-pair tables of the BDF2 study from which the BDF2 summaries are formed.
- `RANS_PostProcessing/`: the 28 CSV exports of the RANS analysis (time series, surface profiles, transfer
  integrals, controls, provenance), its data-section report and its log.
- `generated/`: views derived from the two folders above by `scripts/build_plot_data.py`,
  `scripts/build_rans_data.py` and `scripts/build_views.py`; `generated/rans_numbers.tex` holds text
  numerals computed from the RANS exports.
- `MDFiles/`: the NACA 0012 study report (also the source of the SST floor table) and the writing prompt.
- `audit/`, `closure/`: digests of the released CSV set and the audit reports; `provenance/`: code
  revisions of the case bases and checkpoint digests.
- `submission/`: author details, declarations and the journal package files.
- `References.bib`: the bibliography (49 entries).

`scripts/build_rans_data.py` removes exact duplicate records only from the derived views, rejects
conflicting duplicate keys, aligns physical-event windows, and recomputes the statistics. It never edits
the input exports, and it fails on an unrecognized schema rather than silently discarding data.

After replacing the contents of an input CSV, run `make`: the dependencies regenerate the affected views
and LaTeX plots/tables. Most text numerals are typed in `paper.tex`; `scripts/check_typed_tables.py`
checks the two typed tables against their CSVs, and the scientific interpretation must be reviewed
whenever the data change.

```bash
make rans-audit
make audit
python3 scripts/audit_results.py --check-preservation
```

The last command compares every CSV with `audit/original_csv_sha256.json`, the digest of the released CSV
set; `closure/csv_sha256_baseline.txt` lists the same digests and `audit/rans_input_sha256.json` those of
the inputs read by the RANS builder.

## Evidence and provenance

The Riemann, static and Euler NACA inputs are summary CSVs produced by the analysis scripts of the case
trees; their filenames do not imply complete solver histories. The RANS exports contain reduced time series
and surface profiles. The case trees, VTU/restart fields and simulation code are not included; reproducing
plots from the CSVs is not equivalent to rerunning the CFD or its raw-field extraction (`REPRODUCE.md`).

The ordinary OLS standard errors are not treated as autocorrelation-robust significance tests. RANS pressure
onsets do not establish an acoustic speed.

`scripts/update_paper.py` belongs to a different manuscript layout and is not a build dependency.
