# Mesh-to-mesh transfer paper — reviewed Euler and fixed-body RANS evidence

## Current revision

`paper.tex` and `paper.pdf` incorporate the SA/SST RANS exports into the scientifically reviewed manuscript. The study concerns discrete mesh replacement on fixed physical domains, including remeshing-induced transients about nominally steady airfoil flows. No moving-body/ALE campaign is required or represented as completed.

`RANS_ASSESSMENT_AND_INTEGRATION.md` is the current scientific assessment and change report. The previous `SCIENTIFIC_REVIEW_REPORT.md`, journal strategy, and older execution plans are retained as historical project inputs; statements there that RANS data are absent or moving-body work is required are superseded by this revision.

## Build

Use Python 3.10+ with `numpy` and `pandas` (see `requirements.txt`) and a TeX distribution with `latexmk`, pdfLaTeX, Biber, biblatex, PGFPlots/PGFPlotstable, siunitx and the standard packages in the preamble.

```bash
python -m pip install -r requirements.txt
make
```

No shell escape is needed. Run from the project root. `make clean && make` rebuilds LaTeX; `make distclean && make` also recreates generated data.

## Data and updates

- Existing 30 input CSVs remain under `data/`.
- Existing nine generated static CSVs remain under `generated/`.
- All 14 newly supplied CSVs retain their filenames under `RANS_PostProcessing/`.
- The supplied source report and writing prompt are under `MDFiles/`.
- The canonical bibliography is exactly `References.bib`; its 44 previous entries are retained and two model references are appended.

All 53 original CSV inputs/generated static files were preserved byte-for-byte. New derived RANS views live under `generated/rans_*`. `scripts/build_rans_data.py` removes exact duplicate records only from those views, rejects conflicting duplicate keys, aligns physical-event windows, and recomputes the statistics. It never edits the input exports.

After replacing the contents of an existing CSV, run `make`. Dependencies regenerate the affected views and LaTeX plots/tables. Several key text numerals are also generated. **Scientific interpretations still need review after data changes.** The script intentionally fails on an unrecognized schema or changed campaign structure rather than silently discarding data.

```bash
make rans-audit
make audit
python scripts/audit_results.py --check-preservation
```

The final command checks the *older* static/Euler snapshot. New RANS input preservation is recorded in `audit/rans_input_sha256.json`; `audit/rans_data_audit.json` reports changed hashes against this delivered snapshot. An intentional numerical update should produce a new versioned manifest, not overwrite the historical one.

## Evidence and provenance

Riemann and earlier Euler NACA inputs are report-level summary CSVs; their filenames do not imply complete solver histories. The new RANS archive contains actual reduced time series and surface profiles. The source simulation code, original VTU/restart vectors, and original extractor scripts are not included. Reproducing plots from the available CSVs is not equivalent to rerunning the CFD or its raw-field extraction.

The SST floor table is generated from the supplied Markdown report and labeled accordingly. The ordinary OLS standard errors are not treated as autocorrelation-robust significance tests. RANS pressure onsets do not establish an acoustic speed.

`audit/rans_build_validation.json` records compilation, citation/cross-reference checks, input hashes, derived-data reproducibility, and visual inspection. `audit/rans_manuscript_changes.diff` records the LaTeX changes from the reviewed baseline.

Do not run historical `scripts/update_paper.py` on this revision; it belongs to an earlier manuscript layout and is not a build dependency.
