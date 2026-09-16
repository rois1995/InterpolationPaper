# Reproducing the manuscript

Two levels are available. The first needs only this repository; the second needs the solver
installations, the meshes and the case trees named below.

## Level 1: figures, tables and PDF from the reduced data (this repository)

Dependencies: Python 3.8+ with numpy and pandas (`scripts/build_plot_data.py`,
`scripts/build_rans_data.py`); a TeX distribution with `latexmk`, `pgfplots` >= 1.18 (the manuscript
sets `compat=1.18`), `pgfplotstable`, `siunitx`, `booktabs`, `biblatex`.

    make distclean      # removes generated/*.csv, generated/rans_numbers.tex and the PDF
    make data           # rebuilds generated/ from data/ and RANS_PostProcessing/ (byte-identical to the committed files)
    make                # builds paper.pdf
    make audit          # scripts/audit_results.py: table checks, citation and label consistency
    make rans-audit     # scripts/build_rans_data.py: duplicate policy and input digests (audit/rans_data_audit.json)
    python3 scripts/audit_results.py --check-preservation   # every CSV matches audit/original_csv_sha256.json

`data/` and `RANS_PostProcessing/` hold the analysis outputs of the executed runs and are never written
by the build. `generated/` holds derived views only. `closure/csv_sha256_baseline.txt` and
`audit/original_csv_sha256.json` are the digests of the released CSV set.

## Level 2: the simulations behind the reduced data

Code: MeshAdaptation `Scripts/` at commit 35d865eb (transfer and mesh-motion code; the built case
bases are recorded in `provenance/code_revision_*.txt` with per-file SHA-256), SU2 v8 commit
dcee051d (binary of 2025-12-10). Meshes and steady restart checkpoints: `provenance/checkpoints_sha256.txt`.
The runs need 4 MPI ranks and the case trees under `TestInterp/` of the validation workspace:

| Part | Driver | Analysis / export |
|---|---|---|
| static transfer tests | `TestInterp/RunCases_static.sh`, `DiscontinuousBoundedness.sh` | `Comparison.py`, `DiscontinuousTables.py` |
| Euler NACA 0012 | `SteadyNACA/RunAll_Euler.sh` | `ForceDriftAnalysis.py`, `AcousticDelayAnalysis.py`, `PressureProbeAnalysis.py` |
| BDF2 history study | `SteadyNACA/BDF2/run_bdf2_all.sh` | `bdf2_analysis.py`, `bdf2_identity_double.py`, `bdf2_identity_localize.py`, `PatchTests.py` |
| 2D Riemann | `2D_Riemann/RunAll_Riemann.sh`, `RunSeeds_Riemann.sh`, `RunTier0Refine.sh` | the `2D_Riemann/*.py` analyses, `RiemannTransitionExtract.py`, `RiemannTransitionSummary.py` |
| RANS NACA 0012 | `SteadyNACA/run_rans_campaign.sh`, `run_controls_rans.sh` | `run_rans_postprocessing.sh` |
| all of the above in order | `TestInterp/run_everything.sh` | `TestInterp/ExportPaperData.py` writes `data/`, `RANS_PostProcessing/`, commits and pushes |

Limits: the case trees (about 150 GB with the per-step volume fields) and the Pointwise mesh
projects are not part of this repository; the runs took about eight days on 4 ranks.
