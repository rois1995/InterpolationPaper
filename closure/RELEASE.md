# Release v1.1

Tag `v1.1` on https://github.com/rois1995/InterpolationPaper.git.
Archive: InterpolationPaper_v1.1.tar.gz (git archive of the tag), with its SHA-256 in the
sidecar file InterpolationPaper_v1.1.tar.gz.sha256 next to it.

Contents beyond v1.0: the numerical controls of the RANS transfer comparison (dual-time
iteration caps of 1600 and 6400; physical steps of dt/2 for both closures and dt/4 for SA,
with the archived displacement fields replayed so the mesh events keep their physical times;
the same remeshing window on a finer viscous mesh; no-transfer baselines on three viscous
meshes), the replicated Riemann global-correction counterexample on three deformation
realizations, and the BDF2 history-cancellation and perturbation-amplitude checks.

Release checks (closure plan T16), all run on 2026-09-17 in disposable clones:

1. Clean build in a fresh clone (make distclean; make): PDF builds, 0 LaTeX errors,
   0 undefined references, 25 pages.
2. Regeneration: make data recreates every generated/ file byte-identical to the committed one.
3. Preservation: scripts/audit_results.py --check-preservation hashes all 103 released tables
   (data/ with subdirectories, generated/, RANS_PostProcessing/) against
   audit/original_csv_sha256.json; no file changed or missing.
4. Bad inputs: an exact duplicate row is removed and reported; a conflicting duplicate key
   raises ValueError; an unknown method/order is rejected.
5. Perturbation: changing one value of data/naca_superposition.csv in a disposable copy changes
   the corresponding table entry in the rebuilt PDF; the release data were not altered.
6. Typed tables: scripts/check_typed_tables.py verifies the two tables written literally in
   paper.tex against their CSV sources; it runs as part of `make audit`.
7. Every CSV name and location referenced by paper.tex exists and is declared in the Makefile.
8. REPRODUCE.md gives the commands for the reduced-data level and the solver level.
9. Code revisions (provenance/code_revision_*.txt) and mesh/restart checkpoints
   (provenance/checkpoints_sha256.txt, 28 entries) are recorded; this tag is the immutable
   identifier.

Known limits, stated in the manuscript: the three-mesh sequence is monotone for lift, drag and
pressure drag but not in its asymptotic range; the paired NN2-GP2 load difference is not
time-step converged; the SST inner loop has an omega-residual floor on the moved grids.
