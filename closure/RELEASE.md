# Release v1.0

Commit (see git tag v1.0), tag v1.0, https://github.com/rois1995/InterpolationPaper.git.
Archive: InterpolationPaper_v1.0.tar.gz (git archive of the tag); its SHA-256 is in the sidecar file InterpolationPaper_v1.0.tar.gz.sha256 next to it.

Release checks (closure plan T16), all run on 2026-09-16:

1. Clean build in a fresh clone (make distclean; make): PDF builds, 0 LaTeX errors, 0 undefined references.
2. Regeneration: make data recreates every generated/ file byte-identical to the committed one.
3. Bad inputs: an exact duplicate row is removed and reported; a conflicting duplicate key raises
   ValueError; an unknown method/order raises "Unsupported method/order combination".
4. Perturbation: changing one value of data/naca_superposition.csv in a disposable copy changes the
   corresponding table entry in the rebuilt PDF; the release data were not altered.
5. Every CSV name and location referenced by paper.tex exists; new files were added to the Makefile.
6. Every page and both figures inspected on rendered images: axes, guides, legends, units, table layout.
7. Development narrative removed from the text; no stale settings, markers or superseded values remain.
8. REPRODUCE.md gives the commands for the reduced-data level and the solver level.
9. Code revisions (provenance/code_revision_*.txt) and mesh/restart checkpoints
   (provenance/checkpoints_sha256.txt) are recorded; this tag is the immutable identifier.
