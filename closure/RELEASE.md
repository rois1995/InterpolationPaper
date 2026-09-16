# Release v1.0

Commit 6d59bf6, tag v1.0, https://github.com/rois1995/InterpolationPaper.git.
Archive: InterpolationPaper_v1.0.tar.gz (git archive of the tag), SHA-256 9aa4cfd4962a76780cd10a97c5222038384b4f4d36c1471c79485cc9efc1ee29.

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
