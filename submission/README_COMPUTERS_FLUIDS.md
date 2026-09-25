# Computers & Fluids submission support

This directory accompanies the current Computers & Fluids manuscript.

The manuscript source reads author information from `author_details.tex` and declarations from `author_declarations.tex`. Complete and confirm those fields before final submission.

The repository-level `Makefile` provides manuscript and submission checks. The helper scripts under `scripts/` can perform the mechanical manuscript preflight, export the two vector figures, and prepare a flat Elsevier source package without renaming or modifying the scientific CSV inputs.

The manuscript-facing nomenclature is intentionally independent of raw dataset labels. Raw data filenames and column labels remain unchanged for reproducibility; scientific names are applied only in the paper and exported presentation layer.

The separate figure PDFs and captions correspond to the two data-driven figures in the current manuscript. The final manuscript itself remains the authoritative source for captions and interpretation.

Before submission, confirm author list, affiliations, funding, CRediT roles, competing interests, and the AI-assistance declarations.
