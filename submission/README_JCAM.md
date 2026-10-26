# Journal of Computational and Applied Mathematics submission support

This directory accompanies the manuscript prepared for the Journal of Computational and Applied Mathematics (Elsevier), typeset with the Elsevier CAS single-column class (`cas-sc.cls`, `cas-common.sty`, template v2.4) and numbered references (`elsarticle-num.bst`).

The manuscript source reads author information from `author_details.tex` (CAS front-matter syntax) and declarations from `author_declarations.tex`. Complete and confirm those fields before final submission.

The repository-level `Makefile` provides manuscript and submission checks. The helper scripts under `scripts/` can perform the mechanical manuscript preflight, export the six vector figures, and prepare a flat Elsevier source package (including the CAS class files) without renaming or modifying the scientific CSV inputs.

The manuscript-facing nomenclature is intentionally independent of raw dataset labels. Raw data filenames and column labels remain unchanged for reproducibility; scientific names are applied only in the paper and exported presentation layer.

The submission export includes `supplementary.tex`, its CSV dependencies, and a separately compiled `supplementary.pdf`. Upload the latter as supplementary material alongside the main article.

The separate figure PDFs and captions correspond to the six data-driven figures in the manuscript. The final manuscript itself remains the authoritative source for captions and interpretation.

`highlights.txt` is the separate editable highlights file (3 to 5 bullets, at most 85 characters each); the highlights page of the manuscript repeats it.

Before submission, confirm author list, affiliations, the corresponding author's e-mail address, funding, CRediT roles, competing interests (declarations tool), the data statement, suggested reviewers, and the AI-assistance declarations.
