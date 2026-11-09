LATEXMK ?= latexmk
CLASS_FILES := cas-sc.cls cas-common.sty
SUBMISSION_TEX := submission/author_details.tex submission/author_declarations.tex
# Frozen figure and table inputs are supplied with the manuscript.
MANUSCRIPT_DATA := $(wildcard data/*.csv data/riemann_lc/*.csv generated/*.csv)

.PHONY: all clean distclean
all: paper.pdf supplementary.pdf

paper.pdf: paper.tex $(CLASS_FILES) $(SUBMISSION_TEX) References.bib $(MANUSCRIPT_DATA)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

supplementary.pdf: supplementary.tex $(MANUSCRIPT_DATA)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error supplementary.tex

clean:
	$(LATEXMK) -c paper.tex
	$(LATEXMK) -c supplementary.tex

# Keep the supplied data so the paper remains buildable after cleaning.
distclean: clean
	rm -f paper.pdf supplementary.pdf
