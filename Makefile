PYTHON ?= python3
LATEXMK ?= latexmk
MAIN := paper.tex
PDF := paper.pdf
SUBMISSION_TEX := submission/author_details.tex submission/author_declarations.tex
# Elsevier CAS single-column class files (template v2.4), kept unmodified beside paper.tex.
CLASS_FILES := cas-sc.cls cas-common.sty
DATA := \
	data/Comparison_Errors.csv \
	data/Comparison_Conservation.csv \
	data/Comparison_Orders.csv \
	data/discontinuous_convergence_L2.csv \
	data/discontinuous_conservation_defect.csv \
	data/discontinuous_boundedness.csv \
	data/riemann_refinement_orders.csv \
	data/riemann_refinement_errors.csv \
	data/riemann_global_correction.csv \
	data/riemann_multicycle.csv \
	data/riemann_transition_step399.csv \
	data/riemann_stress_representative.csv \
	data/riemann_frequency_full.csv \
	data/riemann_amplitude_global.csv \
	data/riemann_moment_preservation.csv \
	data/riemann_limiter_activation.csv \
	data/naca_native_convergence.csv \
	data/naca_force_drift.csv \
	data/naca_acoustic_delay.csv \
	data/naca_acoustic_delay_fits.csv \
	data/naca_pressure_propagation.csv \
	data/naca_transfer_defects.csv \
	data/naca_initial_residual_identity.csv \
	data/naca_first_step_response.csv \
	data/naca_first_step_CL_plot.csv \
	data/naca_superposition.csv \
	data/naca_global_correction.csv \
	data/naca_force_history.csv \
	data/naca_identity_wall_split.csv \
	data/naca_initial_residual_identity_pairs.csv \
	data/naca_initial_residual_identity_double.csv \
	data/riemann_amplitude_seeds.csv \
	data/naca_cancellation.csv \
	data/naca_amplitude_sequence.csv
GENERATED_STAMP := generated/.stamp
RANS_DATA := $(wildcard RANS_PostProcessing/*.csv)
RANS_STAMP := generated/.rans_stamp
LC_DATA := $(wildcard data/riemann_lc/*)
LC_STAMP := generated/.lc_stamp

.PHONY: all data audit rans-audit clean distclean cf-author-check cf-package cf-draft-package

all: $(PDF) supplementary.pdf

data: $(GENERATED_STAMP) $(RANS_STAMP) $(LC_STAMP)

$(GENERATED_STAMP): scripts/build_plot_data.py $(DATA)
	$(PYTHON) scripts/build_plot_data.py --data-dir data --out-dir generated
	@touch $(GENERATED_STAMP)

$(RANS_STAMP): scripts/build_rans_data.py scripts/build_views.py $(RANS_DATA) MDFiles/NACA0012_STUDY_REPORT.md
	$(PYTHON) scripts/build_rans_data.py
	$(PYTHON) scripts/build_views.py
	@touch $(RANS_STAMP)

$(LC_STAMP): scripts/build_riemann_lc_data.py provenance/riemann_lc_inputs.json $(LC_DATA)
	$(PYTHON) scripts/build_riemann_lc_data.py
	@touch $(LC_STAMP)

$(PDF): $(MAIN) $(CLASS_FILES) $(SUBMISSION_TEX) References.bib $(GENERATED_STAMP) $(RANS_STAMP) $(LC_STAMP)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN)

supplementary.pdf: supplementary.tex $(GENERATED_STAMP) $(RANS_STAMP) $(LC_STAMP)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error supplementary.tex

clean:
	$(LATEXMK) -c $(MAIN)
	$(LATEXMK) -c supplementary.tex

# Also removes generated tables; the raw CSV files are never deleted.
distclean: clean
	rm -f generated/*.csv generated/*.txt generated/.stamp generated/.rans_stamp generated/.lc_stamp generated/rans_numbers.tex $(PDF) supplementary.pdf

# Recompute manuscript-table checks without altering the numerical data.
audit: $(LC_STAMP)
	$(PYTHON) scripts/audit_results.py
	$(PYTHON) scripts/check_typed_tables.py
	$(PYTHON) -m unittest discover -s tests -p test_riemann_lc_data.py

rans-audit:
	$(PYTHON) scripts/build_rans_data.py

# Export-only targets: source datasets and their locations are not changed.
cf-author-check:
	$(PYTHON) scripts/prepare_cf_submission.py --check-only

cf-package: all cf-figures
	$(PYTHON) scripts/prepare_cf_submission.py

cf-draft-package: all cf-figures
	$(PYTHON) scripts/prepare_cf_submission.py --allow-incomplete

# Mechanical preflight reports author fields separately from manuscript checks.
.PHONY: cf-preflight cf-figures
cf-preflight:
	$(PYTHON) scripts/check_manuscript.py

cf-figures:
	$(PYTHON) scripts/export_manuscript_figures.py
