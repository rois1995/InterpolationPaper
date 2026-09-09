PYTHON ?= python3
LATEXMK ?= latexmk
MAIN := paper.tex
PDF := paper.pdf
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
	data/naca_surface_delta_cp.csv \
	data/naca_gp3_verification.csv
GENERATED_STAMP := generated/.stamp
RANS_DATA := $(wildcard RANS_PostProcessing/*.csv)
RANS_STAMP := generated/.rans_stamp

.PHONY: all data audit rans-audit clean distclean

all: $(PDF)

data: $(GENERATED_STAMP) $(RANS_STAMP)

$(GENERATED_STAMP): scripts/build_plot_data.py $(DATA)
	$(PYTHON) scripts/build_plot_data.py --data-dir data --out-dir generated
	@touch $(GENERATED_STAMP)

$(RANS_STAMP): scripts/build_rans_data.py $(RANS_DATA) MDFiles/NACA0012_STUDY_REPORT.md
	$(PYTHON) scripts/build_rans_data.py
	@touch $(RANS_STAMP)

$(PDF): $(MAIN) References.bib $(GENERATED_STAMP) $(RANS_STAMP)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN)

clean:
	$(LATEXMK) -c $(MAIN)

# Also removes generated tables; the raw CSV files are never deleted.
distclean: clean
	rm -f generated/*.csv generated/*.txt generated/.stamp generated/.rans_stamp generated/rans_numbers.tex $(PDF)

# Recompute manuscript-table checks without altering the numerical data.
audit:
	$(PYTHON) scripts/audit_results.py

rans-audit:
	$(PYTHON) scripts/build_rans_data.py
