PYTHON ?= python3

.PHONY: all figures validate test junction junction-validate junction-figures
all: validate

figures:
	$(PYTHON) scripts/generate_figures.py
	$(PYTHON) scripts/solve_contact_transparency.py
	$(PYTHON) scripts/reciprocal_device.py
	$(PYTHON) scripts/robustness.py
	$(PYTHON) scripts/spectral_tests.py
	$(PYTHON) scripts/hahn_transport.py
	$(PYTHON) scripts/quantum_helix.py

validate: figures
	$(PYTHON) scripts/validate_outputs.py
	$(PYTHON) scripts/test_models.py
	$(PYTHON) scripts/test_junction_transport.py
	$(PYTHON) scripts/validate_junction.py

test:
	$(PYTHON) scripts/test_models.py
	$(PYTHON) scripts/test_junction_transport.py

junction:
	$(PYTHON) scripts/benchmark_junction.py
	$(PYTHON) scripts/validate_junction.py

junction-validate:
	$(PYTHON) scripts/test_junction_transport.py
	$(PYTHON) scripts/validate_junction.py

junction-figures:
	$(PYTHON) scripts/benchmark_junction.py --plots-only
