PYTHON ?= python3

.PHONY: all figures validate test
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

test:
	$(PYTHON) scripts/test_models.py
