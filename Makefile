.PHONY: help clean clean-pyc clean-build list test test-all coverage docs release sdist

help:
	@echo "clean-build - remove build artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "lint - check style with flake8"
	@echo "test - run tests quickly with the default Python"
	@echo "coverage - check code coverage quickly with the default Python"
	@echo "docs - generate Sphinx HTML documentation, including API docs"

clean: clean-build clean-pyc

clean-build:
	rm -fr build/
	rm -fr dist/
	rm -fr *.egg-info

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +

docs:
	rm -f docs/seismic.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ seismic legacy seismic/**/sandbox seismic/ml_classifier seismic/inventory/legacy \
	    seismic/**/example_*.py seismic/ASDFdatabase/ASDF_build*.py seismic/ASDFdatabase/minimus*.py
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	@echo "Run 'xdg-open docs/_build/html/index.html' to view documentation"

lint:
	pytest --junit-xml=test_output/flake8/results.xml \
	    --flake8 -p no:regtest	--cache-clear seismic

test: mpi
	export ELLIPCORR=$PWD/ellip-corr
	pytest --junit-xml=test_output/pytest/results.xml --cache-clear ./tests

coverage: mpi
	pytest --junit-xml=test_output/pytest/results.xml --cov \
	    --cov-report=html:test_output/coverage \
	    --cov-report=term-missing:skip-covered \
	    --cov-fail-under=60 --cache-clear ./tests

mpi:
	mpirun --allow-run-as-root -n 2 pytest tests/test_pyasdf.py
