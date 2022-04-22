.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	coverage run -m pytest
	coverage xml

install: all
	bash install-pplacer.sh
	pip install git+https://github.com/Ecogenomics/CheckM.git@d74bb68d48b2318542eb7137343196d8e12b4fac
	$(PYTHON) setup.py install

dev: all
	bash install-pplacer.sh
	pip install git+https://github.com/Ecogenomics/CheckM.git@d74bb68d48b2318542eb7137343196d8e12b4fac coverage
	pip install -e .

clean: distclean

distclean: ;
