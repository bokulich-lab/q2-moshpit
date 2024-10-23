.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	python -m coverage run --omit "**/tests*,**/__init__.py,**/_version.py,versioneer.py" -m pytest && coverage xml -o coverage.xml

install: all
	$(PYTHON) setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;
