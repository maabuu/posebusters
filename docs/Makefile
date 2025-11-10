# makefile for Sphinx documentation
SPHINXOPTS  ?=
SPHINXBUILD  ?= sphinx-build
SOURCEDIR = source
BUILDDIR = build

# Use uv to run
SPHINX = uv run sphinx-build

# Build docs
.PHONY: html pdf epub
html:
	$(SPHINX) -b html $(SOURCEDIR) $(BUILDDIR)/html

pdf:
	$(SPHINX) -b latexpdf $(SOURCEDIR) $(BUILDDIR)/pdf

epub:
	$(SPHINX) -b epub $(SOURCEDIR) $(BUILDDIR)/epub

# Clean build artifacts
.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/*
