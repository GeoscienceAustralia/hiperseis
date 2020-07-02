# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	sphinx-apidoc -f --ext-autodoc --ext-coverage --ext-viewcode --ext-githubpages -o . ../seismic  ../legacy \
	    ../seismic/**/sandbox ../seismic/ml_classifier ../seismic/inventory/legacy \
	    ../seismic/**/example_*.py ../seismic/ASDFdatabase/ASDF_build*.py ../seismic/ASDFdatabase/minimus*.py

	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)