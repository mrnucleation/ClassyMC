SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help Makefile auto html clean

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Generate auto-docs from Fortran source, then build HTML
html: auto
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo "Build finished. HTML pages are in $(BUILDDIR)/html."

# Run the Fortran source parser to regenerate auto/ RST files
auto:
	@echo "Parsing Fortran source files..."
	@python3 generate_docs.py
	@echo "Auto-documentation generated in auto/"

# Only regenerate auto/ files
generate:
	@python3 generate_docs.py

clean:
	@rm -rf $(BUILDDIR) auto/

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
