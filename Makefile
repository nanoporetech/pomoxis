.PHONY: externals pip_submodules install docs

# Builds a cache of binaries which can just be copied for CI
BINARIES=minimap miniasm racon bwa samtools
BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)


$(BINCACHEDIR)/minimap: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/minimap && make
	cp submodules/minimap/minimap $@

$(BINCACHEDIR)/miniasm: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/miniasm && make
	cp submodules/miniasm/miniasm $@

$(BINCACHEDIR)/racon: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/racon && make modules && make tools && make -j
	cp submodules/racon/bin/racon $@

$(BINCACHEDIR)/bwa: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/bwa && make
	cp submodules/bwa/bwa $@

SAMVER=1.3.1
$(BINCACHEDIR)/samtools: | $(BINCACHEDIR)
	@echo Making $(@F)
	# tar.bz is not a dependency, since that would cause it to be fetched
	#   even when installing from $(BINCACHEDIR)
	if [ ! -e submodules/samtools-${SAMVER}.tar.bz2 ]; then \
	  cd submodules; \
	  wget https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
	fi
	cd submodules && tar -xjf samtools-${SAMVER}.tar.bz2
	cd submodules/samtools-${SAMVER} && make
	cp submodules/samtools-${SAMVER}/samtools $@


venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install numpy # needs to get done before other things

install: venv | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && pip install -r requirements.txt && python setup.py install


# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         =
BUILDDIR      = _build

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

DOCSRC = docs

docs: venv
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll
