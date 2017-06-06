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

install: venv pip_submodules | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && pip install -r requirements.txt && python setup.py install

pip_submodules:
	${IN_VENV} && pip install numpy
	${IN_VENV} && cd submodules/nanonet && pip install .
	${IN_VENV} && cd submodules/fast5_research && pip install .

docs: venv
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd docs && make clean api html
