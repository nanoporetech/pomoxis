.PHONY: install docs
SHELL=/bin/bash
OS := $(shell uname | tr '[:upper:]' '[:lower:]')

# for porechop on travis (or other platform with older gcc)
CXX         ?= g++

CONDA?=~/miniconda3/

# Builds a cache of binaries which can just be copied for CI
BINARIES=minimap2 miniasm racon samtools bcftools seqkit

BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)

MAPVER=2.14
$(BINCACHEDIR)/minimap2: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -e submodules/minimap2-${MAPVER}.tar.bz2 ]; then \
	  cd submodules; \
	  wget https://github.com/lh3/minimap2/releases/download/v${MAPVER}/minimap2-${MAPVER}.tar.bz2; \
	  tar -xjf minimap2-${MAPVER}.tar.bz2; \
	fi
	cd submodules/minimap2-${MAPVER} && make
	cp submodules/minimap2-${MAPVER}/minimap2 $@

ASMVER=0.3
$(BINCACHEDIR)/miniasm: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -e submodules/miniasm-v${ASMVER}.tar.gz ]; then \
	  cd submodules; \
	  wget -O miniasm-v${ASMVER}.tar.gz https://github.com/lh3/miniasm/archive/v${ASMVER}.tar.gz; \
	  tar -xzf miniasm-v${ASMVER}.tar.gz; \
	fi
	cd submodules/miniasm-${ASMVER} && make
	cp submodules/miniasm-${ASMVER}/miniasm $@

RACONVER=1.3.1
$(BINCACHEDIR)/racon: | $(BINCACHEDIR)
	@echo Making $(@F)
	@echo GCC is $(GCC)
	if [ ! -e submodules/racon-v${RACONVER}.tar.gz ]; then \
	  cd submodules; \
	  wget https://github.com/isovic/racon/releases/download/${RACONVER}/racon-v${RACONVER}.tar.gz; \
	  tar -xzf racon-v${RACONVER}.tar.gz; \
	fi
	cd submodules/racon-v${RACONVER} && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
	cd submodules/racon-v${RACONVER}/build && make
	cp submodules/racon-v${RACONVER}/build/bin/racon $@

SAMVER=1.8
$(BINCACHEDIR)/samtools: | $(BINCACHEDIR)
	@echo Making $(@F)
	# tar.bz is not a dependency, since that would cause it to be fetched
	#   even when installing from $(BINCACHEDIR)
	if [ ! -e submodules/samtools-${SAMVER}.tar.bz2 ]; then \
	  cd submodules; \
	  wget https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
	  tar -xjf samtools-${SAMVER}.tar.bz2; \
	fi
	cd submodules/samtools-${SAMVER} && make
	cp submodules/samtools-${SAMVER}/samtools $@

BCFVER=1.7
$(BINCACHEDIR)/bcftools: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -e submodules/bcftools-${BCFVER}.tar.bz2 ]; then \
	  cd submodules; \
	  wget https://github.com/samtools/bcftools/releases/download/${BCFVER}/bcftools-${BCFVER}.tar.bz2; \
	  tar -xjf bcftools-${BCFVER}.tar.bz2; \
	fi
	cd submodules/bcftools-${BCFVER} && make
	cp submodules/bcftools-${BCFVER}/bcftools $@

SEQKITVER=0.8.0
$(BINCACHEDIR)/seqkit: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -e submodules/seqkit_${OS}_amd64.tar.gz ]; then \
	  cd submodules; \
	  wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKITVER}/seqkit_${OS}_amd64.tar.gz; \
	  tar -xzvf seqkit_${OS}_amd64.tar.gz; \
	fi
	cp submodules/seqkit $@	

venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --prompt '(pomoxis) ' --python=python3
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

bwapy: venv
	cd submodules/bwapy && make bwa/libbwa.a 
	${IN_VENV} && cd submodules/bwapy && python setup.py install

install: venv bwapy | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && POMO_BINARIES=1 python setup.py install

PYVER=3.6
IN_CONDA=. ${CONDA}/etc/profile.d/conda.sh
conda:
	${IN_CONDA} && conda remove -n pomoxis --all
	${IN_CONDA} && conda create -y -n pomoxis -c bioconda -c conda-forge porechop \
		samtools=${SAMVER} bcftools=${BCFVER} seqkit=${SEQKITVER} \
		miniasm=${ASMVER} minimap2=${MAPVER} racon=${RACONVER} \
		python=${PYVER}
	grep -v Porechop requirements.txt > conda_reqs.txt
	${IN_CONDA} && conda activate pomoxis && pip install -r conda_reqs.txt
	${IN_CONDA} && conda activate pomoxis && python setup.py install \
		--single-version-externally-managed --record=conda_install.out
	rm conda_reqs.txt

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
