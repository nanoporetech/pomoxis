.PHONY: install docs
SHELL=/bin/bash
OS := $(shell uname | tr '[:upper:]' '[:lower:]')

# for porechop on travis (or other platform with older gcc)
CXX         ?= g++

CONDA?=~/miniconda3/

# Builds a cache of binaries which can just be copied for CI
BINARIES=minimap2 miniasm bwa racon samtools bcftools seqkit

BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)

$(BINCACHEDIR)/minimap2: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/minimap2 && make
	cp submodules/minimap2/minimap2 $@

$(BINCACHEDIR)/miniasm: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/miniasm && make
	cp submodules/miniasm/miniasm $@

$(BINCACHEDIR)/racon: | $(BINCACHEDIR)
	@echo Making $(@F)
	@echo GCC is $(GCC)
	#cd submodules/racon && make modules && make -j $(RACONOS)
	#cp submodules/racon/bin/racon$(RACONTAG) $@
	cd submodules/racon && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
	cd submodules/racon/build && make
	cp submodules/racon/build/bin/racon $@

$(BINCACHEDIR)/bwa: | $(BINCACHEDIR)
	@echo Making $(@F)
	cd submodules/bwa && make
	cp submodules/bwa/bwa $@

SAMVER=1.8
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

BCFVER=1.7
$(BINCACHEDIR)/bcftools: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -e submodules/bcftools-${BCFVER}.tar.bz2 ]; then \
	  cd submodules; \
	  wget https://github.com/samtools/bcftools/releases/download/${BCFVER}/bcftools-${BCFVER}.tar.bz2; \
	fi
	cd submodules && tar -xjf bcftools-${BCFVER}.tar.bz2
	cd submodules/bcftools-${BCFVER} && make
	cp submodules/bcftools-${BCFVER}/bcftools $@

SEQKITVER=0.8.0
$(BINCACHEDIR)/seqkit: | $(BINCACHEDIR)
	@echo Making $(@F)
	if [ ! -e submodules/seqkit_${OS}_amd64.tar.gz ]; then \
	  cd submodules; \
	  wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKITVER}/seqkit_${OS}_amd64.tar.gz; \
	fi
	cd submodules && tar -xzvf seqkit_${OS}_amd64.tar.gz
	cp submodules/seqkit $@	

venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

PYVER=3.6
venv/bin/activate:
	test -d venv || virtualenv venv --prompt '(pomoxis) ' --python=python${PYVER}
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt

bwapy: venv
	cd submodules/bwapy && make bwa/libbwa.a 
	${IN_VENV} && cd submodules/bwapy && python setup.py install

install: venv bwapy | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && POMO_BINARIES=1 python setup.py install

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
