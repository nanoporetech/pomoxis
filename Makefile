SHELL=/bin/bash
OS := $(shell uname | tr '[:upper:]' '[:lower:]')

# for porechop on travis (or other platform with older gcc)
CXX         ?= g++
PYTHON      ?= python3
CONDA       ?= ~/miniconda3/

# Builds a cache of binaries which can just be copied for CI
BINARIES=minimap2 miniasm racon samtools bcftools k8 paftools.js seqkit bedtools bgzip tabix

BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)


BINBUILDDIR=binbuild
$(BINBUILDDIR):
	mkdir -p $(BINBUILDDIR)


#MAPVER=2.14  Changed to 2.23 on 6 December 2021
MAPVER=2.23
$(BINCACHEDIR)/minimap2: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/minimap2-${MAPVER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/lh3/minimap2/releases/download/v${MAPVER}/minimap2-${MAPVER}.tar.bz2; \
	  tar -xjf minimap2-${MAPVER}.tar.bz2; \
	fi
	cd ${BINBUILDDIR}/minimap2-${MAPVER} && make
	cp ${BINBUILDDIR}/minimap2-${MAPVER}/minimap2 $@


ASMVER=0.3
$(BINCACHEDIR)/miniasm: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/miniasm-v${ASMVER}.tar.gz ]; then \
	  cd ${BINBUILDDIR}; \
	  wget -O miniasm-v${ASMVER}.tar.gz https://github.com/lh3/miniasm/archive/v${ASMVER}.tar.gz; \
	  tar -xzf miniasm-v${ASMVER}.tar.gz; \
	fi
	cd ${BINBUILDDIR}/miniasm-${ASMVER} && make
	cp ${BINBUILDDIR}/miniasm-${ASMVER}/miniasm $@


RACONVER=1.4.13
$(BINCACHEDIR)/racon: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	@echo GCC is $(GCC)
	if [ ! -e ${BINBUILDDIR}/racon-v${RACONVER}.tar.gz ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/lbcb-sci/racon/releases/download/${RACONVER}/racon-v${RACONVER}.tar.gz; \
	  tar -xzf racon-v${RACONVER}.tar.gz; \
	fi
	cd ${BINBUILDDIR}/racon-v${RACONVER} && mkdir build && cd build && \
		cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-mno-avx2 " ..
	cd ${BINBUILDDIR}/racon-v${RACONVER}/build && make
	cp ${BINBUILDDIR}/racon-v${RACONVER}/build/bin/racon $@


#SAMVER=1.8  Changed to 1.14 on 6 December 2021
SAMVER=1.14
$(BINCACHEDIR)/samtools: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	# tar.bz is not a dependency, since that would cause it to be fetched
	#   even when installing from $(BINCACHEDIR)
	if [ ! -e ${BINBUILDDIR}/samtools-${SAMVER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/samtools/samtools/releases/download/${SAMVER}/samtools-${SAMVER}.tar.bz2; \
	  tar -xjf samtools-${SAMVER}.tar.bz2; \
	fi
	# make all-htslib to get bgzip and tabix
	cd ${BINBUILDDIR}/samtools-${SAMVER} && make all all-htslib
	cp ${BINBUILDDIR}/samtools-${SAMVER}/samtools $@

$(BINCACHEDIR)/tabix: | $(BINCACHEDIR)/samtools
	cp ${BINBUILDDIR}/samtools-${SAMVER}/htslib-${SAMVER}/$(@F) $@


$(BINCACHEDIR)/bgzip: | $(BINCACHEDIR)/samtools
	cp ${BINBUILDDIR}/samtools-${SAMVER}/htslib-${SAMVER}/$(@F) $@

K8VER=0.2.5
$(BINCACHEDIR)/k8: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/k8-${K8VER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/attractivechaos/k8/releases/download/${K8VER}/k8-${K8VER}.tar.bz2; \
	  tar -xjf k8-${K8VER}.tar.bz2; \
	fi
	cp ${BINBUILDDIR}/k8-${K8VER}/k8-Linux $@


$(BINCACHEDIR)/paftools.js: | $(BINCACHEDIR)/minimap2
	cp ${BINBUILDDIR}/minimap2-${MAPVER}/misc/$(@F) $@


#BCFVER=1.7  changed to 1.14 on 6 December 2021
BCFVER=1.14
$(BINCACHEDIR)/bcftools: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/bcftools-${BCFVER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/samtools/bcftools/releases/download/${BCFVER}/bcftools-${BCFVER}.tar.bz2; \
	  tar -xjf bcftools-${BCFVER}.tar.bz2; \
	fi
	cd ${BINBUILDDIR}/bcftools-${BCFVER} && make
	cp ${BINBUILDDIR}/bcftools-${BCFVER}/bcftools $@


#SEQKITVER=0.8.0 changed to 2.1.0 on 13 December 2021
SEQKITVER=2.1.0
$(BINCACHEDIR)/seqkit: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/seqkit_${OS}_amd64.tar.gz ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKITVER}/seqkit_${OS}_amd64.tar.gz; \
	  tar -xzvf seqkit_${OS}_amd64.tar.gz; \
	fi
	cp ${BINBUILDDIR}/seqkit $@	


#BEDTOOLSVER=2.29.0 Changed to 2.30.0 on 6 December 2021
BEDTOOLSVER=2.30.0
$(BINCACHEDIR)/bedtools: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/bedtools-2.29.0.tar.gz	]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLSVER}/bedtools-${BEDTOOLSVER}.tar.gz; \
	  mkdir bedtools-${BEDTOOLSVER}; \
	  tar -xzvf bedtools-${BEDTOOLSVER}.tar.gz --directory bedtools-${BEDTOOLSVER}; \
	fi
	cd ${BINBUILDDIR}/bedtools-${BEDTOOLSVER}/bedtools2/ && make
	cp ${BINBUILDDIR}/bedtools-${BEDTOOLSVER}/bedtools2/bin/bedtools $@


.PHONY: venv
venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate
venv/bin/activate:
	test -d venv || $(PYTHON) -m venv venv --prompt '(pomoxis) '
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -r requirements.txt


.PHONY: install
install: venv | $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	${IN_VENV} && POMO_BINARIES=1 python setup.py install


.PHONY: build
build: pypi_build/bin/activate
IN_BUILD=. ./pypi_build/bin/activate
pypi_build/bin/activate:
	test -d pypi_build || $(PYTHON) -m venv pypi_build --prompt "(pypi) "
	${IN_BUILD} && pip install pip --upgrade
	${IN_BUILD} && pip install --upgrade pip setuptools twine wheel readme_renderer[md]


.PHONY: sdist
sdist: pypi_build/bin/activate
	${IN_BUILD} && python setup.py sdist


# You can set these variables from the command line.
SPHINXOPTS  =
SPHINXBUILD = sphinx-build
PAPER       =
BUILDDIR    = _build
DOCSRC      = docs

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: docs
docs: venv
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && ./prog_docs.py > $(DOCSRC)/programs.rst
	${IN_VENV} && cd $(DOCSRC) && $(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(DOCSRC)/$(BUILDDIR)/html."
	touch $(DOCSRC)/$(BUILDDIR)/html/.nojekyll
