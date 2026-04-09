SHELL=/bin/bash
OS := $(shell uname | tr '[:upper:]' '[:lower:]')
CXX         ?= g++
PYTHON      ?= python3

# Builds a cache of binaries which can just be copied for CI
BINARIES=minimap2 samtools bcftools paftools.js seqkit bedtools bgzip tabix k8

BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)


BINBUILDDIR=binbuild
$(BINBUILDDIR):
	mkdir -p $(BINBUILDDIR)

#MAPVER=2.23 upgraded to 2.30 on 5 December 2025
MAPVER=2.30
$(BINCACHEDIR)/minimap2: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/minimap2-${MAPVER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/lh3/minimap2/releases/download/v${MAPVER}/minimap2-${MAPVER}.tar.bz2; \
	  tar -xjf minimap2-${MAPVER}.tar.bz2; \
	fi
	cd ${BINBUILDDIR}/minimap2-${MAPVER} && make
	cp ${BINBUILDDIR}/minimap2-${MAPVER}/minimap2 $@

#SAMVER=1.14 upgraded to 1.22 on 5 December 2025
SAMVER=1.22
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

$(BINCACHEDIR)/paftools.js: | $(BINCACHEDIR)/minimap2
	cp ${BINBUILDDIR}/minimap2-${MAPVER}/misc/$(@F) $@

# K8VER=0.2.5 upgraded to v1.2 on 5 January 2026
# k8 is required for paftools.js
K8VER=1.2
$(BINCACHEDIR)/k8: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/k8-${K8VER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/attractivechaos/k8/releases/download/v${K8VER}/k8-${K8VER}.tar.bz2; \
	  tar -xjf k8-${K8VER}.tar.bz2; \
	fi
	cp ${BINBUILDDIR}/k8-${K8VER}/k8*Linux $@

#BCFVER=1.14 upgraded to 1.22 on 5 December 2025
BCFVER=1.22
$(BINCACHEDIR)/bcftools: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/bcftools-${BCFVER}.tar.bz2 ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/samtools/bcftools/releases/download/${BCFVER}/bcftools-${BCFVER}.tar.bz2; \
	  tar -xjf bcftools-${BCFVER}.tar.bz2; \
	fi
	cd ${BINBUILDDIR}/bcftools-${BCFVER} && make
	cp ${BINBUILDDIR}/bcftools-${BCFVER}/bcftools $@


#SEQKITVER=2.1.0 upgraded to 2.12.0 on 10 December 2025
SEQKITVER=2.12.0
$(BINCACHEDIR)/seqkit: | $(BINCACHEDIR) $(BINBUILDDIR)
	@echo Making $(@F)
	if [ ! -e ${BINBUILDDIR}/seqkit_${OS}_amd64.tar.gz ]; then \
	  cd ${BINBUILDDIR}; \
	  wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKITVER}/seqkit_${OS}_amd64.tar.gz; \
	  tar -xzvf seqkit_${OS}_amd64.tar.gz; \
	fi
	cp ${BINBUILDDIR}/seqkit $@	


#BEDTOOLSVER=2.12.0 upgraded to 2.31.1 on 5 Dec 2025
BEDTOOLSVER=2.31.1
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
	test -d venv || $(PYTHON) -m venv venv --prompt 'pomoxis'
	${IN_VENV} && pip install pip --upgrade
	${IN_VENV} && pip install -e .


.PHONY: install
install: venv $(addprefix $(BINCACHEDIR)/, $(BINARIES)) | copy_binaries

test: install
	${IN_VENV} && python -m unittest discover

# We copy the binaries in Make as we don't officially package them in the wheel
.PHONY: copy_binaries
copy_binaries: venv $(addprefix $(BINCACHEDIR)/, $(BINARIES))
	cp $(addprefix $(BINCACHEDIR)/, $(BINARIES)) venv/bin/


# Building and distributing wheels is done from a separate build environment
.PHONY: build
build: venv_pypi/bin/activate
IN_BUILD=. ./venv_pypi/bin/activate
venv_pypi/bin/activate:
	test -d venv_pypi || $(PYTHON) -m venv venv_pypi --prompt "pypi"
	${IN_BUILD} && pip install --upgrade pip
	${IN_BUILD} && pip install .[packaging]

.PHONY: sdist
sdist: venv_pypi/bin/activate
	${IN_BUILD} && python -m build
