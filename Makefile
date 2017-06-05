.PHONY: pip_submodules install docs

# Builds a cache of binaries which can just be copied for CI
BINARIES=minimap miniasm racon bwa samtools
BINCACHEDIR=bincache
$(BINCACHEDIR):
	mkdir -p $(BINCACHEDIR)
define build_rule
    $(BINCACHEDIR)/$1: $(BINCACHEDIR) submodules/$1/$1
	cp submodules/$1/$1 $(BINCACHEDIR)/$1
endef
$(foreach f,$(BINARIES),$(eval $(call build_rule,$f)))
CACHEDBIN=$(addprefix $(BINCACHEDIR)/, $(BINARIES))



submodules/minimap/minimap:
	cd submodules/minimap && make

submodules/miniasm/miniasm:
	cd submodules/miniasm && make

submodules/racon/racon:
	cd submodules/racon && make modules && make tools && make -j
	cp submodules/racon/bin/racon $@

submodules/bwa/bwa:
	cd submodules/bwa && make

submodules/samtools/samtools:
	cd submodules && wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && tar -xjf samtools-1.3.1.tar.bz2
	mv submodules/samtools-1.3.1 submodules/samtools
	cd submodules/samtools && make


venv: venv/bin/activate
IN_VENV=. ./venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=python3
	${IN_VENV} && pip install pip --upgrade

install: venv pip_submodules | $(CACHEDBIN)
	${IN_VENV} && pip install -r requirements.txt && python setup.py install

pip_submodules:
	${IN_VENV} && pip install numpy
	${IN_VENV} && cd submodules/nanonet && pip install .
	${IN_VENV} && cd submodules/fast5_research && pip install .

docs: venv
	${IN_VENV} && pip install sphinx sphinx_rtd_theme sphinx-argparse
	${IN_VENV} && cd docs && make clean api html
