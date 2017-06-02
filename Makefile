.PHONY: develop docs minimap miniasm racon bwa samtools nanonet

venv: venv/bin/activate

venv/bin/activate:
	test -d venv || virtualenv venv --python=python3
	. ./venv/bin/activate && pip install pip --upgrade
	#touch venv/bin/activate

install: venv externals
	. ./venv/bin/activate && pip install -r requirements.txt && python setup.py install

externals: minimap miniasm racon bwa samtools nanonet fast5_research

minimap: submodules/minimap/minimap
submodules/minimap/minimap:
	cd submodules/minimap && make

miniasm: submodules/miniasm/miniasm
submodules/miniasm/miniasm:
	cd submodules/miniasm && make

racon: submodules/racon/racon
submodules/racon/racon:
	cd submodules/racon && make modules && make tools && make -j
	cp submodules/racon/bin/racon submodules/racon/racon

bwa: submodules/bwa/bwa
submodules/bwa/bwa:
	cd submodules/bwa && make

samtools: submodules/samtools/samtools
submodules/samtools/samtools:
	cd submodules && wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && tar -xjf samtools-1.3.1.tar.bz2
	mv submodules/samtools-1.3.1 submodules/samtools
	cd submodules/samtools && make

fast5_research: venv
	. ./venv/bin/activate && cd submodules/fast5_research && pip install .

nanonet: venv
	. ./venv/bin/activate && pip install numpy
	. ./venv/bin/activate && cd submodules/nanonet && pip install .


docs: venv
	. ./venv/bin/activate && pip install sphinx sphinx_rtd_theme sphinx-argparse
	. ./venv/bin/activate && cd docs && make clean api html
