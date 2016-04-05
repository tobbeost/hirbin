Hirbin
=======

Hirbin (High-resolution binning) is a python package for binning and sub-binning of metagenomes at the functional level. 
It combines supervised functional annotation using e.g. PFAM, TIGRFAM or COGs, with an unsupervised clustering step 
where the sequences in each bin are clustered using `Uclust`.


Installation
-------------

Hirbin depends on the following python packages:
	
	numpy (1.8.1)
	biopython (1.66)
	multiprocessing (2.6.2.1)


It is recommended to create a new conda environment using [miniconda](http://conda.pydata.org/miniconda.html) and install
the required packages inside:

	conda create -n hirbin numpy biopython
	source activate hirbin

Alternatively install the required python packages first by for example typing `pip install biopython`.

After installing the python dependencies, clone hirbin using git and install the hirbin package:

	git clone https://github.com/tobbeost/hirbin.git
	pip install hirbin/


The four main hirbin programs are:

    functionalAnnotation.py
	clusterBinsToSubbins.py
	mappingReads.py
	statisticalAnalysis.py

After installation type e.g. `functionalAnnotation.py -h` to see the options for the program.


### Install external programs
