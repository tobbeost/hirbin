Hirbin
=======

Hirbin (High-resolution binning) is a python package for binning and sub-binning of metagenomes at the functional level. 
It combines supervised functional annotation using e.g. PFAM, TIGRFAM or COGs, with an unsupervised clustering step 
where the sequences in each bin are clustered using `Uclust`.

The four main hirbin programs are:

    functionalAnnotation.py
   	mappingReads.py
    clusterBinsToSubbins.py
  	statisticalAnalysis.py

Documentation
-------------
[HirBin documentation](https://github.com/cmbio/hirbin/wiki)


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

After installing the python dependencies, download the hirbin source code, extract the files and install the hirbin package:

	wget https://github.com/cmbio/hirbin/archive/v0.1.tar.gz
	pip install v0.1.tar.gz/

After installation type e.g. `functionalAnnotation.py -h` to see the options for the program.



### Install external programs

Depending on which part of Hirbin you want to use you need to install some external
programs and add them to the path environmental variable, so that Hirbin can find them.

The external programs needed are:
* [hmmer3](http://hmmer.org/) (used by `functionalAnnotation.py`)
* [EMBOSS transeq v6.5.7](ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/) (used by `functionalAnnotation.py`)
* [USEARCH (uclust) ](http://www.drive5.com/usearch/download.html) (used by `clusterBinsToSubbins.py`)
* [R](https://cran.r-project.org/) (used by `statisticalAnalysis.py`)
* Bowtie2 (used by mappingReads.py)
* bedtools (used by mappingReads.py)
* samtools (used by mappingReads.py)
