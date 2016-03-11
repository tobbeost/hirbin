HirBin tutorial
=======

This tutorial gives an example on how to use `HirBin` for sub-binning and statistical analysis. 
In this tutorial the supervised annotation of the reference sequences and mapping of the reads to the reference sequences has already been done, but it can also be done using HirBin functions `functionalAnnotation.py` and `mappingReads.py`

Getting started
----------------

Start with having a look at the tutorial data. The example data used here is a subset of the Qin et al. (2012) dataset comparing the gut microbiota of 5 patients with t2d and 5 healthy patients.
The input files are:

Step 1: Prepare the metadata file
----------------------------------

HirBin uses the same metadata file throughout all analysis. In the metadata files the relative file paths to the different files are defined together with groups used for the statistical analysis. Create a file called `metadata.txt` in the tutorial directory and paste the following information:

  name  group  reference  annotation
  DLM005  t2d  tutorial_data/contigs/DLM005.fasta  tutorial_data/TIGRFAM/DLM005.hmmout
  DOM005  t2d  tutorial_data/contigs/DOM005.fasta  tutorial_data/TIGRFAM/DOM005.hmmout
  DOM012  t2d  tutorial_data/contigs/DOM012.fasta  tutorial_data/TIGRFAM/DOM012.hmmout
  DOM017  t2d  tutorial_data/contigs/DOM017.fasta  tutorial_data/TIGRFAM/DOM017.hmmout
  DOM025  t2d  tutorial_data/contigs/DOM025.fasta  tutorial_data/TIGRFAM/DOM025.hmmout
  NLM021  control  tutorial_data/contigs/NLM021.fasta  tutorial_data/TIGRFAM/NLM021.hmmout  
  NLM031  control  tutorial_data/contigs/NLM031.fasta  tutorial_data/TIGRFAM/NLM031.hmmout
  NOM008  control  tutorial_data/contigs/NOM008.fasta  tutorial_data/TIGRFAM/NOM008.hmmout
  NOM017  control  tutorial_data/contigs/NOM017.fasta  tutorial_data/TIGRFAM/NOM017.hmmout
  NOM026  control  tutorial_data/contigs/NOM026.fasta  tutorial_data/TIGRFAM/NOM026.hmmout

Note: Here the reads from each sample have been mapped to sample-specific referent files (per sample assembly). If a common reference is used, then the same file should be repeated in the reference column.

Step 2: Unsupervised clustring (sub-binning)
----------------------------------------------

