#!/usr/bin/env python
from setuptools import setup

setup(name='hirbin',
      version='0.1',
      description='High-resolution binning (hirbin)',
      long_description = open('README.md').read(),
      url='http://github.com/cmbio/hirbin',
      author='Tobias Osterlund',
      author_email='tobiaso@chalmers.se',
      packages=setuptools.find_packages(),
      package_data={'hirbin': ['README.md', 'hirbin/scripts/statistical_analysis.R']
                   },
      install_requires=['Biopython', 
                        'multiprocessing'],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['hirbin/functionalAnnotation.py', 
               'hirbin/clusterBinsToSubbins.py', 
               'hirbin/statisticalAnalysis.py',
               'hirbin/mappingReads.py'
              ],
      zip_safe=False)
