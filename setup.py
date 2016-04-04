#!/usr/bin/env python
from setuptools import setup

setup(name='hirbin',
      version='0.1',
      description='High-resolution binning (hirbin)',
      long_description = open('README.md').read(),
      url='http://github.com/tobbeost/hirbin',
      author='Tobias Osterlund',
      author_email='tobiaso@chalmers.se',
      license='MIT',
      packages=['hirbin'],
      install_requires=['Biopython', 
                        'multiprocessing'],
      classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['hirbin/functionalAnnotation.py', 
               'hirbin/clusterBinsToSubbins.py', 
               'hirbin/statisticalAnalysis.py'
              ],
      zip_safe=False)
