#!/usr/bin/env python3

from distutils.core import setup

setup(name='cnp',
      version='0.1',
      description='Assorted utilities for working with chromatic planar graphs (fun!)',
      author='Oisin Carroll',
      author_email='hi@imois.in',
      install_requires=[
          "python-sat",
          "sympy",
          "networkx",
          "scipy",
          "matplotlib"
      ],
      # url='',
      packages=['cnp'],
     )
