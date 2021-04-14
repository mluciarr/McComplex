#!/usr/bin/env python

from setuptools import setup

setup(name='McCrocomplex',
	version='1.0',
	description='This program is able to reconstruct biological macrocomplexes of protein-protein interactions as well as protein-DNA/RNA interactions given a set of binary interactions and the desired number of chains of the target complex.',
	author='Maria Lucía Romero, Ferran Pegenaute, Ipek Yaren',
	author_email='ferran.pegenaute01@estudiant.upf.edu',
	long_description=open('README.md').read(),
	install_requires=['biopython >= 1.73.0','argparse >= 1.1.0', 'pysimplelog' ],
	license='LICENSE.txt',
	url='https://github.com/ferranpgp/McCrocomplex',
	scripts='McComplex_builder.py')
