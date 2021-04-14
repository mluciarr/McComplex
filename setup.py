#!/usr/bin/env python

from distutils.core import setup
import setuptools 

setup(name='McComplex',
	version='1.0',
	description="""This program reconstructs macrocomplexes of protein-protein 
	and protein-(DNA/RNA) from a list of files of binary interactions of its chains""",
	author='Maria Lucía Romero, Ferran Pegenaute, Ipek Yaren',
	author_email='ferran.pegenaute01@estudiant.upf.edu',
	long_description=open('README.md').read(),
	install_requires=['biopython >= 1.73.0','argparse >= 1.1.0', 'pysimplelog'],
	packages=['McComplex', 'McComplex.functions'],
	license='LICENSE.txt',
	url='https://github.com/ferranpgp/McCrocomplex')
