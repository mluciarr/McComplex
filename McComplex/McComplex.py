import Bio.PDB
import sys
import string
import os
import argparse
import timeit
import logging
import re
from functions.parser_logger import *
from functions.McComplex_functions import *
from functions.class_error import *
from functions.optimization import *


beginning = timeit.default_timer()
arguments = parser.parse_args()


# Mandatory arguments 

try:
	if os.path.exists(arguments.indir):
		files = sorted(list(filter(lambda x: x.endswith(".pdb"), os.listdir(arguments.indir))))
		arguments.indir = os.path.abspath(arguments.indir)
except:
	raise InputError(arguments.indir)

os.chdir(arguments.indir + "/../")

# Optional arguments

if arguments.outdir == None:		
	arguments.outdir = arguments.indir + "_output"		
else:
	pass
if not os.path.exists(arguments.outdir):		
	os.mkdir(arguments.outdir)	
	arguments.outdir = os.path.abspath(arguments.outdir)
else:
	arguments.outdir = os.path.abspath(arguments.outdir)

# Change to output directory
os.chdir(arguments.outdir)		

# Setting up the Log system
l.set_log_file_basename("McComplex_mylog")

	
if arguments.verbose:
	l.set_log_to_stdout_flag(True)
else:
	l.set_log_to_stdout_flag(False)

l.debug('...STARTING...')

# Checking files

if re.search("\d",files[0]):		
	files = sorted(files, key=lambda x: str("".join([i for i in x if i.isdigit()])))
else:
	files = sorted(files)

l.info("""Parameters used are:\n  - Number of chains: %d\n  - RMSD threshold: 
%.4f\n  - Clashes threshold %d""" % (arguments.number_chains, 
arguments.rmsd_threshold, arguments.clashes))

# Using the first file as the reference for the McComplex 	
file_path = arguments.indir + "/" + files[0]

# Read PDB. First file is assing as reference structure.
pdb_parser = Bio.PDB.PDBParser(QUIET = True)	
ref_structure = pdb_parser.get_structure("reference", file_path)	
l.info("The initial complex has %d chains and are the following:" % (ref_structure[0].__len__()))
for ID in [chain.get_id() for chain in ref_structure[0].get_chains()]:		
	l.info("Chain %s", ID)		#prints the ID

# Calling the main function for the first time.
ref_structure = McComplex(ref_structure = ref_structure, 
files_list = files, iteration = 0, parsed_no_chains = 0, command_arguments = arguments)

# Function finished. Output formats in function of the lenght

if len(list(ref_structure[0].get_atoms())) > 99999 or len(list(ref_structure[0].get_chains())) > 62:		
	io = Bio.PDB.MMCIFIO()								
	io.set_structure(ref_structure[0])					
	io.save("McComplex.cif")							
	l.info("Output files %s saved in %s" %("""McComplex.cif and 
	McComplex.log""",os.path.abspath(arguments.outdir)))
else: 													
	io = Bio.PDB.PDBIO()								
	io.set_structure(ref_structure[0])					
	io.save("McComplex.pdb")							
	l.info("Output files %s saved in %s" %("""McComplex.pdb and 
	McComplex.log""",os.path.abspath(arguments.outdir)))

	if arguments.optimisation:
		try:
			Optimizemodel("McComplex.pdb")
		except Exception:
			l.warn("Optimization not available, check Modeller")
			pass


end = timeit.default_timer()
l.info("The program has finished. It took %f seconds" % (end - beginning))
