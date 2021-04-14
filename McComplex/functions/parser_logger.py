from __future__ import print_function
from pysimplelog import Logger
import argparse

parser = argparse.ArgumentParser(description = 
	"This program is able to reconstruct biological macrocomplexes of protein-protein interactions as well as protein-DNA/RNA interactions given a set of binary interactions and the desired number of chains of the target complex.")

requiredNamed = parser.add_argument_group('required arguments')

requiredNamed.add_argument('-i', '--indir',			#INPUT FOLDER argument
							dest = "indir",
							action = "store",
							required=True,
							help = "Input folder (or path) containing all PDB files with the protein binary interactions. It is a required argument.")

parser.add_argument('-nc', '--number_chains',		#NUMBER OF CHAINS argument
							dest = "number_chains",
							action = "store",
							type = int,
							default = 100,
							help = "Number of chains desired for the target complex. This is an optional argument.")

parser.add_argument('-o', '--outdir',			#OUTPUT FOLDER argument
					dest = "outdir",
					action = "store",
					default = None,
					help = "If set, all the models generated in each iteration, the final macrocomplex structure in PDB format and the log file will be saved in this folder. By default, the output folder will be named as the input folder + \"_output\".")

parser.add_argument('-v', '--verbose',			#VERBOSE argument
					dest = "verbose",
					action = "store_true",
					default = False,
					help = "If set, the progression log printed in standard output file.")

parser.add_argument('-pi', '--pdb_iterations',		#PDB FILES ITERATIONS argument
					dest = "pdb_iterations",
					action = "store_true",
					default = False,
					help = "If set, each time a chain is added to the complex, a new PDB file will be saved.")

parser.add_argument('-rmsd', '--rmsd_threshold',		
					dest = "rmsd_threshold",
					action = "store",
					default = 0.3,
					type = float,
					help = "If set, the RMSD threshold for considering a superimposition as correct will take this value. If not, it will be 0.3 by default. The output of the program is very sensitive to this value, we advise to be careful when modifying it.")

parser.add_argument('-cl', '--clashes_threshold',		
					dest = "clashes",
					action = "store",
					default = 30,
					type = int,
					help = "If set, the threshold of the number of clashes will take this value. If not, it will be 30 by default. The output of the program is very sensitive to this value, we advise to be careful when modifying it.")
parser.add_argument('-opt', '--optimisation',     
					dest= "optimisation", 
					action= "store_true", 
					default= False, 
					help = "If set, it runs an optimisation on your output structure and saves it in a separate PDB file. It is only available for structures with less than 99.999 residues"
					)					


l = Logger()