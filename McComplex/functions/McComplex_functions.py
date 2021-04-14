import Bio.PDB
import sys
import string
import os
import argparse
import timeit
import logging
import re
import sys
import imp
sequence_data = imp.load_source('sequence_data', 'data/sequence_data.py')
import subprocess
import sys
from McComplex.functions.parser_logger import l



def get_key_atoms(chain):
	"""
	This function returns the important atoms atoms from an input chain for 
	superimposition:
	     CA: proteins  
		 C4': nucleic acids,
	And a string with the type of 
	molecule: DNA, RNA or Protein
	"""
	# Lists with all possible nucleic acids letters
	nucleic_acids = sequence_data.nucleic_acids_PDB		
	RNA = sequence_data.RNA_PDB											
	DNA = sequence_data.DNA_PDB						
	atoms = []

	# Loop through all residues of the chain, check its type and append the C4' 
	# or Calphas
	for residue in chain:		
		residue_name = residue.get_resname()[0:3].strip()		
		if residue.get_id()[0] == " " and residue_name not in nucleic_acids:		
			if 'CA' not in residue:
				l.warning("Residue %d %s does not have CAlpha \
					atom" % (residue.get_id()[1], residue_name))
			else:
				atoms.append(residue['CA'])
				molecule_type = 'PROTEIN'	

		## Append the C4's and set the molecule type to DNA or RNA 
		elif residue.get_id()[0] == " " and residue_name in nucleic_acids:	
			if residue_name in DNA:			
				molecule_type = 'DNA'		
			elif residue_name in RNA:		
				molecule_type = 'RNA'		
			atoms.append(residue['C4\''])	
	return(atoms, molecule_type)		


def rmsd_best(i,best_RMSD,prev_RMSD, ref_atoms, query_atoms, rmsd_threshold):
	"""
	Superimpose lists of atoms and retrieve the best RMSD. Returs the best RMSD 
	and a superimposition instance, along with the rest of the input variables 
	to preserve recursivity
	"""
	sup = Bio.PDB.Superimposer()
	sup.set_atoms(ref_atoms, query_atoms)
	RMSD = sup.rms  
	if RMSD > rmsd_threshold:
		pass
	if prev_RMSD is True or RMSD < prev_RMSD:
		best_RMSD = RMSD  
		prev_RMSD = RMSD
	return(i, best_RMSD,prev_RMSD, ref_atoms, query_atoms, rmsd_threshold, sup)




def dict_superimposition(ref_structure, query_structure, rmsd_threshold,i):
	"""
	This funcion returns:
		superimpositions_dict: dictionary of the chain superimposed IDs as a keys 
							  and their superimposition instance as values. 
							  Sorted by rmsd
		superimposed_chains: boolean variables-> True: presence of superimposition; 
												 False: No superimposition
	 	best_RMSD: values of the best RMSD
	 	i: counter of iteration of the recursive function. 
	"""
	
	ref_model = ref_structure[0]
	query_model = query_structure[0]
	superimpositions_dict = {}
	best_RMSD = 0  
	prev_RMSD = True  # To know if we are in the first pair of chains
	superimposed_chains = False
	for ref_chain in ref_model:
		print("Processing reference chain %s" % ref_chain.id)
		l.info("Processing reference chain %s" % ref_chain.id)
		ref_atoms, ref_molecule_type = get_key_atoms(ref_chain)
		## loops through all chains in the sample model ##
		for query_chain in query_model:
			print("Processing query chain %s" % query_chain.id)
			l.info("Processing query chain %s" % query_chain.id)
			query_atoms, query_molecule_type = get_key_atoms(query_chain)
			# Check for equal molecule type and length, superimpose if both true
			if ref_molecule_type != query_molecule_type:
				l.warning(f"Cannot superimpose. Chain {ref_chain.id} is {ref_molecule_type} and {query_chain.id} is {query_molecule_type}")
			elif len(query_atoms) != len(ref_atoms):  
				l.warning(f"Cannot superimpose. Chain {ref_chain.id} and {query_chain.id} differ in length")
			else:
				i,best_RMSD,prev_RMSD,ref_atoms,query_atoms, rmds_threshold,sup\
					= rmsd_best(i,best_RMSD,prev_RMSD, 
					ref_atoms, query_atoms, 0.3)
				superimpositions_dict[(ref_chain.id, query_chain.id)] = sup
				superimposed_chains = True
				l.info(f"The RMSD between fixed chain {ref_chain.id} and query chain {query_chain.id} is {sup.rms}""")
	if superimposed_chains is True:
		superimpositions_dict =sorted(superimpositions_dict.items(),key=lambda x:x[1].rms)

	return(superimpositions_dict, superimposed_chains, best_RMSD, i)

def id_maker(IDs, ID):
	"""
	This function returns a new string ID from a list of IDs. It is limited to
	748 new IDs
	"""
	upperc = list(string.ascii_uppercase)	
	lowerc = list(string.ascii_lowercase)
	digit = list(string.digits)
	alphabet = upperc + lowerc + digit	

	if len(IDs) < 62:				
		if ID not in IDs:			
			return ID
		elif ID in IDs:				
			for i in range(0, len(alphabet)): 
				if alphabet[i] not in IDs:			
					return alphabet[i]
				else:								
					continue
	elif len(IDs) >= 62:	
		for char1 in alphabet:
			for char2 in alphabet:
				ID = char1 + char2	
				if ID not in IDs:	
					return ID
				else:				
					continue


def save_iterations(ref_structure,outdir):
	""" 
	This function produces a PDB file or a mmCIF File depending on the number of
	atoms.Requires a structure instance and the output directory
	"""
	if len(list(ref_structure[0].get_atoms())) > 99999 or \
		len(list(ref_structure[0].get_chains())) > 62:
						io = Bio.PDB.MMCIFIO()
						io.set_structure(ref_structure[0])
						io.save("It_McComplex_%d_chains.cif" %(ref_structure[0].__len__()))	
						l.info("saving It_McComplex_%d_chains.cif in %s" %\
							(ref_structure[0].__len__(),os.path.abspath(outdir)))
	else:
		io = Bio.PDB.PDBIO()
		io.set_structure(ref_structure[0])
		io.save("It_McComplex_%d_chains.pdb" %(ref_structure[0].__len__()))
		l.info("Saving It_McComplex_%d_chains.pdb in %s" %(ref_structure[0].__len__(),os.path.abspath(outdir)))

def update_loop(files_list, i, n):
	"""
	This function updates the i and n counter variables and the files list in 
	each iteration.
	It sums up 1 to the i and n counter vairables and substract the first file 
	of the files list and adds that file at the end of the list. 
	"""
	file = files_list.pop(0)
	files_list.append(file)
	i += 1
	n += 1
	return(files_list, i, n)


def McComplex(ref_structure, files_list, iteration,
	parsed_no_chains, command_arguments):
	"""
	This function superimposes the most similar chain of a binary interaction 
	PDB file with a reference structure and adds the rotated chain to the 
	building complex. Itis a recursive function

	ref_structure: structure instance on which the Mccomplex is going to get 
	build

	files_list: A list off all the binary interaction pdb files storing 
					   the different subunits or chains that form the complex
	iteration: Counter of iterations. Int
	parsed_no_chains: this is a counter that keeps track of the files that have 
	been parsed but no chains were added to the complex. Int
	command_arguments:
			RMSD: this is the RMSD threshold. If the RMSD of a 
			superimposition > than this = wrong superimposition. float

			clashes: this is the clashes or contacts threshold. number of 
			clashes > The chain in question will be dismissed. int

			number_chains: this is the numbers of chains that the complex 
			must have. int

			indir: Input directory relative path. str

			outdir: Output directory relative path. str

			pdb_iterations(boolean):Defined by the user in the command line 
			args, by default False. It will keep a PDB of the reference 
			structure at each iteration. boolean 
	"""
	
	iteration_count = iteration 		
	n = parsed_no_chains
	n_of_chains = command_arguments.number_chains
	clashes_threshold = command_arguments.clashes
	RMSD_threshold = command_arguments.rmsd_threshold
	indir = command_arguments.indir
	outdir = command_arguments.outdir
	pdb_iterations = command_arguments.pdb_iterations

	chains = ref_structure[0].__len__()
	# Prints in McComplex.log file the iteration and number of chains
	l.info("This is the iteration #%d of the recursive function" % iteration_count )
	l.info("The complex has %d chains at this point" % chains)

	
	if chains == n_of_chains or n > len(files_list): 
		l.info("The whole macrocomplex has been build")
		l.info("The final complex has %d chains" % chains)
		l.info("We have arrived to iteration %d" %(iteration_count))
		return 	ref_structure	# Final end or the recursive function

	# The first file of the list is selected to be analyzed. 
	sample = files_list[0]		
	l.info("We are processing the file %s" % (sample))

	# Path of the query file
	file_path = indir + "/" + sample 	
	pdb_parser = Bio.PDB.PDBParser(QUIET = True)
	#saves the Structure object of the sample PDBParser object
	query_structure = pdb_parser.get_structure("sample", file_path)	
	#first model of the sample structure 
	query_model = query_structure[0]
	
	superimpositions_dict, superimposed_chains, best_RMSD, iteration_count = dict_superimposition(ref_structure, query_structure, RMSD_threshold, iteration_count)

	#No superimposed chains or RMSD higer the threshold -> call McComplex function.
	if superimposed_chains is False or best_RMSD > RMSD_threshold:
		files_list, iteration_count, n = update_loop(files_list, iteration_count, n)
		return McComplex(ref_structure=ref_structure,
			files_list = files_list, iteration = iteration_count, parsed_no_chains = n, \
				command_arguments = command_arguments)

	else: 
		for chains, superimp in superimpositions_dict:
			l.info(f"Superimposition of ref chain {chains[0]} with sample chain {chains[1]} with an RMSD of {superimp.rms}")
			if superimp.rms > RMSD_threshold:
				l.info("This superimposition of ref chain %s with sample chain %s has an RMSD bigger than the threshold, therefore it is skipped" % (chains[0],chains[1]))
				continue
			# Rot and transl. matrices to all the atoms of the query model
			superimp.apply(query_model.get_atoms())
			# Gets the chain that was not superimposed with the reference 
			# chain, the putative chain to add ##
			putative_chain  = [chain for chain in query_model.get_chains() \
				if chain.get_id() != chains[0]][0]
			
			# to check whether the putative chain in the complex. 
			present_chain = False
			query_atoms, query_molecule = get_key_atoms(putative_chain )
			l.info("Putative chain to add is %s" % putative_chain .id)


			## Loops through all the chains from the reference structure ##
			for chain in ref_structure[0].get_chains():					
				ref_atoms, ref_molecule = get_key_atoms(chain)		
				## Neighbor Search to look for clashes between the chain 
				# to add and the chains from the reference structure ##
				Neighbor = Bio.PDB.NeighborSearch(ref_atoms)			
				clashes = []
				#loops through the list of atoms of putative_chain 		
				for atom in query_atoms:								
					atoms_clashed = Neighbor.search(atom.coord,5)		
					if len(atoms_clashed) > 0:
						clashes.extend(atoms_clashed)
				if len(clashes) > clashes_threshold:
					present_chain = True
					l.info("The number of clashes between the chain to add %s and reference chain %s is %d, therefore the chain is skipped" % (putative_chain .id, chain.id,len(clashes)))
					break
				elif len(clashes) <= clashes_threshold:
					l.info("The number of clashes between the chain to add %s and reference chain %s is %d, it is under the threshold" % (putative_chain .id, chain.id,len(clashes)))
					continue


			# Rotated chain to add is not a chain already in the ref structure, 
			# so it adds it
			if present_chain is False:
				l.info("Chain %s superimposed with chain %s yields rotated chain %s which is not in the complex" %(chains[0],chains[1],putative_chain .id))
				chain_ids = [chain.id for chain in ref_structure[0].get_chains()]	
				#list of IDs of all chains in ref structure
				ID = id_maker(chain_ids, putative_chain .id)
				putative_chain .id = ID
				ref_structure[0].add(putative_chain )
				l.info("Added Chain %s" % ID)
				## Checks whether the user provided the iterations argument, 
				# then save each iteration of the current complex in a PDB file 
				if pdb_iterations:
					save_iterations(ref_structure, outdir)
				file = files_list.pop(0)
				files_list.append(file)
				iteration_count += 1
				n = 0
				return McComplex(ref_structure=ref_structure,
					files_list = files_list, iteration = iteration_count, parsed_no_chains = n, \
						command_arguments = command_arguments)

	files_list, iteration_count, n = update_loop(files_list, iteration_count, n)

	return McComplex(ref_structure = ref_structure, \
		files_list = files_list, iteration = iteration_count, parsed_no_chains = n, \
			command_arguments = command_arguments)



