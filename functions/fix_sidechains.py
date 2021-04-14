import subprocess
import sys
import nglview, ipywidgets
from biobb_model.model.fix_side_chain import FixSideChain


# Fix sidechain
prop = { 'use_modeller': True }
FixSideChain(input_pdb_path='macrocomplex.pdb', 
    output_pdb_path='fixed_macrocomplex.pdb', 
    properties = prop).launch()


