from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import re, glob
import sys
from pdbfixer import PDBFixer
'''
>>> python fix_pdb.py file.pdb
Fixer options: 
	1. find and add missing residues
	2. find and replace Nonstandard residues
	3. add missing hydrogens
	4. add missing atoms, 
'''

def fix_pdbfile(pdb_file, keepIds=False, add_H=False):
	'''
	@Fixes and return pdb file, adding missing atoms, residues, replace nonstandard residues

	pdb_file - file.pdb which will be processed
	keepIds - if False - remain residues' ids, if True - renumerate ids. It can be used for cases with wrong numeration
	add_H - add missing hydrogens at pH=7. Since Gromacs adds hydrogens, don't use it
	'''
	fixer = PDBFixer(filename=pdb_file)

	fixer.findMissingResidues()
	missing_res = fixer.missingResidues
	print(f'Missed residues {missing_res}')

	fixer.findNonstandardResidues()
	nonstandard_residues = fixer.nonstandardResidues
	fixer.replaceNonstandardResidues()
	print(f'Nonstandard res {nonstandard_residues}')

	if add_H:
		fixer.addMissingHydrogens(pH=7.0)

	fixer.findMissingAtoms()
	missingAtoms = fixer.missingAtoms
	missingTerminals = fixer.missingTerminals
	fixer.addMissingAtoms()
	print(f'Missing atoms {missingAtoms} and terminal atoms {missingTerminals}')

	PDBFile.writeFile(fixer.topology, fixer.positions, open(re.sub('.pdb', '_fixed.pdb', pdb_file), 'w'), keepIds=keepIds)

if __name__ == '__main__':
	pdb_file = sys.argv[1]
	fix_pdbfile(pdb_file)

# 	pdb_files = glob.glob('./*.pdb')
# 	pdb_files = ['1btl.pdb', '1fcc.pdb', '3fqm.pdb', '2oob.pdb', '2lzm.pdb', '2cg9.pdb',
# 				'1x75.pdb', '1nd4.pdb', '1gvp.pdb', '1gnx.pdb', '1dct.pdb', '3coq.pdb', '1rvx.pdb']
# 
# 	for pdb_file in pdb_files:
# 		print(pdb_file)
# 		fix_pdbfile(pdb_file)
