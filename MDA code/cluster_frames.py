import os
import sys

import MDAnalysis as mda
import MDAnalysis.transformations as trans
from pdb_class import *

import MDAnalysis.analysis.encore as encore
from MDAnalysis.analysis.encore.clustering import ClusteringMethod as clm

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
# pip install --upgrade MDAnalysisTests
## from MDAnalysis.tests.datafiles import PDB, XTC


# Input data, same as PDB, XTC from tests.datafiles
# https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/data/adk_oplsaa.pdb
# https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/data/adk_oplsaa.xtc
PDB = 'adk_oplsaa.pdb' # 
XTC = 'adk_oplsaa.xtc'
PHI = 'adk_oplsaa.phi'
# phifile = file1.replace('.pdb', '.phi')

u = mda.Universe(PDB, XTC)
prot = u.select_atoms('protein')
print(f'Number of frames (time steps) - {u.trajectory.n_frames}')

###############
###Centering###
###############
# print(prot.positions[0:3]) # just for comparison

# https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/
# for transformations can be used custom atom group
transforms = (trans.unwrap(u.atoms), # u.atoms
              trans.center_in_box(prot, wrap=True), 
              trans.wrap(u.atoms))

u.trajectory.add_transformations(*transforms)
# print(u.select_atoms('protein').positions[0:3]) # just for comparison


# ################
# ###Clustering###
# ################
# # https://docs.mdanalysis.org/stable/documentation_pages/analysis/encore/clustering.html
km1 = clm.KMeans(2,  # no. clusters
                 init = 'k-means++',  
                 algorithm="auto")

# https://docs.mdanalysis.org/1.0.0/documentation_pages/analysis/encore/confdistmatrix.html
# calculates the conformational distance (RMSD) matrix. 
# The distance matrix is calculated between all the frames of all the Universe objects given as input.
rmsd_matrix = encore.get_distance_matrix(u)

cluster_collection = encore.cluster(u, method=km1, distance_matrix=rmsd_matrix)
centroids = []
for cluster in cluster_collection.clusters:
    centroids.append(cluster.centroid)
    print(f'Cluster {cluster.id} has {cluster.size} elements. Frame {cluster.centroid} is cluster centroid')


frames_dict = {}
for ts in u.trajectory:
#     frames_dict[ts.frame] = u.atoms.positions
#     break
# # print(frames_dict)
    if ts.frame in centroids:
        frames_dict[ts.frame] = u.atoms.positions



# if __name__ == '__main__':
prepare_secondary_structure_file(PDB, PHI)
prepare_rsa_file(PDB, PHI)

def create_net_files(ts, frames_dict):
    '''
    ts - time step of mda trajectories
    frames_dict - data with positions for time step

    for each ts we run whole processing creating net file
    '''
    file1 = f'adk_{ts}.pdb'
    u = mda.Universe(PDB, frames_dict[ts])
    u1 = u.select_atoms('protein') 

    u1_atoms = u1.atoms
    u1_resnames, u1_resids, u1_segids, u1_names = u1_atoms.resnames, u1_atoms.resids, u1_atoms.segids, u1_atoms.names


    residues = [(i+':'+str(j)+':'+k) for i, j, k in zip(u1_atoms.resnames, u1_atoms.resids, u1_segids)]
    coords = {res+':'+atom : pos for res, atom, pos in zip(residues, u1_names, u1.atoms.positions)}
    chain = {} 
    for res, atom in zip(residues, u1_names):
        if atom in ["N","O","C","CA","HA2","HA3"]:
            if res[0:3] == 'GLY' and atom in ["O","CA","HA2","HA3","N","NH"]:
                chain[res+':'+atom] = 'SC'
            else:
                chain[res+':'+atom] = 'MC'
        else:
            chain[res+':'+atom] = 'SC'

    acids_class = [AminoAcid(res) for res in u1.residues]

    
    # Writing files for PIPI,PICation, SB, Disulfide bonds
    net = open(file1.replace('.pdb', '_bonds'), 'w')
    for i in range(len(acids_class)):
        for j in range(i+1, len(acids_class)):
            bonds = find_pipi_bonds(acids_class[i], acids_class[j])
            if bonds: 
                for bond in bonds:
                    net.write(bond)
            bonds = find_pication_bonds(acids_class[i], acids_class[j])
            if bonds: 
                for bond in bonds:
                    net.write(bond)
            # saltBrdiges = []
            bonds = find_salt_bridges(acids_class[i], acids_class[j])
            if bonds: 
                for bond in bonds:
                    net.write(bond)
            bonds = find_disulfide_bonds(acids_class[i], acids_class[j])
            if bonds: 
                for bond in bonds:
                    net.write(bond)
    net.close()

    # Find hydrogen bonds
    h = hbonds(u, selection1='protein', selection2= 'protein', distance=2.5, distance_type='hydrogen', angle=120) 
    h.run()
    h.generate_table() 
    find_hydrogen_bonds(h, file1, u1_segids, coords, chain)

    # Find vanderWaals bonds
    vdw_coords = u1.atoms.positions
    pairs, distances = mda.lib.distances.capped_distance(vdw_coords, vdw_coords, max_cutoff=4.25, min_cutoff=0)
    find_vdw_bonds(pairs, file1, u1_resnames, u1_resids, u1_segids, u1_names, coords, chain)

    # this function uses for appropriate writing into file
    VandWaals_awk_replacement(file1.replace('.pdb','')+'_vdw')
    VandWaals_awk_replacement(file1.replace('.pdb','')+'_vdw2')

    # Find metal bonds
    select_metalls = 'resname {}'.format(' '.join(list(prs.metals)))
    metalls = u.select_atoms(select_metalls)
    find_metal_bonds(file1, metalls, acids_class)

    # Concatenate created files into one 
    pdb = file1.replace('.pdb', '')
    cmd = f'cat {pdb}_bonds {pdb}_hb {pdb}_vdw_noRepulse {pdb}_vdw_noRepulse2 {pdb}_metal > {pdb}_net'
    os.system(cmd) 

    # Delete unneccessary files
    cmd = f'rm {pdb}_bonds {pdb}_hb {pdb}_vdw_noRepulse {pdb}_vdw_noRepulse2 {pdb}_metal {pdb}_vdw {pdb}_vdw2'
    os.system(cmd)



for ts in frames_dict.keys(): # iterate over most valuable time steps
    create_net_files(ts, frames_dict) # call function to create file with interactions, then run network analysis
    cmd = "sed -i 's/SYSTEM/-/g' adk_{}_net".format(ts)
    os.system(cmd)


# to run need terminalAtoms file
# for ts in dict_frames.keys():
#     cmd = "python energetics_main.py adk_oplsaa.pdb adk_{}".format(ts)
#     print(cmd)
#     os.system(cmd)

# # need uniqueAtoms file
# for ts in dict_frames.keys():
#     cmd = "python centroid_main.py adk_oplsaa.pdb adk_{}".format(ts)
#     print(cmd)
#     os.system(cmd)

# for ts in dict_frames.keys():
#     cmd = "python finalscore.py adk_{}".format(ts)
#     print(cmd)
#     os.system(cmd)