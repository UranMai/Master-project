import os
import sys
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import MDAnalysis.analysis.encore as encore
from MDAnalysis.analysis.encore.clustering import ClusteringMethod as clm

from pdb_class import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
# pip install --upgrade MDAnalysisTests
from MDAnalysis.tests.datafiles import PDB, XTC

# Input data, same as PDB, XTC from tests.datafiles
# https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/data/adk_oplsaa.pdb
# https://github.com/MDAnalysis/mdanalysis/blob/develop/testsuite/MDAnalysisTests/data/adk_oplsaa.xtc
# PDB = '../Test_data/adk_oplsaa.pdb'  
# XTC = '../Test_data/adk_oplsaa.xtc'
# PHI = 'adk_oplsaa.phi'

file1 = '../Test_data/adk_oplsaa.pdb' 
# PHI = PDB.replace('.pdb', '.phi')
PHI = '../Test_data/adk_oplsaa.phi'
phifile = open(PHI, 'r').readlines()
phi_data = [line for line in phifile if 'ASG' in line]  

u = mda.Universe(PDB, XTC)
prot = u.select_atoms('protein')
total_frames = u.trajectory.n_frames # total number of frames 
print(f'Number of frames (time steps) - {total_frames}')

###############
###Centering###
###############
# https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/
# for transformations can be used custom atom group
transforms = (trans.unwrap(u.atoms), # u.atoms
              trans.center_in_box(prot, wrap=True), 
              trans.wrap(u.atoms))

u.trajectory.add_transformations(*transforms)


# ################
# ###Clustering###
# ################
# # https://docs.mdanalysis.org/stable/documentation_pages/analysis/encore/clustering.html
# https://scikit-learn.org/stable/modules/classes.html#module-sklearn.cluster
# methods = [clm.AffinityPropagation(),   
#            clm.AffinityPropagationNative(),
#            clm.DBSCAN()]

def cluster_frames(universe, n_clusters, method='k-means'):
    '''
    @description
        Calculate distance (RMSD) matrix between al the frames of the Universe object and cluster by KMeans method based on RMSD matrix. 
        Find cluster centroid and weights and return their dictionaries
    @input
        universe - MDAnalysis Universe object
        n_clusters - select number of clusters 
        method - KMeans method and other methods        
    '''
    # Define number of clusters 
    if method == 'k-means':
        km1 = clm.KMeans(n_clusters=n_clusters, algorithm="auto")

    # https://docs.mdanalysis.org/1.0.0/documentation_pages/analysis/encore/confdistmatrix.html
    # calculates the conformational distance (RMSD) matrix. 
    # The distance matrix is calculated between all the frames of all the Universe objects given as input.
    rmsd_matrix = encore.get_distance_matrix(ensemble=u)

    cluster_collection = encore.cluster(ensembles=u, method=km1, distance_matrix=rmsd_matrix)
    centroids = {}
    weights = {}
    for cluster in cluster_collection.clusters:
        centroids[cluster.id] = cluster.centroid
        weights[cluster.id] = cluster.size/total_frames
        print(f'Cluster {cluster.id} has {cluster.size} elements. Frame {cluster.centroid} is cluster centroid')
    return centroids, weights

def create_Kth_frames(k, total_frames):
    k_frames = list(range(total_frames))[::k]
    return k_frames


'''
For selected time frames create dict of atoms positions
Write atoms positions for considered frames 
    'centroid' - select frames from clustering function
    'all' - select all frames
    'k_frame' - select each k-th frame from trajectory
'''
select = 'k_frame'
frames_dict = {}
if select == 'centroid':
    centroids, weights = cluster_frames(universe=u, n_clusters=2, method='k-means')
    for ts in u.trajectory:
        if ts.frame in centroids.values():
            frames_dict[ts.frame] = u.atoms.positions
elif select == 'k_frame':
    k_frames = create_Kth_frames(5, total_frames)
    for ts in u.trajectory:
        if ts.frame in k_frames:
            frames_dict[ts.frame] = u.atoms.positions
elif select == 'all':
    frames_dict[ts.frame] = u.atoms.positions


print(frames_dict.keys())

prepare_secondary_structure_file(file1, phi_data)
prepare_rsa_file(file1, phi_data)

def create_net_files(ts, frames_dict):
    '''
    @description
        Create final data with all interactions on each time step
        -Define Universe and other common variables for each time step
        -
    @input
        ts - time step of MD trajectory
        frames_dict - data with atoms positions of MD trajectory
    '''
    test = '../Test_data/adk_oplsaa.pdb' 
    file1 = test.replace('.pdb', f'_{ts}.pdb')
    print(file1)
    u = mda.Universe(PDB, frames_dict[ts])
    protein = u.select_atoms('protein') 

    hoh = u.select_atoms('resname HOH')
    metalls = u.select_atoms('resname {}'.format(' '.join(list(prs.metals))))
    # ligands = allatoms - protein - hoh_atoms - metalls
    ligands = u.select_atoms('not protein and not resname HOH') - metalls

    prot = protein.atoms #define for each atom of protein
    prot_resnames, prot_resids, prot_segids, prot_atoms = prot.resnames, prot.resids, prot.segids, prot.names

    allatoms = [(i+':'+str(j)+':'+k+':'+l) for i, j, k, l in zip(prot_resnames, prot_resids, prot_segids, prot_atoms)] # HIS:26:A:N

    residues = [(res+':'+str(rid)+':'+sid) for res, rid, sid in zip(prot_resnames, prot_resids, prot_segids)] # HIS:26:A
    # coords = {res+':'+atom : pos for res, atom, pos in zip(residues, prot_atoms, prot.positions)} # {HIS:26:A:N : array([2.61,1.454,10.018]}
    coords = {atom : pos for atom, pos in zip(allatoms, prot.positions)}
    chain = {} #like {'HIS:26:A:N': 'MC'}
    for res, atom in zip(residues, prot_atoms):
        if atom in ["N","O","C","CA","HA2","HA3"]:
            if res[0:3] == 'GLY' and atom in ["O","CA","HA2","HA3","N","NH"]:
                chain[res+':'+atom] = 'SC'
            else:
                chain[res+':'+atom] = 'MC'
        else:
            chain[res+':'+atom] = 'SC'

    acids_class = [AminoAcid(res) for res in prot.residues]
    
    saltBridges = []

    # Writing files for PIPI,PICation, SB, Disulfide bonds
    with open(file1.replace('.pdb', '_bonds'), 'w') as net:
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
                # bonds = find_salt_bridges(acids_class[i], acids_class[j])
                # if bonds: 
                #     net.write(bonds[0]) #why we write only one connection b/w acids
                bonds = find_disulfide_bonds(acids_class[i], acids_class[j])
                if bonds: 
                    for bond in bonds:
                        net.write(bond)

    # Find hydrogen bonds
    # find_hydrogen_bonds(file1, u, prot_segids, coords, chain)

    # Find vanderWaals bonds
    find_vdw_bonds(file1, prot, prot_resnames, prot_resids, prot_segids, prot_atoms, coords, chain)

    find_metal_bonds(file1, metalls, acids_class)

    # LIGANDS
    ligand_names = [i+':'+str(j)+':'+k for i, j, k in zip(ligands.residues.resnames, ligands.residues.resids, ligands.residues.segids)]
    ligand_centroids = dict(zip(ligand_names, ligands.center_of_geometry(compound='group')))

    # CENTROIDS
    # Calc weighted center of residues, where centroid_data.masses are weights
    centroid_data = u.select_atoms('protein and not resname DG DC DT DA and not backbone or (resname GLY and not name N C O)') 
    center_coords = centroid_data.center(centroid_data.masses, compound='group')
    centroid = centroid_data.residues
    centroid_names = [i+':'+str(j)+':'+k for i, j, k in zip(centroid.resnames, centroid.resids, centroid.segids)]
    centroid_coords = dict(zip(centroid_names, center_coords))
    
    # find_ligand_atom_bonds(file1, allatoms, ligand_centroids, chain, coords)
    # find_centroid_bonds(centroid_coords)
    # find_centroid_ligand_bonds(file1, centroid_coords, ligand_centroids)
   

# for ts in frames_dict.keys(): # iterate over most valuable time steps
    # create_net_files(ts, frames_dict) # call function to create file with interactions, then run network analysis
    # cmd = f"sed -i 's/SYSTEM/-/g' adk_oplsaa_{ts}_net"
    # os.system(cmd)

