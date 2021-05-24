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
from MDAnalysis.tests.datafiles import PDB, XTC
# Input data, same as PDB, XTC from tests.datafiles


###############
###Centering###
###############
# https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/
# for transformations can be used custom atom group
# transforms = (trans.unwrap(u.atoms), # u.atoms
#               trans.center_in_box(prot, wrap=True), 
#               trans.wrap(u.atoms))

# u.trajectory.add_transformations(*transforms)


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
    print('These frames: ', *k_frames, 'will be processed')
    return k_frames

# these files now changes tpr file
# prepare_secondary_structure_file(file1, phi_data)
# prepare_rsa_file(file1, phi_data)

if __name__ == '__main__':  
    tpr = sys.argv[1] # topology file
    xtc = sys.argv[2] # trajectory file
    # maybe must be third arg, PDB file

    file1 = 'run_1.tpr' 
    # PHI = PDB.replace('.pdb', '.phi')
    # PHI = '1BTL.phi'
    # phifile = open(PHI, 'r').readlines()
    # phi_data = [line for line in phifile if 'ASG' in line]  

    u = mda.Universe(tpr, xtc)
    prot = u.select_atoms('protein')
    total_frames = u.trajectory.n_frames # total number of frames 
    print(f'Number of frames (time steps) - {total_frames}')

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
        k_frames = create_Kth_frames(100, total_frames)
        for ts in u.trajectory:
            if ts.frame in k_frames:
                frames_dict[ts.frame] = u.atoms.positions
    elif select == 'all':
        frames_dict[ts.frame] = u.atoms.positions

    # these files now changes tpr file
    # prepare_secondary_structure_file(file1, phi_data)
    # prepare_rsa_file(file1, phi_data)

    for ts in frames_dict.keys(): # iterate over most valuable time steps
        file1 = f'1btl_{ts}.pdb'
        print(f'Process  {ts}th frame')
        u = mda.Universe(tpr, frames_dict[ts])
        main(file1, u)
        print(f'End process of  {ts}th frame')

    # create_net_files(ts, frames_dict) # call function to create file with interactions, then run network analysis
    # cmd = f"sed -i 's/SYSTEM/-/g' adk_oplsaa_{ts}_net"
    # os.system(cmd)