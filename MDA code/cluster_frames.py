import os
import sys
import pandas as pd
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

def cluster_frames(u, n_clusters, method='k-means'):
    '''
    @description
        Calculate distance (RMSD) matrix between al the frames of the Universe object and cluster by KMeans method based on RMSD matrix. 
        Find cluster centroid and weights and return their dictionaries
    @input
        universe - MDAnalysis Universe object
        n_clusters - select number of clusters 
        method - KMeans method and other methods        
    '''
    total_frames = u.trajectory.n_frames
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
        # weights[cluster.id] = cluster.size/total_frames
        weights[cluster.centroid] = cluster.size/total_frames
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
    pdb = sys.argv[1] # 1BTL.pdb
    tpr = sys.argv[2] # topology file
    xtc = sys.argv[3] # trajectory file
    # maybe must be third arg, PDB file
    os.system(f"stride {pdb} > {pdb.replace('pdb', 'phi')}") # stride must be in usr/local/bin
    PHI = pdb.replace('pdb', 'phi')

    if os.path.exists(pdb) and os.path.exists(PHI):
        phifile = open(PHI, 'r').readlines()
        # print(phifile)
        phi_data = [line for line in phifile if 'ASG' in line] 
        print('PDB and PHI files are ready')
    else:
        print('Point appropriate PDB and PHI files')
        exit()

    prepare_secondary_structure_file(pdb, phi_data)
    # print(phi_data)
    prepare_rsa_file(pdb, phi_data)

    # file1 = 'run_1.tpr' 
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
        centroids, weights = cluster_frames(u, n_clusters=10, method='k-means')
        # return weights dict of cluster centoids and weights
        for ts in u.trajectory:
            if ts.frame in centroids.values():
                frames_dict[ts.frame] = u.atoms.positions
    elif select == 'k_frame':
        k_frames = create_Kth_frames(500, total_frames)
        for ts in u.trajectory:
            if ts.frame in k_frames:
                frames_dict[ts.frame] = u.atoms.positions
    elif select == 'all':
        frames_dict[ts.frame] = u.atoms.positions
    
    # print(weights)
    print(frames_dict.keys())
    # these files now changes tpr file
    # prepare_secondary_structure_file(file1, phi_data)
    # prepare_rsa_file(file1, phi_data)
    dataframes = []
    for ts in frames_dict.keys(): # iterate over most valuable time steps
        file1 = f'../Test_data/1btl_{ts}.pdb'
        # file1 = f'1btl_{ts}.pdb'
        print(f'Process  {ts}th frame')
        u = mda.Universe(tpr, frames_dict[ts])
        main(file1, u)
        print(f'End process of  {ts}th frame')

        cmd = f"python energetics_np.py {pdb} {file1}"
        os.system(cmd)

        cmd = f"python centroid_np.py {pdb} {file1}"
        os.system(cmd)

        # print(file1)
        cmd = f"python finalscore.py {file1}"
        print(cmd)
        os.system(cmd)

        tmp_df = pd.read_csv(file1.replace('.pdb', '_FinalSum'), sep='\t', header=None, names=['acid', f'score{ts}'])
        if select == 'centroid':
            tmp_df[f'score{ts}'] = weights[ts]*tmp_df[f'score{ts}']
        dataframes.append(tmp_df)
        # tmp_df = pd.merge(tmp_df, tmp_df, on='acid')
        # print(tmp_df)
        # break
    from functools import reduce
    data = reduce(lambda left,right: pd.merge(left,right,on='acid'), dataframes)
    data['acid'] = data['acid'].apply(lambda x: ':'.join(x.split(':')[:2]))
    data['mean'] = data.mean(axis=1)

    # Read data of mutational tolerance on 1btl
    mt_data = pd.read_csv('../Test_data/1btl_MT')
    out_data = pd.merge(mt_data, data, on='acid')
    print(out_data)
    # print(data.corr(method='pearson'))
    # print(data['acid'].to_list())
    print('Correlation: ', out_data['mt_score'].corr(out_data['mean']))
    # create_net_files(ts, frames_dict) # call function to create file with interactions, then run network analysis
    # cmd = f"sed -i 's/SYSTEM/-/g' adk_oplsaa_{ts}_net"
    # os.system(cmd)