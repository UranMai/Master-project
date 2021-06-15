import os
import sys
import pandas as pd
import logging
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import MDAnalysis.analysis.encore as encore
from MDAnalysis.analysis.encore.clustering import ClusteringMethod as clm
from MDAnalysis.analysis.encore.clustering.ClusteringMethod import encode_centroid_info
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from pdb_class import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
from MDAnalysis.tests.datafiles import PDB, XTC
import matplotlib.pyplot as plt


'''Notes
Clustering methods
    https://towardsdatascience.com/cheat-sheet-to-implementing-7-methods-for-selecting-optimal-number-of-clusters-in-python-898241e1d6ad
ELBOW method
    Using longest distance to points on plot, angles not good method

'''

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("output4.log", mode='w'),
        logging.StreamHandler()
    ]
)
logging.getLogger("MDAnalysis").setLevel(logging.WARNING)

###Centering###
# for transformations can be used custom atom group
# transforms = (trans.unwrap(u.atoms), # u.atoms
#               trans.center_in_box(prot, wrap=True), 
#               trans.wrap(u.atoms))
# u.trajectory.add_transformations(*transforms)


###Clustering###
def create_rmsd_matrix(u):
    """
    Calculate rmsd_matirx as TriangularMatrix class
    And save as numpy binary file, .npz
    IF this binary file exists, create Triangular class and LOAD file

    Save rmsd_matrix as png image
    # need to point the pathDir
    """
    if not os.path.exists('rmsd_matrix.npz'):
        rmsd_matrix = encore.get_distance_matrix(ensemble=u, select='name CA', save_matrix='rmsd_matrix.npz')
        # rmsd_matrix = rmsd_matrix.as_array()
        # np.savez('matrix.npz', rmsd=rmsd_matrix)
        logging.info(f'Created RMSD matrix with shape {rmsd_matrix.as_array().shape} and save as npz file')
    else:
        rmsd_matrix = encore.utils.TriangularMatrix(u.trajectory.n_frames)
        rmsd_matrix.loadz('rmsd_matrix.npz')
        logging.info(f'Load RMSD matrix saved as numpy binary file')
    # plt.imshow(rmsd_matrix.as_array(), cmap='viridis')
    # plt.xlabel('Frame')
    # plt.ylabel('Frame')
    # plt.colorbar(label=r'RMSD ($\AA$)')
    # plt.savefig('RMSD.png', format='png')
    return rmsd_matrix

def segments_gain(p1, v, p2):
    # From https://datascience.stackexchange.com/questions/6508/k-means-incoherent-behaviour-choosing-k-with-elbow-method-bic-variance-explain
    # calculate angle b/w 3 points
    vp1 = np.linalg.norm(p1 - v)
    vp2 = np.linalg.norm(p2 - v)
    p1p2 = np.linalg.norm(p1 - p2)
    return np.arccos((vp1 ** 2 + vp2 ** 2 - p1p2 ** 2) / (2 * vp1 * vp2)) / np.pi

def optimal_number_of_clusters(matrix):
    '''
    https://jtemporal.com/kmeans-and-elbow-method/
    https://www.linkedin.com/pulse/finding-optimal-number-clusters-k-means-through-elbow-asanka-perera/
    Compute cluster sum of squares and distance b/w each point and line
    '''
    matrix = matrix.as_array()
    wcss = [] # within cluster sum of squares
    distortions = []

    K = range(2, 15) 
    for k in K:
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(matrix)
        wcss.append(kmeans.inertia_)
        distortions.append(sum(np.min(cdist(matrix, kmeans.cluster_centers_, 'euclidean'), axis=1)) / matrix.shape[0])

    x1, y1 = 2, wcss[0]
    x2, y2 = 20, wcss[len(wcss)-1]
    distances = []
    for i in range(len(wcss)):
        x0 = i+2
        y0 = wcss[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        distances.append(numerator/denominator)
    return distances.index(max(distances)) + 2


# def define_n_clusters(matrix):
#     '''Use method to define the smallest angle'''
#     # input is matrix as MDAnalysis.TriangularMatrix, so convert to np.array
#     matrix = matrix.as_array()
#     distortions = []

#     K = range(2, 15) # number of clusters. ELBOW method
#     for k in K:
#         kmeans = KMeans(n_clusters=k)
#         kmeans.fit(matrix)
#         distortions.append(sum(np.min(cdist(matrix, kmeans.cluster_centers_, 'euclidean'), axis=1)) / matrix.shape[0])

#     # Normalize the data
#     criterion = np.array(distortions)
#     criterion = (criterion - criterion.min()) / (criterion.max() - criterion.min())

#     # Compute the angles 
#     seg_gains = np.array([0, ] + [segments_gain(*
#                                                 [np.array([K[j], criterion[j]]) for j in range(i - 1, i + 2)]
#                                                 ) for i in range(len(K) - 2)] + [np.nan, ])
#     # Get the first index satisfying the threshold
#     seg_threshold = 0.99  # Set this to your desired target

#     kIdx = np.argmax(seg_gains > seg_threshold) # number of clusters
#     logging.info(f'Define N cluster in KMeans using Elbow method, it equals to {kIdx}')
#     return kIdx

# def Elbow_method(matrix):  
#     from sklearn.cluster import KMeans
#     from sklearn import metrics
#     from scipy.spatial.distance import cdist
#     import numpy as np
#     import matplotlib.pyplot as plt
#     matrix = matrix.as_array()
#     distortions = []
#     inertias = []
#     mapping1 = {}
#     mapping2 = {}
#     K = range(2, 25)
#     for k in K:
#         kmeans = KMeans(n_clusters=k)
#         kmeans.fit(matrix)
#         distortions.append(sum(np.min(cdist(matrix, kmeans.cluster_centers_, 'euclidean'), axis=1)) / matrix.shape[0])
#         inertias.append(kmeans.inertia_)
     
#         mapping1[k] = sum(np.min(cdist(matrix, kmeans.cluster_centers_, 'euclidean'), axis=1)) / matrix.shape[0]
#         mapping2[k] = kmeans.inertia_
#     plt.plot(K, distortions, 'bx-')
#     # plt.xlabel('Values of K')
#     # plt.ylabel('Distortion')
#     # plt.title('The Elbow Method using Distortion')
#     # plt.show()
#     plt.savefig('ELBOW.png')

def heirarchy_clustering(u):
    '''
    Use linkage clustering and plot dendrogram to choose number of clusters
    '''
    from scipy.cluster.hierarchy import dendrogram, linkage
    rmsd_matrix = create_rmsd_matrix(u)
    # plt.figure(figsize=(100, 100))
    plt.figure()
    linkage = linkage(rmsd_matrix.as_array(), method='single')
    dendrogram(linkage)
    plt.savefig('dendrogram.png', format='png', bbox_inches='tight')


class MeanShift(clm.ClusteringMethod): # create mda.ClusteringMethod from sklearn function MeanShift
    #https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html
    def __init__(self):
        from sklearn.cluster import MeanShift
        self.meanshift = MeanShift()

    def __call__(self, distance_matrix):
        clusters = self.meanshift.fit_predict(distance_matrix.as_array())
        cluster_representatives = np.unique(clusters, return_index=True)[1]
        clusters = encode_centroid_info(clusters,
                                            cluster_representatives)
        return clusters



def cluster_frames(u, method='meanshift'): # k-means, affinity, 
    '''out cluster stats
    @description
        Calculate distance (RMSD) matrix between al the frames of the Universe object and cluster by KMeans method based on RMSD matrix. 
        Find cluster centroid and weights and return their dictionaries
    @input
        universe - MDAnalysis Universe object
        n_clusters - select number of clusters 
        method - KMeans method and other methods        
    '''
    total_frames = u.trajectory.n_frames
    rmsd_matrix = create_rmsd_matrix(u)

    # n_clusters = define_n_clusters(rmsd_matrix)
    n_clusters = optimal_number_of_clusters(rmsd_matrix)
    # Define number of clusters 
    if method == 'k-means':
        km1 = clm.KMeans(n_clusters=n_clusters, algorithm="auto")
    elif method == 'affinity': # https://scikit-learn.org/stable/modules/clustering.html#affinity-propagation
        km1 = clm.AffinityPropagationNative()
    elif method == 'meanshift':
        km1 = MeanShift()


    # https://docs.mdanalysis.org/1.0.0/documentation_pages/analysis/encore/confdistmatrix.html
    # calculates the conformational distance (RMSD) matrix. 
    # The distance matrix is calculated between all the frames of all the Universe objects given as input.
    
    cluster_collection = encore.cluster(ensembles=u, method=km1, distance_matrix=rmsd_matrix)
    centroids = {}
    weights = {}
    for cluster in cluster_collection.clusters:
        centroids[cluster.id] = cluster.centroid
        weights[cluster.centroid] = cluster.size/total_frames
        logging.info(f'Cluster {cluster.id} has {cluster.size} elements. Frame {cluster.centroid} is cluster centroid')
    return centroids, weights

def create_Kth_frames(K, total_frames):
    # K - process each Kth frame from trajectory
    k_frames = list(range(total_frames))[4::K]
    logging.info(f'Process each K={K}th frame of trajectory: {k_frames}. You can change this K value')
    return k_frames

# def create_dict_frames(u, method='centroid'):
#     '''
#     For selected time frames create dict of atoms positions
#     Write atoms positions for considered frames 
#         'centroid' - select frames from clustering function
#         'all' - select all frames
#         'k_frame' - select each k-th frame from trajectory
#     '''    
#     logging.info(f'Choose {method} method to select frames and its coordinates')

#     frames_dict = {}

#     if method == 'centroid':
#         centroids, weights = cluster_frames(u, method='k-means')
#         # return weights dict of cluster centoids and weights
#         for ts in u.trajectory:
#             if ts.frame in centroids.values():
#                 frames_dict[ts.frame] = u.atoms.positions
#     elif method == 'k_frame':
#         k_frames = create_Kth_frames(5, total_frames)
#         for ts in u.trajectory:
#             if ts.frame in k_frames:
#                 frames_dict[ts.frame] = u.atoms.positions
#     elif method == 'all':
#         logging.info('Process all frames of trajectory')
#         for ts in u.trajectory:
#             frames_dict[ts.frame] = u.atoms.positions
#     logging.info('Created dictionary of frames atoms positions')
#     return frames_dict

if __name__ == '__main__':  
    tpr = sys.argv[1] # topology file
    xtc = sys.argv[2] # trajectory file
    logging.info(f'Start trajectory frames clustering on {tpr} and {xtc} files')
    pathDir = os.path.dirname(tpr)


    u = mda.Universe(tpr, xtc)
    prot = u.select_atoms('protein')
    total_frames = u.trajectory.n_frames # total number of frames 
    logging.info(f'Number of frames (time steps) - {total_frames}')
    


    # frames_dict = create_dict_frames(u, method='centroid')
    method = 'centroid'
    logging.info(f'Choose {method} method to select frames and its coordinates')
    frames_dict = {}
    if method == 'centroid':
        centroids, weights = cluster_frames(u, method='k-means')
        # return weights dict of cluster centoids and weights
        for ts in u.trajectory:
            if ts.frame in centroids.values():
                frames_dict[ts.frame] = u.atoms.positions
    elif method == 'k_frame':
        k_frames = create_Kth_frames(5, total_frames)
        for ts in u.trajectory:
            if ts.frame in k_frames:
                frames_dict[ts.frame] = u.atoms.positions
    elif method == 'all':
        logging.info('Process all frames of trajectory')
        for ts in u.trajectory:
            frames_dict[ts.frame] = u.atoms.positions

    dataframes = []
    for ts in frames_dict.keys():
        tmp_u = mda.Universe(tpr, frames_dict[ts])
        allatoms = tmp_u.select_atoms('all')

        frame = os.path.join(pathDir, f'frame_{ts}.pdb')
        allatoms.write(frame)
        logging.info(f'Traverse interested frames and write PDB files:  {frame} ')    

        # Read each frame and find for it all interactions, bonds, ligands
        u = mda.Universe(frame)
        
        phi_fname = frame.replace('pdb', 'phi')
        os.system(f"stride {frame} > {phi_fname}")
        logging.info(f"Run stride {frame} > {phi_fname}")


        main(PDB=frame, u=u)
        logging.info('Created NET file and find ligands. Lets compute centroid and energetics')

        logging.info('Run energetics, centroid scripts. Input files are initial pdb file and frame pdb file')
        cmd = f"python energetics_np.py {frame}"
        os.system(cmd)
        logging.info(cmd)

        cmd = f"python centroid_np.py {frame}"
        os.system(cmd)
        logging.info(cmd)

        logging.info('Run calculation of final score')
        cmd = f"python finalscore.py {frame}"
        os.system(cmd)
        logging.info(cmd + '\n')

        tmp_df = pd.read_csv(frame.replace('.pdb', '_FinalSum'), sep='\t', header=None, names=['acid', f'score{ts}'])
        
        if method == 'centroid':
            tmp_df[f'score{ts}'] = weights[ts]*tmp_df[f'score{ts}']
            print(tmp_df)
        dataframes.append(tmp_df)

    from functools import reduce
    data = reduce(lambda left,right: pd.merge(left,right,on='acid'), dataframes)
    data['acid'] = data['acid'].apply(lambda x: ':'.join(x.split(':')[:2]))
    # data['mean'] = data.mean(axis=1)
    data.to_csv(pathDir+'/k_means_scores_weighted', sep='\t', index=False)
    # print(pathDir + '/Score_data4')

    # # Read data of mutational tolerance on 1btl
    # data = pd.read_csv(pathDir+'/Score_data4', sep='\t')
    # mt_data = pd.read_csv('../Test_data/1btl_MT')
    # out_data = pd.merge(mt_data, data, on='acid')
    # print(out_data)
    # # # print(data.corr(method='pearson'))
    # # # print(data['acid'].to_list())
    # print('Correlation: ', out_data['mt_score'].corr(out_data['mean']))