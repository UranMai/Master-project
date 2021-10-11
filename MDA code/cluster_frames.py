import pandas as pd
import logging
import MDAnalysis.analysis.encore as encore
from MDAnalysis.analysis.encore.clustering import ClusteringMethod as clm
from MDAnalysis.analysis import align
from pdb_class import *
from functools import reduce
import warnings
import argparse
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)


def create_rmsd_matrix(u, pdb_name):
    """
    Calculate rmsd_matrix as mda.TriangularMatrix class and save as numpy binary file, .npz
    If this binary file exists, create Triangular class and load it
    """
    if not os.path.exists(f'rmsd_matrix_{pdb_name}.npz'):
        rmsd_matrix = encore.get_distance_matrix(ensemble=u, select='name CA', save_matrix=f'rmsd_matrix_{pdb_name}.npz')
        # rmsd_matrix = rmsd_matrix.as_array()
        # np.savez('matrix.npz', rmsd=rmsd_matrix)
        logging.info(f'Created RMSD matrix with shape {rmsd_matrix.as_array().shape} and save as npz file')
    else:
        rmsd_matrix = encore.utils.TriangularMatrix(u.trajectory.n_frames)
        rmsd_matrix.loadz(f'rmsd_matrix_{pdb_name}.npz')
        logging.info(f'Load RMSD matrix saved as numpy binary file')
    return rmsd_matrix


def cluster_frames(u, pdb_name, method):
    """
    @description
        Calculate distance (RMSD) matrix between al the frames of the Universe object and cluster by KMeans method based on RMSD matrix.
        Find cluster centroid and weights and return their dictionaries
    @input
        u - MDAnalysis Universe object
        pdb - select appropriate pdb file name and its rmsd numpy binary file
        method - KMeans method or other methods: affinity, dbscan
    """
    total_frames = u.trajectory.n_frames
    rmsd_matrix = create_rmsd_matrix(u, pdb_name)

    if method == 'k-means':
        km1 = clm.KMeans(n_clusters=7, algorithm="auto", random_state=42)
        logging.info('Cluster frames using KMeans')
    elif method == 'affinity':
        km1 = clm.AffinityPropagationNative(damping=2, max_iter=300, convergence_iter=50, add_noise=True)
        logging.info('Cluster frames using AffinityPropagationNative')
    elif method == 'dbscan':
        # km1 = clm.DBSCAN(eps=1, min_samples=100, algorithm='auto', leaf_size=30) # 1bt;
        km1 = clm.DBSCAN(eps=1, min_samples=20, leaf_size=30)  # 1gvp, 2oob
        logging.info('Cluster frames using DBSCAN')
    else:
        return None

    # https://docs.mdanalysis.org/1.0.0/documentation_pages/analysis/encore/confdistmatrix.html
    cluster_collection = encore.cluster(ensembles=u, method=km1, distance_matrix=rmsd_matrix)
    centroids = {}
    weights = {}
    for cluster in cluster_collection.clusters:
        centroids[cluster.id] = cluster.centroid
        weights[cluster.centroid] = cluster.size / total_frames
        logging.info(f'Cluster {cluster.id} has {cluster.size} elements. Frame {cluster.centroid} is cluster centroid with weight {cluster.size / total_frames}')

    return centroids, weights


def create_Kth_frames(k, total_frames):
    """
    Define what frames to be proceed, select each k-th frame in the trajectory among total frames
    """
    frames = list(range(total_frames))[::k]
    logging.info(f'Process each {k}th frame of trajectory: {frames}. You can change this value')
    return frames


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("LOGFILE_new_xtc", mode='w'),
            logging.StreamHandler()
        ]
    )
    logging.getLogger("MDAnalysis").setLevel(logging.WARNING)

    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--topol", required=True, help='path to input topology file')
    ap.add_argument("-t", "--traj", required=True, help='path to input trajectory file centered by gromacs')
    ap.add_argument("-w", "--weighted", default=True, required=False, help='whether use weighted scores or not')
    ap.add_argument("-m", "--method", choices=['centroid', 'k_frame', 'all'], help='method to select frames for processing')
    ap.add_argument("-a", "--algo", required=False, help='if centroid method selected, need to choose clustering algo')
    args = vars(ap.parse_args())

    weighted = args['weighted']
    method = 'k_frame'  # args['method']
    algo = args['algo']
    tpr = args['topol']
    xtc = args['traj']

    pathDir = os.path.dirname(tpr)  # path/to/tpr
    pdb = os.path.basename(tpr)  # run_1_1btl.tpr
    pdb = pdb.replace('.tpr', '').replace('run_1_', '')  # run_1_1btl.tpr -> 1btl
    logging.info(f'Process {tpr} and {xtc} files, method to choose frames is {method, algo}')

    u = mda.Universe(tpr, xtc)
    aligner = align.AlignTraj(u, u, select='name CA', in_memory=True).run()  # align to the first frame
    total_frames = u.trajectory.n_frames
    logging.info(f'Number of frames (time steps) - {total_frames}')

    frames_dict = {}
    if method == 'centroid':
        centroids, weights = cluster_frames(u, pdb, method=args['algo'])
        for ts in u.trajectory:
            if ts.frame in centroids.values():
                frames_dict[ts.frame] = u.atoms.positions
    elif method == 'k_frame':
        k_frames = create_Kth_frames(50, total_frames)
        for ts in u.trajectory:
            if ts.frame in k_frames:
                frames_dict[ts.frame] = u.atoms.positions
    elif method == 'all':
        for ts in u.trajectory:
            frames_dict[ts.frame] = u.atoms.positions

    dataframes_weighted = []  # weighted dataframes_weighted for centroid method
    dataframes = []  # unweighted dataframes_weighted

    for ts in frames_dict.keys():  # iterate each frame, write pdb file and for it find
        tmp_u = mda.Universe(tpr, frames_dict[ts])
        allatoms = tmp_u.select_atoms('protein or (resname SOL and around 3 protein)')  # 'all' - here can select what atoms to be proceed
        frame = os.path.join(pathDir, f'frame_{ts}.pdb')
        allatoms.write(frame)
        logging.info(f'Traverse interested frames and write PDB files: {frame} ')

        logging.info('Find protein bonds. Run energetics, centroid scripts and calculate final score')
        os.system(f'python pdb_class.py {frame}')
        os.system(f"python energetics_np.py {frame}")
        os.system(f"python centroid_np.py {frame}")
        os.system(f"python finalscore.py {frame}")

        tmp_df = pd.read_csv(frame.replace('.pdb', '_FinalSum'), sep='\t', header=None, names=['acid', f'score{ts}'])
        dataframes.append(tmp_df)

        cmd = f"cd {pathDir}; \
                rm -rf *_net *_ligands *.rsa *_phi *_secondaryStructure2 *_secondaryStructure *_centroidNetLigand \
                rm -rf *_scoresEnergetics *_scoresEnergeticsZ *_scoresCentroid *_scoresCentroidZ *_centroidNetSC *.phi"
        os.system(cmd)

        # if weighted == 'True' and method == 'centroid':
        #     tmp_df1 = pd.read_csv(frame.replace('.pdb', '_FinalSum'), sep='\t', header=None, names=['acid', f'score{ts}'])
        #     # tmp_df1[f'score{ts}'] = weights[ts] * tmp_df1[f'score{ts}']
        #     dataframes_weighted.append(tmp_df1)
        break

    data = reduce(lambda left, right: pd.merge(left, right, on='acid'), dataframes)
    data['acid'] = data['acid'].apply(lambda x: ':'.join(x.split(':')[:2]))
    data.to_csv(pathDir + '/scores_frames', sep='\t', index=False)
    print(data)

    # data = reduce(lambda left,right: pd.merge(left,right,on='acid'), dataframes_weighted)
    # data['acid'] = data['acid'].apply(lambda x: ':'.join(x.split(':')[:2]))
    # data.to_csv(os.path.join(pathDir, f'Scores_{weighted}_{method}_{algo}'), sep='\t', index=False)
