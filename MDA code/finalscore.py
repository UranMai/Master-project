import os
import sys
import pandas as pd


def load_monomer_data(pdb_name):
    """
    Load monomer data for test structure (1btl.pdb, 1gvp.pdb)
    """
    data1 = pd.read_csv(pdb_name + '_scoresEnergetics', sep='\t', index_col='AA')
    data1Z = pd.read_csv(pdb_name + '_scoresEnergeticsZ', sep='\t', index_col='AA')
    data1Z = data1Z.fillna(value=0)

    data2 = pd.read_csv(pdb_name + '_scoresCentroid', sep='\t', index_col='AA')
    data2Z = pd.read_csv(pdb_name + '_scoresCentroidZ', sep='\t', index_col='AA')
    data2Z = data2Z.fillna(value=0)

    data2 = data2[data2.index.isin(data1.index)]
    data2Z = data2Z[data2Z.index.isin(data1Z.index)]

    data2Z.columns = data2Z.columns + 'CENTROIDSC'

    return data1, data1Z, data2, data2Z


def load_multimer_data(pdb_name):
    """
    Load multimer data. For test `2oob` structure consider ubiquitin ligase ligand as HETATM
    and run scripts for it with (multimer data) and without (monomer data) ligase ligand
    dir_path - points out the path to folder with multimer data
    for monomer data it's same folder as for function `load_monomer_data`
    """
    dir_path = '.'
    pdb_name = os.path.join(dir_path, pdb_name)
    data1Multimer = pd.read_csv(pdb_name + '_scoresEnergetics', sep='\t', index_col='AA')
    data1MultimerZ = pd.read_csv(pdb_name + '_scoresEnergeticsZ', sep='\t', index_col='AA')
    data1MultimerZ = data1MultimerZ.fillna(value=0)

    data2Multimer = pd.read_csv(pdb_name + '_scoresCentroid', sep='\t', index_col='AA')
    data2MultimerZ = pd.read_csv(pdb_name + '_scoresCentroidZ', sep='\t', index_col='AA')
    data2MultimerZ = data2MultimerZ.fillna(value=0)

    data2Multimer = data2Multimer[data2Multimer.index.isin(data1Multimer.index)]
    data2MultimerZ = data2MultimerZ[data2MultimerZ.index.isin(data1MultimerZ.index)]

    data1MultimerZ.columns = data1MultimerZ.columns + 'MULTIMER'
    data2MultimerZ.columns = data2MultimerZ.columns + 'MULTIMERCENTROIDSC'

    return data1Multimer, data1MultimerZ, data2Multimer, data2MultimerZ


if __name__ == "__main__":
    """
    Final Score of amino acid is result of :
    1. Second Order Intermodular Degree - two degrees of separation between residues from different higher-order structures, 
    as an average of classical secondary structure and Walktrap definitions.
    2. Node-Edge Betweenness - the summed frequency that a node’s edges were utilized as a shortest path between
    all pairs of nodes in the network
    3. Euclidean Distance from Centroid to Ligand - the distance b/w residue’s centroid and the center of mass of the protein’s ligand.
    If multimer exists, 1. and 3. calculated for only multimer conformation, 2. for monomer and multimer
    The path to multimer is considered in function `load_multimer_data`
    """
    pdb = sys.argv[1]
    pdb_name = os.path.splitext(pdb)[0]

    data1, data1Z, data2, data2Z = load_monomer_data(pdb_name)
    data1Multimer, data1MultimerZ, data2Multimer, data2MultimerZ = load_multimer_data(pdb_name)

    data1 = pd.merge(data2[['RSA']], data1, on='AA')
    dataZ = pd.merge(data1Z, data2Z, on='AA')

    # monomer acids and multimer acid can differ, so consider missing nodes (monomeric)
    acidsMonomer = dataZ.index.to_list()
    acidsMultimer = data1MultimerZ.index.to_list()
    allAcids = list(set(acidsMonomer + acidsMultimer))
    missingNodes = [acid for acid in allAcids if acid not in acidsMultimer]

    # Calculate SecondOrderIntermodularDegree as the average of the multimer only
    SecondOrderIntermodularDegree_ALL = (data1MultimerZ["SecondOrderIntermodularDegreeSTRIDE_sidechainMULTIMER"] +
                                         data1MultimerZ["SecondOrderIntermodularDegreeSTRIDEMULTIMER"] +
                                         data2MultimerZ['SecondOrderIntermodularDegreeSTRIDEMULTIMERCENTROIDSC'] +
                                         data2MultimerZ["SecondOrderIntermodularDegreeWALKTRAPMULTIMERCENTROIDSC"]) / 4

    SecondOrderIntermodularDegree_AVERAGE = SecondOrderIntermodularDegree_ALL.to_dict()
    for node in missingNodes:
        value = (data1Z[data1Z.index == node]["SecondOrderIntermodularDegreeSTRIDE_sidechain"] +
                 data1Z[data1Z.index == node]["SecondOrderIntermodularDegreeSTRIDE"] +
                 data2Z[data2Z.index == node]["SecondOrderIntermodularDegreeSTRIDECENTROIDSC"] +
                 data2Z[data2Z.index == node]["SecondOrderIntermodularDegreeWALKTRAPCENTROIDSC"]) / 4
        SecondOrderIntermodularDegree_AVERAGE[node] = value

    # NodeEdgeBetweenness is the max of the monomer or multimer
    NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER = (dataZ["NodeEdgeBetweennessSTRIDE"] + dataZ["NodeEdgeBetweennessSTRIDE_sidechain"]) / 2
    NodeEdgeBetweennessSTRIDE_AVERAGE_MULTIMER = (data1MultimerZ["NodeEdgeBetweennessSTRIDEMULTIMER"] +
                                                  data1MultimerZ["NodeEdgeBetweennessSTRIDE_sidechainMULTIMER"]) / 2

    NodeEdgeBetweennessSTRIDE_sidechain_ALL = NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER.append(NodeEdgeBetweennessSTRIDE_AVERAGE_MULTIMER)
    df = NodeEdgeBetweennessSTRIDE_sidechain_ALL
    NodeEdgeBetweennessSTRIDE_sidechain_MAX = {acid: df[df.index == acid].max() for acid in allAcids}

    # Ligand is the min of the multimer only
    LigandMULTIMERCENTROIDSC_ALL = data2MultimerZ["LigandMULTIMERCENTROIDSC"]
    df = LigandMULTIMERCENTROIDSC_ALL
    LigandMULTIMERCENTROIDSC_MIN = {acid: df[df.index == acid].min() for acid in acidsMultimer}

    for node in missingNodes:
        value = dataZ['LigandCENTROIDSC']
        LigandMULTIMERCENTROIDSC_MIN[node] = value

    # Calculations of FINAL SCORE
    final_score = {}
    for node in allAcids:
        score = (SecondOrderIntermodularDegree_AVERAGE[node] + NodeEdgeBetweennessSTRIDE_sidechain_MAX[node] - LigandMULTIMERCENTROIDSC_MIN[node])
        final_score[node] = score

    pd.DataFrame(final_score.items()).sort_values(by=0).to_csv(pdb_name + '_FinalSum', sep='\t', index=False, header=False)
