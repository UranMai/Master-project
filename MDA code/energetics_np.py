import sys
import warnings
from igraph import Graph
from functools import reduce
from score_functions import *
from pandas.core.common import SettingWithCopyWarning

warnings.filterwarnings(action="ignore", category=DeprecationWarning)
warnings.filterwarnings(action="ignore", category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


def process_data(data_net):
    """
    @description
        Consider reverse interactions of edges and merge direct and reverse data (bi-directed interactions)
        Ex. for given (THR:271:A:C, MET:272:A:CB) consider reverse bond (MET:272:A:CB, THR:271:A:C) changing atoms places and chain types
        For merged data add info about whether acids atoms are terminal(1) or not terminal(0) using terminalAtoms dict
        In merged data for groups of 2 common acids, count total Energy (sumWeight), count chains and number of bonds b/w 2 acids
    @input
        data_net - np.array data of interactions ('*net' file after running pdb_class.py script) in format:
            [acid1, acid2, Chain types, Energy_value, Bond_type, atom1, atom2]
            [THR:271:A,	MET:272:A,	MCSC,	3.8923571582256375,	VDW, C,   CB]
    @output
        dict of np.arrays, where
            'out' - data of unique processed acids with not negative Energy, WeightSum
            'influence_all' - data with at least one terminal atom
    """
    reverse_data = data_net.copy()
    mask1 = np.where(reverse_data[:, 2] == 'MCSC')
    mask2 = np.where(reverse_data[:, 2] == 'SCMC')
    reverse_data[mask1, 2] = 'SCMC'
    reverse_data[mask2, 2] = 'MCSC'
    if data_net.shape[1] == 7:
        reverse_data[:, [0, 1, 5, 6]] = reverse_data[:, [1, 0, 6, 5]]
    elif data_net.shape[1] == 9:
        reverse_data[:, [0, 1, 5, 6, 7, 8]] = reverse_data[:, [1, 0, 6, 5, 8, 7]]

    data = np.vstack((data_net, reverse_data))
    data = np.unique(data, axis=0)

    if data.shape[1] == 7:  # Check if atoms are terminal or not
        codes1 = []
        codes2 = []
        for res, atom in data[:, [0, 5]]:
            codes1.append(1 if atom in terminalAtoms[res[:3]] else 0)
        for res, atom in data[:, [1, 6]]:
            codes2.append(1 if atom in terminalAtoms[res[:3]] else 0)
        data = np.insert(data, 7, codes1, axis=1)  # add terminal info of atoms to data
        data = np.insert(data, 8, codes2, axis=1)
        out_details = data
    else:
        out_details = data

    # Consider interactions with at least one terminal atom out of two and non negative energies
    influence_all = out_details[(out_details[:, 3].astype(float) > 0) & np.bitwise_or(out_details[:, 7] == '1', out_details[:, 8] == '1')]

    # Group data by common 2 acids, calculate WeightSUM, percentMC
    bonds = []
    processed = []
    for res1, res2 in data[:, [0, 1]]:
        if [res1, res2] not in processed:
            processed.append([res1, res2])
            grouped_data = data[(data[:, 0] == res1) & (data[:, 1] == res2)]
            edgeSum = sum(grouped_data[:, 3].astype('float'))
            MCMC_count, MCSC_count, SCSC_count, SCMC_count = 0, 0, 0, 0
            for tp in grouped_data[:, 2]:
                if tp == 'MCMC' or tp == 'MCSC':
                    MCMC_count += 1
                    MCSC_count += 1
                elif tp == 'SCMC':
                    SCMC_count += 1
                elif tp == 'SCSC':
                    SCSC_count += 1
            percentMC = (MCMC_count + MCSC_count) / (MCMC_count + MCSC_count + SCMC_count + SCSC_count)
            values = [grouped_data[0, 0], grouped_data[0, 1], 'mixed', edgeSum, percentMC, sum(grouped_data[:, 7].astype(int)),
                      sum(grouped_data[:, 8].astype(int)), len(data[data[:, 0] == grouped_data[0, 0]])]
            bonds.append(values)

    bonds1 = np.array(bonds)
    bonds1 = bonds1[bonds1[:, 3].astype(float) > 0]  # Consider only interactions with positive weightSum

    results = {'out': bonds1, 'influence_all': influence_all}
    return results


def removeRedundancy(data):
    """
    @description
        Select unique amino acids and for them define sub_data from data of interactions
        Then calculate the sum of energies for sub_data of these acids
        Here, bonds (HIS:1:A, PRO:2:A) and (PRO:2:A, HIS:1:A) process equally
    @input
        data - np.array of data of interactions, where first 2 columns - amino acids
    @output
        np.array in format [acid1, acid2, SumWeightEnergy]
    """
    pairs = np.sort(data[:, :2], axis=1)
    pairs = np.unique(pairs, axis=0)
    bonds = []
    for i in range(len(pairs)):
        sub_data = data[(data[:, 0] == pairs[i, 0]) & (data[:, 1] == pairs[i, 1])]
        bonds.append([pairs[i, 0], pairs[i, 1], sum(sub_data[:, 3].astype(float))])
    return np.array(bonds)


def make_directed_interactions(data):
    """
    @description
        Select unique amino acids and for them define sub_data from data of interactions
        Then calculate the sum of energies for sub_data of these acids
        Here, bonds (HIS:1:A, PRO:2:A) and (PRO:2:A, HIS:1:A) process differently
    @input
        data - np.array of data of interactions, where first 2 columns - amino acids
    @output
        np.array in format [acid1, acid2, SumWeightEnergy]
    """
    pairs = np.unique(data[:, :2], axis=0)
    bonds = []
    for i in range(len(pairs)):
        sub_data = data[((data[:, 0] == pairs[i, 0]) & (data[:, 1] == pairs[i, 1]))]
        bonds.append([sub_data[0, 0], sub_data[0, 1], sum(sub_data[:, 3].astype(float))])
    return np.array(bonds)


if __name__ == '__main__':
    pdb = sys.argv[1]  # input pdb structure, 1btl.pdb
    try:
        secStructure = np.loadtxt(pdb.replace('.pdb', '_secondaryStructure'), dtype=np.str)
        rsa = np.loadtxt(pdb.replace('.pdb', '.rsa'), dtype=np.str)
        data_all = np.loadtxt(pdb.replace('.pdb', '_net'), dtype=np.str)
    except FileNotFoundError:
        print('File not accessible')

    col0 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 0])))
    col1 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 1])))
    data_all[:, 0] = col0[:, 0]
    data_all[:, 1] = col1[:, 0]
    data_all = np.concatenate((data_all, col0[:, [1]]), axis=1)
    data_all = np.concatenate((data_all, col1[:, [1]]), axis=1)
    # data_all must be in format ['LYS:6:A' 'ASP:9:A' 'SCSC' '10' 'SB' 'NZ' 'OD2'] after processing
    # where columns 0,1-amino acids; 2-chain type; 3-energy value; 4-bond type; 5,6-atoms

    originalNodes = nodes = np.unique(data_all[:, 0:2])  # Define nodes (amino acids) that are interacting in net file

    basedata = process_data(data_all)
    basedata_noPP = removeRedundancy(basedata['out'])
    # Consider only side chain interactions
    # basedata_allsidechain = make_undirected_interactions(data_all[data_all[:,2]!='MCMC'])

    influencedata_all = basedata['influence_all']  # data with terminal atoms
    influencedata_collapse = process_data(influencedata_all)
    influencedata = removeRedundancy(influencedata_collapse['out'])

    influencenet = Graph.TupleList(edges=influencedata[:, :2], directed=False)
    influencenet_weights = 1 / influencedata[:, 2].astype(float)

    # Define directed graph network
    influencedata_directed = make_directed_interactions(influencedata_all[influencedata_all[:, 7] == '1'])
    influencenet_directed = Graph.TupleList(edges=influencedata_directed[:, [1, 0]], directed=True)
    influencenet_directed_weights = 1 / influencedata_directed[:, 2].astype(float)

    # Fill secStructure and rsa data with missing nodes and values
    if not any(np.isin(nodes, secStructure[:, 0])):
        for node in nodes:
            if not (node in secStructure[:, 0]):
                secStructure = np.append(secStructure, [[node, 'xxx', '360']], axis=0)
            if not (node in rsa[:, 0]):
                rsa = np.append(rsa, [[node, 0]], axis=0)

    # Define undirected network
    net = Graph.TupleList(edges=basedata_noPP[:, :2], directed=False)
    net_weights = 1 / basedata_noPP[:, 2].astype(float)

    """ WALKTRAP CLUSTERING
    Using basedata with all connections create GRAPH. Use Walktrap algo to define community
    clusters for nodes and write it into `net_community_vec` dict. And calculate network attributes """
    net_community = net.community_walktrap(net_weights, steps=4)
    cluster = net_community.as_clustering()
    net_community_vec = {name: member for name, member in zip(cluster.graph.vs['name'], cluster.membership)}
    nodes = net.vs['name']

    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_weights)
    node_edge_betweenness_wt = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_wt = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_wt, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_wt, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)
    secondOrder_node_intermodular_degree_wt, secondOrder_node_intermodular_degree_sidechain_wt = \
        node_intermodular_dgr_2order(node_intermodular_degree_wt, node_intermodular_degree_wt,
                                     node_modules, node_modules_sidechain, method='energetics')

    firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain = \
        first_second_order_degree(influencenet, influencenet_directed)

    """SECONDARY STRUCTURE CLUSTERING
    Use Stride created secStructure file to define community clusters for nodes.
    Based on type of secondary structure (structs) create `net_community_vec` dict. And calculate network attributes """
    factors = sorted(np.unique(secStructure[:, 1]))
    factors_dict = {factor: num for num, factor in enumerate(factors, 1)}
    net_community_vec = {node: factors_dict[f] for node, f in zip(secStructure[:, 0], secStructure[:, 1])}
    nodes = influencenet.vs['name']

    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_weights)
    node_edge_betweenness_stride = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_stride = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_stride, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_stride, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)
    secondOrder_node_intermodular_degree_stride, secondOrder_node_intermodular_degree_sidechain_stride = \
        node_intermodular_dgr_2order(node_intermodular_degree_stride, node_intermodular_degree_stride,
                                     node_modules, node_modules_sidechain, method='energetics')

    """ Create dataframes from dictionaries and merge them """
    firstOrderDegree_df = create_df(firstOrderDegree, nodes)
    firstOrderDegree_sidechain_df = create_df(firstOrderDegree_sidechain, nodes)

    secondOrderDegree_df = create_df(secondOrderDegree, nodes)
    secondOrderDegree_sidechain_df = create_df(secondOrderDegree_sidechain, nodes)

    node_edge_betweenness_stride_df = create_df(node_edge_betweenness_stride, nodes)
    node_edge_betweenness_sidechain_stride_df = create_df(node_edge_betweenness_sidechain_stride, nodes)

    node_intermodular_degree_stride_df = create_df(node_intermodular_degree_stride, nodes)
    node_intermodular_degree_sidechain_stride_df = create_df(node_intermodular_degree_sidechain_stride, nodes)

    secondOrder_node_intermodular_degree_sidechain_stride_df = create_df(secondOrder_node_intermodular_degree_sidechain_stride, nodes)
    secondOrder_node_intermodular_degree_stride_df = create_df(secondOrder_node_intermodular_degree_stride, nodes)

    node_edge_betweenness_wt_df = create_df(node_edge_betweenness_wt, nodes)
    node_edge_betweenness_sidechain_wt_df = create_df(node_edge_betweenness_sidechain_wt, nodes)

    node_intermodular_degree_wt_df = create_df(node_intermodular_degree_wt, nodes)
    node_intermodular_degree_sidechain_wt_df = create_df(node_intermodular_degree_sidechain_wt, nodes)

    secondOrder_node_intermodular_degree_wt_df = create_df(secondOrder_node_intermodular_degree_wt, nodes)
    secondOrder_node_intermodular_degree_sidechain_wt_df = create_df(secondOrder_node_intermodular_degree_sidechain_wt, nodes)

    out = [firstOrderDegree_df, firstOrderDegree_sidechain_df,
           secondOrderDegree_df, secondOrderDegree_sidechain_df,
           node_edge_betweenness_stride_df, node_edge_betweenness_sidechain_stride_df,
           node_intermodular_degree_stride_df, node_intermodular_degree_sidechain_stride_df,
           secondOrder_node_intermodular_degree_stride_df, secondOrder_node_intermodular_degree_sidechain_stride_df,
           node_edge_betweenness_wt_df, node_edge_betweenness_sidechain_wt_df,
           node_intermodular_degree_wt_df, node_intermodular_degree_sidechain_wt_df,
           secondOrder_node_intermodular_degree_wt_df, secondOrder_node_intermodular_degree_sidechain_wt_df]

    final_df = reduce(lambda left, right: pd.merge(left, right, on=0), out)
    colnames = ['AA', "Degree", "Degree_sidechain", "SecondOrderDegree", "SecondOrderDegree_sidechain",
                "NodeEdgeBetweennessSTRIDE", "NodeEdgeBetweennessSTRIDE_sidechain",
                "IntermodularDegreeSTRIDE", "IntermodularDegreeSTRIDE_sidechain",
                "SecondOrderIntermodularDegreeSTRIDE", "SecondOrderIntermodularDegreeSTRIDE_sidechain",
                "NodeEdgeBetweennessWALKTRAP", "NodeEdgeBetweennessWALKTRAP_sidechain",
                "IntermodularDegreeWALKTRAP", "IntermodularDegreeWALKTRAP_sidechain",
                "SecondOrderIntermodularDegreeWALKTRAP", "SecondOrderIntermodularDegreeWALKTRAP_sidechain"]
    final_df.columns = colnames

    # Add nodes that were not networked in influencenet graph
    zeroNodes = [node for node in originalNodes if node not in nodes]
    zeroMat = pd.DataFrame(0, index=np.arange(len(zeroNodes)), columns=final_df.columns[0:])
    zeroMat['AA'] = zeroNodes

    # Concat final dataset with zero dataframe, replace None with 0, and add 1 to everything to avoid zeros, except ligand
    out = pd.concat([final_df, zeroMat])
    out1 = out.fillna(value=0)
    out1 = out1[out1.columns[1:]] + 1
    out = pd.concat([out[['AA']], out1], axis=1)

    for col in out.columns[1:]:
        print(col, out[col].sum())

    outZ = Scaling(out[out.columns[1:]])
    outZ = pd.concat([out[['AA']], outZ], axis=1)

    out.to_csv(pdb.replace('.pdb', '_scoresEnergetics'), sep='\t', index=False)
    outZ.to_csv(pdb.replace('.pdb', '_scoresEnergeticsZ'), sep='\t', index=False)
