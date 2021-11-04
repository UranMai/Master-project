import sys
from igraph import Graph
import warnings
from functools import reduce
from pandas.core.common import SettingWithCopyWarning
from score_functions import *

warnings.filterwarnings(action="ignore", category=DeprecationWarning)
warnings.filterwarnings(action="ignore", category=UserWarning)
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


def add_terminal_atoms(data_net):
    """
    @description
        Create data of reverse interactions, Merge directed and reversed data (bi-directed interactions)
        Consider [HIS:26:A, PRO:27:A] and [ PRO:27:A, HIS:26:A] interactions
    @input
        data_net - np.array data of interactions ('centroidNetSC' file from pdb_class.py script)
    @output
        data - np.array of processed bi-directed data
    """
    reverse_data = data_net.copy()
    reverse_data[:, [0, 1]] = reverse_data[:, [1, 0]]
    data = np.vstack((data_net, reverse_data))
    data = np.unique(data, axis=0)
    return data


def make_undirected_interactions(data):
    """
    @description
        Consider undirected interactions: from input data sort acids by name and select unique pairs,
        Select `sub_data` grouped by unique pair of two acids and calculate sumWeight
    @input
        data - np.array data of interactions
    @output
        np.array in format [Acid1, Acid2, SumWeight]
    """
    pairs = np.sort(data[:, :2], axis=1)
    pairs = np.unique(pairs, axis=0)
    output = []
    for i in range(len(pairs)):
        sub_data = data[np.bitwise_and(data[:, 0] == pairs[i, 0], data[:, 1] == pairs[i, 1]) | np.bitwise_and(data[:, 1] == pairs[i, 0], data[:, 0] == pairs[i, 1])]
        output.append([pairs[i, 0], pairs[i, 1], 0.5 * sum(sub_data[:, 3].astype(float))])
    return np.array(output)


def exclude_neighbour_nodes(data, cutoff=1):
    """
    @description
        For each bond in data exclude neighbored acids which resid diff > cutoff,
        like [HIS:26:A, PRO:27:A] for cutoff = 1
    @input
        data - np.array data of interactions
    @output
        np.array of non-neighbored interactions
    """
    if data.shape[0] == 0:
        return data

    bonds = [line for line in data if np.abs(int(line[0].split(':')[1]) - int(line[1].split(':')[1])) > cutoff]
    return np.array(bonds)


if __name__ == '__main__':
    """
    - data_all ('centroidNetSC') must be in format [HIS:26:A, PRO:27:A, SCSC, 1, CENTROID, CENTROID, CENTROID, 4.70];
    where columns 0,1-amino acids; 2-chain type; 3-energy value; 7-distance b/w acids centroids
    - ligand_data ('ligands') in format [HIS:26:A:CB, SO4:291:A, 28.64];
    where 0-acid atom, 1-ligand, 2-distance b/w acid atom and ligand centroid
    - ligandCentroidFile ('centroidNetLigand') in format [HIS:26:A,	SO4:291:A, 30.527385116085433] 
    where 0-amino acid, 1-ligand, 2-distance b/w acid and ligand centroids
    """
    pdb = sys.argv[1]  # input pdb structure
    try:
        data_all = np.loadtxt(pdb.replace('.pdb', '_centroidNetSC'), dtype=np.str)
        secStructure = np.loadtxt(pdb.replace('.pdb', '_secondaryStructure'), dtype=np.str)
        rsa = np.loadtxt(pdb.replace('.pdb', '.rsa'), dtype=np.str)
        ligand_data = np.loadtxt(pdb.replace('.pdb', '_ligands'), dtype=np.str)
        ligandCentroidFile = np.loadtxt(pdb.replace('.pdb', '_centroidNetLigand'), dtype=np.str)

        ligand_mode = False  # whether ligands exist or not
        if len(ligand_data) != 0:
            ligand_mode = True
    except FileNotFoundError:
        print('Files not found')

    # Fill rsa, secondary_structure with ligands; if ligand_mode is False set up empty list
    if ligand_mode:
        ligands = np.unique(ligand_data[:, 1])
        for ligand in ligands:
            rsa = np.append(rsa, [[ligand, '-1']], axis=0)
            secStructure = np.append(secStructure, [[ligand, 'ligand', '360']], axis=0)
    else:
        ligands = []

    originalNodes = np.unique(data_all[:, 0:2])  # Define nodes (amino acids) that are interacting in net file

    # Process data of interactions
    basedata = add_terminal_atoms(data_all)
    basedata = make_undirected_interactions(basedata)
    influencedata = exclude_neighbour_nodes(basedata)
    influencenet = Graph.TupleList(edges=influencedata[:, :2], directed=False)
    influencenet_weight = 1 / influencedata[:, 2].astype(float)
    nodes = influencenet.vs['name']

    influencedata_directed = np.array([[], []], dtype='<U18').T  # empty np.array
    influencenet_directed = Graph.TupleList(edges=influencedata_directed[:, :2], directed=True)

    if not any(np.isin(nodes, secStructure[:, 0])):
        for node in nodes:
            if not (node in secStructure[:, 0]):
                secStructure = np.append(secStructure, [[node, 'xxx', '360']], axis=0)
            if not (node in rsa[:, 0]):
                rsa = np.append(rsa, [[node, 0]], axis=0)

    """WALKTRAP CLUSTERING
    Using basedata with all connections create GRAPH. Use Walktrap algo to define community
    clusters for nodes and write it into `net_community_vec` dict. And calculate network attributes """
    net = Graph.TupleList(edges=basedata[:, :2], directed=False)
    net_weight = 1 / basedata[:, 2].astype(float)
    net_community = net.community_walktrap(weights=net_weight)
    cluster = net_community.as_clustering()
    net_community_vec = {node: cluster for node, cluster in zip(cluster.graph.vs['name'], cluster.membership)}

    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_weight)
    node_edge_betweenness_wt = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_wt = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_wt, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_wt, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)
    secondOrder_node_intermodular_degree_wt, secondOrder_node_intermodular_degree_sidechain_wt = \
        node_intermodular_dgr_2order(node_intermodular_degree_wt, node_intermodular_degree_wt,
                                     node_modules, node_modules_sidechain, method='centroid')

    firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain = \
        first_second_order_degree(influencenet, influencenet_directed)

    """SECONDARY STRUCTURE CLUSTERING
    Use Stride created secStructure file to define community clusters for nodes.
    Based on type of secondary structure (structs) create `net_community_vec` dict. And calculate network attributes """
    structs = sorted(np.unique(secStructure[:, 1]))
    structs_dict = {factor: num for num, factor in enumerate(structs, 1)}
    net_community_vec = {node: structs_dict[f] for node, f in zip(secStructure[:, 0], secStructure[:, 1])}

    node_edge_betweenness_stride = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_stride = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_stride, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_stride, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)
    secondOrder_node_intermodular_degree_stride, secondOrder_node_intermodular_degree_sidechain_stride = \
        node_intermodular_dgr_2order(node_intermodular_degree_stride, node_intermodular_degree_stride, node_modules, node_modules_sidechain, method='centroid')

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
    secondOrder_node_intermodular_degree_sidechain_wt_df = create_df(secondOrder_node_intermodular_degree_sidechain_wt, nodes)  # ???

    out = [firstOrderDegree_df, firstOrderDegree_sidechain_df,
           secondOrderDegree_df, secondOrderDegree_sidechain_df,
           node_edge_betweenness_stride_df, node_edge_betweenness_sidechain_stride_df,
           node_intermodular_degree_stride_df, node_intermodular_degree_sidechain_stride_df,
           secondOrder_node_intermodular_degree_stride_df, secondOrder_node_intermodular_degree_sidechain_stride_df,
           node_edge_betweenness_wt_df, node_edge_betweenness_sidechain_wt_df,
           node_intermodular_degree_wt_df, node_intermodular_degree_sidechain_wt_df,
           secondOrder_node_intermodular_degree_wt_df, secondOrder_node_intermodular_degree_sidechain_wt_df]

    final_df = reduce(lambda left, right: pd.merge(left, right, on=0), out)
    colnames = ['AA', 'Degree', 'Degree_uniqueAtoms', 'SecondOrderDegree', 'SecondOrderDegree_uniqueAtoms',
                'NodeEdgeBetweennessSTRIDE', 'NodeEdgeBetweennessSTRIDE_uniqueAtoms', 'IntermodularDegreeSTRIDE',
                'IntermodularDegreeSTRIDE_uniqueAtoms', 'SecondOrderIntermodularDegreeSTRIDE',
                'SecondOrderIntermodularDegreeSTRIDE_uniqueAtoms', 'NodeEdgeBetweennessWALKTRAP',
                'NodeEdgeBetweennessWALKTRAP_uniqueAtoms', 'IntermodularDegreeWALKTRAP',
                'IntermodularDegreeWALKTRAP_uniqueAtoms', 'SecondOrderIntermodularDegreeWALKTRAP',
                'SecondOrderIntermodularDegreeWALKTRAP_uniqueAtoms']
    final_df.columns = colnames

    # Add nodes that were not networked in influencenet graph
    zeroNodes = [node for node in originalNodes if node not in nodes]
    zeroMat = pd.DataFrame(0, index=np.arange(len(zeroNodes)), columns=final_df.columns[0:])
    zeroMat['AA'] = zeroNodes
    output_df = pd.concat([final_df, zeroMat])

    # Check whether each acid has min distance to ligand
    if ligand_mode is True:
        ligandvec = pd.DataFrame(ligandCentroidFile)
        ligandvec[2] = ligandvec[2].astype(float)
        ligandvec = ligandvec.groupby(0)[2].min().to_dict()
    elif ligand_mode is False:
        ligandvec = {node: 0 for node in nodes}

    # Extract rsa info only for nodes that are in output_df
    rsa_list = []
    for acid in output_df['AA'].to_list():
        if acid in rsa[:, 0]:
            rsa_list.append(rsa[rsa[:, 0] == acid][0])
    rsa_df = pd.DataFrame(rsa_list, columns=['AA', 'RSA'])
    rsa_df['RSA'] = rsa_df['RSA'].astype(float)

    ligand_list = []
    for acid in output_df['AA'].to_list():
        if acid in ligandvec:
            ligand_list.append([acid, float(ligandvec[acid])])
        else:
            ligand_list.append([acid, 0])
    ligand_df = pd.DataFrame(ligand_list, columns=['AA', 'Ligand'])

    # Add ligand and rsa dataframes
    output_df = pd.merge(output_df, ligand_df, on='AA')
    output_df = pd.merge(output_df, rsa_df, on='AA')

    # Add 1 to everything to avoid zeros, except ligand
    output_df = output_df.fillna(value=0)
    output_df.loc[output_df['Ligand'] == 0, 'Ligand'] = 150  # set for Nan ligand 150 value
    output_df1 = output_df[output_df.columns[1:-1]] + 1

    final_out = pd.concat([output_df[['AA', 'RSA']], output_df1], axis=1)
    for col in final_out.columns[1:]:
        print(col, final_out[col].sum())

    outZ = Scaling(final_out[final_out.columns[1:]])
    outZ = pd.concat([final_out[['AA']], outZ], axis=1)

    final_out.to_csv(pdb.replace('.pdb', '_scoresCentroid'), sep='\t', index=False)
    outZ.to_csv(pdb.replace('.pdb', '_scoresCentroidZ'), sep='\t', index=False)
