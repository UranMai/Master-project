import numpy as np
import pandas as pd


def node_edge_btwns(data, edge_btwns, nodes, community):
    """
    @description
        Edge_betweenness - the number of the shortest paths that go through an edge.
        Calculate edge_betweenness for nodes (acids) for `data`. Nodes must be from different community clusters
    @input
        data - np.array data of interactions
        edge_btwns - calculated igraph.edge_betweenness for network
        nodes - amino acids
        net_community_vec - nodes belong to clusters defined using Walktrap algo or secStructure clusters
    @output
        node_edge_betweenness - dict of nodes and edge_betweenness values
    """
    node_edge_betweenness = {node: 0 for node in nodes}
    for i in range(len(edge_btwns)):
        acid1, acid2 = data[i, 0], data[i, 1]
        if acid1 in community and acid2 in community:
            if community[acid1] != community[acid2]:
                node_edge_betweenness[acid1] = node_edge_betweenness[acid1] + edge_btwns[i]
                node_edge_betweenness[acid2] = node_edge_betweenness[acid2] + edge_btwns[i]
    return node_edge_betweenness


def node_intermodular_degree(data, nodes, community):
    """
    @description
        node_intermodular_degree - number of interactions b/w acids from different clusters
        Select one node, for it select sub_data of interactions with other nodes,
        and count the number of nodes from other clusters differed from the current node's cluster
    @input
        data - np.array data of interactions
        nodes - amino acids
        net_community_vec - nodes belong to clusters defined using Walktrap algo or secStructure clusters
    @output
        node_intermodular_degree - dict of nodes and their number of inter-clustered nodes
        node_modules - dict of nodes and their inter-clustered nodes with cluster number
    """
    node_intermodular_dgr = {node: 0 for node in nodes}
    node_modules = {}
    for node in nodes:
        sub_data = data[(data[:, 0] == node) | (data[:, 1] == node)]
        if len(sub_data) == 0:
            node_intermodular_dgr[node] = 0

        bound = np.unique(sub_data[:, 0:2])
        bound = bound[bound != node]
        bound_modules = {}
        for bond in bound:
            if bond in community and node in community:
                if community[bond] != community[node]:
                    bound_modules[bond] = community[bond]
        node_modules[node] = bound_modules
        node_intermodular_dgr[node] = len(bound_modules)
    return node_intermodular_dgr, node_modules


def node_edge_btwns_SC(data_directed, data, edge_btwns, nodes, community):
    """
    @description
        Calculate edge_betweenness for nodes (acids) weighted by side chain.
        Nodes must be from different community clusters.
    @input
        data_directed, data - np.array data of interactions
        edge_btwns - calculated igraph.edge_betweenness for network
        nodes - amino acids
        net_community_vec - nodes belong to clusters defined using Walktrap algo or secStructure clusters
    @output
        node_edge_betweenness_sidechain - dict of nodes and edge_betweenness weighted by side chain
    """
    node_edge_betweenness_sidechain = {node: 0 for node in nodes}
    for i in range(len(edge_btwns)):
        acid1, acid2 = data[i, 0], data[i, 1]
        if acid1 in community and acid2 in community:
            if community[acid1] != community[acid2]:
                sub_data = data_directed[(data_directed[:, 0] == acid1) & (data_directed[:, 1] == acid2)]
                if len(sub_data) == 0:
                    weight1 = 0
                else:
                    weight1 = sub_data[0, 2].astype(float)

                sub_data = data_directed[(data_directed[:, 0] == acid2) & (data_directed[:, 1] == acid1)]
                if len(sub_data) == 0:
                    weight2 = 0
                else:
                    weight2 = sub_data[0, 2].astype(float)

                if weight1 == 0 and weight2 == 0:
                    node_edge_betweenness_sidechain[acid1] = node_edge_betweenness_sidechain[acid1] + 0
                    node_edge_betweenness_sidechain[acid2] = node_edge_betweenness_sidechain[acid2] + 0
                else:
                    node_edge_betweenness_sidechain[acid1] = node_edge_betweenness_sidechain[acid1] + weight1 / (weight1 + weight2) * edge_btwns[i]
                    node_edge_betweenness_sidechain[acid2] = node_edge_betweenness_sidechain[acid2] + weight2 / (weight1 + weight2) * edge_btwns[i]
    return node_edge_betweenness_sidechain


def node_intermodular_degree_SC(data_directed, nodes, community):
    """
    @description
        node_intermodular_degree - number of interactions b/w acids from different clusters
        Select one node, for it select sub_data of interactions with other nodes,
        and count the number of nodes from other clusters differed from the current node's cluster
    @input
        data_directed - np.array data of interactions
        nodes - amino acids
        net_community_vec - nodes belong to clusters defined using Walktrap algo or secStructure clusters
    @output
        node_intermodular_degree_sidechain - dict of nodes and their number of inter-clustered nodes
        node_modules_sidechain - dict of nodes and their inter-clustered nodes with cluster number
    """
    node_intermodular_degree_sidechain = {node: 0 for node in nodes}
    node_modules_sidechain = {}
    for node in nodes:
        sub_data = data_directed[data_directed[:, 0] == node]
        if len(sub_data) == 0:
            node_intermodular_degree_sidechain[node] = 0

        bound = np.unique(sub_data[:, 0:2])
        bound = bound[bound != node]
        bound_modules = {}
        for bond in bound:
            if bond in community and node in community:
                if community[bond] != community[node]:
                    bound_modules[bond] = community[bond]
        node_modules_sidechain[node] = bound_modules
        node_intermodular_degree_sidechain[node] = len(bound_modules)
    return node_intermodular_degree_sidechain, node_modules_sidechain


def node_intermodular_dgr_2order(node_intermodular_degree, node_intermodular_degree_sidechain,
                                 node_modules, node_modules_sc, method='energetics'):
    """
    @description
        Find number of secondary order interactions b/w acids from different clusters
        Select node, go through its inter-clustered nodes. For these inter-clustered nodes
        count 2nd order inter-clustered nodes
    @input
        node_intermodular_degree - dict of nodes and their number of inter-clustered nodes
        node_intermodular_degree_sidechain - dict of nodes and their number of inter-clustered nodes (SideChained)
        node_modules - dict of nodes and their inter-clustered nodes with cluster number
        node_modules_sidechain - dict of nodes and their inter-clustered nodes with cluster number (SideChained)
        method - ['centroid', 'energetics']
    @output
        secondOrder_node_intermodular_degree - dict of nodes and number(degree) of inter-inter-clustered nodes
        secondOrder_node_intermodular_degree_sidechain - dict of nodes and number(degree) of inter-inter-clustered nodes for side chain
    """
    secondOrder_node_intermodular_degree = {node: 0 for node in node_intermodular_degree.keys()}
    secondOrder_node_intermodular_degree_sidechain = {node: 0 for node in node_intermodular_degree_sidechain.keys()}

    for node in node_intermodular_degree.keys():
        if node_intermodular_degree[node] == 0:
            continue
        secondOrder_degree = []
        for node1 in node_modules[node]:
            for node2 in node_modules[node1]:
                secondOrder_degree.append(node2)
        if method == 'centroid':
            secondOrder_node_intermodular_degree[node] = len(secondOrder_degree)
        elif method == 'energetics':
            secondOrder_node_intermodular_degree[node] = len(secondOrder_degree) - 1

        secondOrder_degree_sidechain = []
        for node1 in node_modules_sc[node]:
            for node2 in node_modules_sc[node1]:
                secondOrder_degree_sidechain.append(node2)
        if method == 'centroid':
            secondOrder_node_intermodular_degree_sidechain[node] = len(secondOrder_degree_sidechain)
        elif method == 'energetics':
            secondOrder_node_intermodular_degree_sidechain[node] = len(secondOrder_degree_sidechain) - 1
    return secondOrder_node_intermodular_degree, secondOrder_node_intermodular_degree_sidechain


def first_second_order_degree(net, net_directed):  # same output for energetics and centroid versions
    """
    @description
        Calculate nodes degree (number of adjacent nodes) as graph.degree() function; mode='in' for in-degrees
        in-degrees means that the edges must be oriented towards the root.
        For each node iterate through its first-order neighboring nodes, then second-order neighboring nodes
        And calculate the number of second-order interactions
    @input
        net, net_directed - igraph.Graph objects
        nodes - amino acids
    @output
        firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain -
        dicts of 1st-order, 2nd-order degrees (also for SideChain)
    """
    firstOrderDegree = {name: degree for name, degree in zip(net.vs['name'], net.degree())}
    firstOrderDegree_sidechain = {name: degree for name, degree in zip(net_directed.vs['name'], net_directed.degree(mode='in'))}

    secondOrderDegree = {}
    secondOrderDegree_sidechain = {}
    nodes = net.vs['name']
    for num, node in enumerate(nodes):
        first_order = net.neighbors(node)
        second_order = []
        for neigh in first_order:
            second_order = np.append(second_order, net.neighbors(neigh))
            second_order = np.unique(second_order)
            second_order = second_order[second_order != num]
            secondOrderDegree[node] = len(second_order)
        if node in net_directed.vs['name']:
            first_order_sidechain = net_directed.neighbors(node, mode='in')
            num1 = net_directed.vs['name'].index(node)
        else:
            first_order_sidechain = []

        second_order = []
        for neigh in first_order_sidechain:
            second_order = np.append(second_order, net_directed.neighbors(neigh, mode='in'))
            second_order = np.unique(second_order)
            second_order = second_order[second_order != num1]
            secondOrderDegree_sidechain[node] = len(second_order)

    return firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain


def create_df(dictionary, nodes):
    """
    Create dataframe or add rows for nodes not in dictionary
    """
    # Fill data with None values for missing nodes
    if len(dictionary) == 0:
        dct = {node: None for node in nodes}
        df = pd.DataFrame.from_dict(data=dct.items())
    else:
        df = pd.DataFrame.from_dict(data=dictionary.items())
        for node in nodes:
            if node not in df[0].values:
                dct = {0: node, 1: None}
                df = df.append(dct, ignore_index=True)
    return df


def Scaling(df):
    """Apply standard normalization (optional), ddof=1"""
    for col in df.columns:
        std = df[col].std()
        mean = df[col].mean()
        df[col] = (df[col] - mean) / std
    return df


# used for energetics calculations
terminalAtoms = {
    'ALA': ['CB', 'HB1', 'HB2', 'HB3'],
    'ARG': ['NE', 'HE', 'CZ', 'NH1', 'NH2', 'HH11', 'HH12', 'HH21', 'HH22', 'NH1/2'],
    'ASN': ['OD', 'ND2', 'HD21', 'HD22'],
    'ASP': ['OD1', 'OD2', 'HD2', 'OD1/2'],
    'CYS': ['SG', 'HG'],
    'GLN': ['OE1', 'NE2', 'HE21', 'HE22'],
    'GLU': ['OE1', 'OE2', 'HE2', 'OE1/2'],
    'GLY': ['CA', 'HA2', 'HA3', 'O'],
    'HIS': ['ND1', 'CE1', 'NE2', 'HD1', 'HE1', 'HE2'],
    'ILE': ['CD1', 'HD1'],
    'LEU': ['CD1', 'CD2', 'HD1', 'HD2'],
    'LYS': ['NZ', 'HD1', 'HD2', 'HD3'],
    'MET': ['CE', 'HE1', 'HE2', 'HE3'],
    'PHE': ['CZ', 'HZ1', 'HZ2', 'HZ3'],
    'PRO': ['CG', 'CD', 'HG2', 'HG3', 'HD2', 'HD3'],
    'SER': ['OG', 'HG'],
    'THR': ['OG1', 'HG1', 'CG2', 'HG1', 'HG2', 'HG3'],
    'TRP': ['CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2'],
    'TYR': ['OH', 'HH'],
    'VAL': ['CG1', 'CG2', 'HG1', 'HG2']}
