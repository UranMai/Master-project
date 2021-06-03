import os 
import sys
import igraph
import numpy as np
import pandas as pd
# combine functions for calculate graph attributes
# 

def node_edge_btwns(data, edge_btwns, nodes, community):
    '''
    @input
        influencedata - processed data of interactions
        edge_betweenness - calculated igraph.edge_betweenness for network
        nodes - acids
        net_community_vec - nodes belong to clusters defined using walktrap or secStructure clusters
    @description
        calculate edge_betweenness for nodes(acids) from influencedata
        nodes must be from different clusters(communities)
    @return 
        dict of nodes and edge_betweenness
    '''
    node_edge_betweenness = {node:0 for node in nodes}
    for i in range(len(edge_btwns)):
        if community[data[i,0]] != community[data[i,1]]:
            node_edge_betweenness[data[i,0]] = node_edge_betweenness[data[i,0]] + edge_btwns[i]
            node_edge_betweenness[data[i,1]] = node_edge_betweenness[data[i,1]] + edge_btwns[i]
    return node_edge_betweenness

def node_intermodular_degree(data, nodes, community):
    '''
    @input
        influencedata - processed data of interactions
        nodes - acids
        net_community_vec - nodes belong to clusters defined using walktrap or secStructure clusters
    @description
        node_intermodular_degree - number of interactions b/w acids from different clusters
        select one node, for it select subdata of interactions with other nodes, 
        and count the number of nodes from other clusters different from the node's one
    @return
        node_intermodular_degree - dict of nodes and their number of interclustered nodes
        node_modules - dict of nodes and their interclustered nodes with cluster number
    ''' 
    node_intermodular_degree = {node:0 for node in nodes}
    node_modules = {}
    for node in nodes:
        subdata = data[(data[:,0]==node) | (data[:,1]==node)]
        if len(subdata)==0:
            node_intermodular_degree[node]=0
        
        bound = np.unique(subdata[:,0:2])
        bound = bound[bound!=node]        
        bound_modules = {}
        for bond in bound:
            if community[bond] != community[node]:
                bound_modules[bond] = community[bond]
        node_modules[node] = bound_modules
        node_intermodular_degree[node] = len(bound_modules)
    return node_intermodular_degree, node_modules


def node_edge_btwns_SC(data_directed, data, edge_btwns, nodes, community):
    '''
    @input
        influencedata_directed - processed data of interactions SideChain
        influencedata - processed data of interactions 
        edge_betweenness - calculated igraph.edge_betweenness for network
        nodes - acids
        net_community_vec - nodes belong to clusters defined using walktrap or secStructure clusters
    @description
        calculate edge_betweenness for nodes(acids) from influencedata weighted by sidechain
        nodes must be from different clusters(communities)
    @return 
        dict of nodes and edge_betweenness weighted by sidechain
    '''
    node_edge_betweenness_sidechain = {node:0 for node in nodes}
    for i in range(len(edge_btwns)):
        if community[data[i,0]] != community[data[i,1]]:
            subdata = data_directed[(data_directed[:,0]==data[i,0]) & (data_directed[:,1]==data[i,1])]
            if len(subdata)==0:
                weight1 = 0
            else:
                weight1 = subdata[0,2].astype(float) #change 3 to 2

            subdata = data_directed[(data_directed[:,0]==data[i,1]) & (data_directed[:,1]==data[i,0])]
            if len(subdata)==0:
                weight2 = 0
            else:
                weight2 = subdata[0,2].astype(float) # change 3 to 2
            # it is important to process first zero weights than non_zero weights why?
            # how to handle with division zero by zero
            if weight1 == 0 and weight2 == 0:
                node_edge_betweenness_sidechain[data[i,0]] = node_edge_betweenness_sidechain[data[i,0]] + 0
                node_edge_betweenness_sidechain[data[i,1]] = node_edge_betweenness_sidechain[data[i,1]] + 0 
            else:
                node_edge_betweenness_sidechain[data[i,0]] = node_edge_betweenness_sidechain[data[i,0]]+ weight1/(weight1+weight2) * edge_btwns[i]
                node_edge_betweenness_sidechain[data[i,1]] = node_edge_betweenness_sidechain[data[i,1]]+ weight2/(weight1+weight2) * edge_btwns[i]
    return node_edge_betweenness_sidechain

def node_intermodular_degree_SC(data_directed, nodes, community):
    '''
    @input
        influencedata_directed - processed data of interactions SideChain
        nodes - acids
        net_community_vec - nodes belong to clusters defined using walktrap or secStructure clusters
    @description
        node_intermodular_degree - number of interactions b/w acids from different clusters
        select one node, for it select subdata of interactions with other nodes, 
        and count the number of nodes from other clusters different from the node's one
    @return
        node_intermodular_degree_sidechain - dict of nodes and their number of interclustered nodes
        node_modules_sidechain - dict of nodes and their interclustered nodes with cluster number
    ''' 
    node_intermodular_degree_sidechain = {node:0 for node in nodes}
    node_modules_sidechain = {}
    for node in nodes:
        subdata = data_directed[data_directed[:,0]==node]
        if len(subdata) == 0:
            node_intermodular_degree_sidechain[node] = 0
        
        bound = np.unique(subdata[:,0:2])
        bound = bound[bound!=node]
        bound_modules = {}
        for bond in bound:
            if community[bond] != community[node]:
                bound_modules[bond] = community[bond]
        node_modules_sidechain[node] = bound_modules
        node_intermodular_degree_sidechain[node] = len(bound_modules)
    return node_intermodular_degree_sidechain, node_modules_sidechain

def node_intermodular_degr_2order(node_intmod_degr, node_intmod_degr_sc, node_modules, node_modules_sc, method='centroid'):
    '''
    @input
        node_intermodular_degree - dict of nodes and their number of interclustered nodes
        node_intermodular_degree_sidechain - dict of nodes and their number of interclustered nodes (SideChain)
        node_modules - dict of nodes and their interclustered nodes with cluster number
        node_modules_sidechain - dict of nodes and their interclustered nodes with cluster number (SideChain)
    @description
        Find number of secondary interactions b/w acids from different clusters
        Select node, go through its interclustered nodes and count interinterclustered nodes
    @return
        dict of nodes and number(degree) of inter-interclustered nodes
        and also for SideChain
    '''
    secondOrder_node_intermodular_degree = {node:0 for node in node_intmod_degr.keys()}
    secondOrder_node_intermodular_degree_sidechain = {node:0 for node in node_intmod_degr_sc.keys()}

    for node in node_intmod_degr.keys():
        if node_intmod_degr[node]==0:
            continue
        secondOrder_degree = []
        for node1 in node_modules[node]:
            for node2 in node_modules[node1]:
                secondOrder_degree.append(node2)
        if method == 'centroid':
            secondOrder_node_intermodular_degree[node] = len(secondOrder_degree) 
        else:
            secondOrder_node_intermodular_degree[node] = len(secondOrder_degree) - 1

        secondOrder_degree_sidechain = []
        for node1 in node_modules_sc[node]:
            for node2 in node_modules_sc[node1]:
                secondOrder_degree_sidechain.append(node2)
        if method == 'centroid':
            secondOrder_node_intermodular_degree_sidechain[node] = len(secondOrder_degree_sidechain) 
        else:
            secondOrder_node_intermodular_degree_sidechain[node] = len(secondOrder_degree_sidechain) - 1

    return secondOrder_node_intermodular_degree, secondOrder_node_intermodular_degree_sidechain

# def calc_secondorder_degree(influencenet, influencenet_directed, nodes): # for energetics
#     '''
#     @input
#         influencenet - 
#         influencenet_directed - 
#         nodes - acids
#     @description
#         Calculate nodes' degree (number of adjacent nodes) as graph.degree() function; mode='in' means ...
#         For each node iterate through its first-order neighboring nodes, then second-order neighboring nodes
#         And cacluate the number of second-order interactions
#     @return
#         dicts of 1st-order, 2nd-order degree (also for SideChain)
#     '''
#     firstOrderDegree = {name:degree for name, degree in zip(influencenet.vs['name'], influencenet.degree())}
#     firstOrderDegree_sidechain = {name:degree for name, degree in zip(influencenet_directed.vs['name'], influencenet_directed.degree(mode='in'))}

#     secondOrderDegree = {}
#     secondOrderDegree_sidechain = {}
#     nodes = influencenet.vs['name']
#     for num, node in enumerate(nodes):
#         firstorder = np.array(influencenet.neighbors(node))
#         secondorder = []
#         for neigh in firstorder:
#             secondorder = np.append(secondorder, influencenet.neighbors(neigh))
#             secondorder = np.unique(secondorder)
#             secondorder = secondorder[secondorder!=num]
#             secondOrderDegree[node] = len(secondorder)

#         firstorder = np.array(influencenet_directed.neighbors(node, mode='in'))
#         secondorder = []
#         for neigh in firstorder:
#             secondorder = np.append(secondorder, influencenet_directed.neighbors(neigh, mode='in'))
#             secondorder = np.unique(secondorder)
#             num1 = influencenet_directed.vs['name'].index(node)
#             secondorder = secondorder[secondorder!=num1]
#             secondOrderDegree_sidechain[node] = len(secondorder)

#     return firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain

def first_second_order_degree(net, net_directed): # same output for energetics and centroid versions
    '''
    @input
        influencenet - 
        influencenet_directed - 
        nodes - acids
    @description
        Calculate nodes' degree (number of adjacent nodes) as graph.degree() function; mode='in' means ...
        For each node iterate through its first-order neighboring nodes, then second-order neighboring nodes
        And cacluate the number of second-order interactions
    @return
        dicts of 1st-order, 2nd-order degree (also for SideChain)
    '''
    firstOrderDegree = {name:degree for name, degree in zip(net.vs['name'], net.degree())}
    firstOrderDegree_sidechain = {name:degree for name, degree in zip(net_directed.vs['name'], net_directed.degree(mode='in'))}

    secondOrderDegree = {}
    secondOrderDegree_sidechain = {}
    nodes = net.vs['name']
    for num, node in enumerate(nodes):
        firstorder = net.neighbors(node)
        secondorder = []
        for neigh in firstorder:
            secondorder = np.append(secondorder, net.neighbors(neigh))
            secondorder = np.unique(secondorder)
            secondorder = secondorder[secondorder!=num]
            secondOrderDegree[node] = len(secondorder)
        if node in net_directed.vs['name']:
            firstorder_sidechain = net_directed.neighbors(node, mode='in')
            num1 = net_directed.vs['name'].index(node)
        else:
            firstorder_sidechain = []
        secondorder = []
        for neigh in firstorder_sidechain:
            secondorder = np.append(secondorder, net_directed.neighbors(neigh, mode='in'))
            secondorder = np.unique(secondorder)
            secondorder = secondorder[secondorder!=num1]
            secondOrderDegree_sidechain[node] = len(secondorder)

    return firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain


def Nanfill_df(dictionary, nodes):
    # Fill data with None values for missing nodes
    # nodes = influencenet.vs['name']
    df = pd.DataFrame.from_dict(data=dictionary.items())
    for node in nodes:
        if node not in df[0].values:
            dict1 = {0:node, 1:None}
            df = df.append(dict1, ignore_index=True)
    return df

def Scaling(df):
    # Apply standard normalization, ddof=1
    cols = df.columns
    for col in cols:
        std = np.std(df[col], ddof=1)
        mean = np.mean(df[col])
        df[col] = (df[col] - mean)/std
    return df

# used for centroid calculations
# uniqueAtoms = {
#                 'ALA': ['CB,HB1,HB2,HB3'],
#                 'ARG': ['CD,HD2,HD3,NE,HE,CZ,NH1,NH2,HH11,HH12,HH21,HH22,NH1/2'],
#                 'ASN': ['OD1,ND2,HD21,HD22', 'CG'],
#                 'ASP': ['OD1,OD2,HD2,OD1/2', 'CG'],
#                 'CYS': ['SG,HG'],
#                 'GLN': ['OE1,NE2,HE21,HE22,CD'],
#                 'GLU': ['OE1,OE2,HE2,OE1/2', 'CD'],
#                 'GLY': ['CA,HA2,HA3,O,N,C'],
#                 'HIS': ['ND1,CE1,NE2,HD1,HE1,HE2,CD2,HD2'],
#                 'ILE': ['CD1,HD11,HD12,HD13,CG2,HG21,HG22,HG23'],
#                 'LEU': ['CD1,CD2,HD11,HD12,HD13,HD21,HD22,HD23'],
#                 'LYS': ['CE,HE2,HE3,NZ,HZ1,HZ2,HZ3'],
#                 'MET': ['CE,HE1,HE2,HE3,SD'],
#                 'PHE': ['CZ,HZ,CE1,CE2,HE1,CE2,HE2'],
#                 'PRO': ['CG,CD,CB,HG2,HG3,HD2,HD3,HB2,HB3'],
#                 'SER': ['OG,HG,CB,HB2,HB3'],
#                 'THR': ['OG1,HG1,CG2,HG21,HG22,HG23'],
#                 'TRP': ['CD2,NE1,HE1,CE2,CE3,HE3,CZ2,HZ2,CZ3,HZ3,CH2,HH2'],
#                 'TYR': ['OH,HH,CE1,HE1,CZ,CE2,HE2'],
#                 'VAL': ['CG1,CG2,HG11,HG12,HG13,HG21,HG22,HG23']
#                 }
uniqueAtoms = {
    'ALA': ['CB', 'HB1', 'HB2', 'HB3'],
    'ARG': ['CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'NH2', 'HH11', 'HH12', 'HH21', 'HH22', 'NH1/2'],
    'ASN': ['OD1', 'ND2', 'HD21', 'HD22', 'CG'],
    'ASP': ['OD1', 'OD2', 'HD2', 'OD1/2', 'CG'],
    'CYS': ['SG', 'HG'],
    'GLN': ['OE1', 'NE2', 'HE21', 'HE22', 'CD'],
    'GLU': ['OE1', 'OE2', 'HE2', 'OE1/2', 'CD'],
    'GLY': ['CA', 'HA2', 'HA3', 'O', 'N', 'C'],
    'HIS': ['ND1', 'CE1', 'NE2', 'HD1', 'HE1', 'HE2', 'CD2', 'HD2'],
    'ILE': ['CD1', 'HD11', 'HD12', 'HD13', 'CG2', 'HG21', 'HG22', 'HG23'],
    'LEU': ['CD1', 'CD2', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
    'LYS': ['CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'HZ2', 'HZ3'],
    'MET': ['CE', 'HE1', 'HE2', 'HE3', 'SD'],
    'PHE': ['CZ', 'HZ', 'CE1', 'CE2', 'HE1', 'CE2', 'HE2'],
    'PRO': ['CG', 'CD', 'CB', 'HG2', 'HG3', 'HD2', 'HD3', 'HB2', 'HB3'],
    'SER': ['OG', 'HG', 'CB', 'HB2', 'HB3'],
    'THR': ['OG1', 'HG1', 'CG2', 'HG21', 'HG22', 'HG23'],
    'TRP': ['CD2', 'NE1', 'HE1', 'CE2', 'CE3', 'HE3', 'CZ2', 'HZ2', 'CZ3', 'HZ3', 'CH2', 'HH2'],
    'TYR': ['OH', 'HH', 'CE1', 'HE1', 'CZ', 'CE2', 'HE2'],
    'VAL': ['CG1', 'CG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23']}
# used for energetics calculations
terminalAtoms = {
    'ALA': ['CB', 'HB1', 'HB2', 'HB3'],
    'ARG': ['NE', 'HE', 'CZ', 'NH1', 'NH2', 'HH11', 'HH12', 'HH21', 'HH22'],
    'ASN': ['OD', 'ND2', 'HD21', 'HD22'],
    'ASP': ['OD1', 'OD2', 'HD2'],
    'CYS': ['SG', 'HG'],
    'GLN': ['OE1', 'NE2', 'HE21', 'HE22'],
    'GLU': ['OE1', 'OE2', 'HE2'],
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