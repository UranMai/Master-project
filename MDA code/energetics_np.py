# before create frame_1000.pdb and its files
# run >> python energetics_np.py ../Test_data/frame_1000.pdb
# it reads frame_1000_net, frame_1000_secondaryStructure2, frame_1000.rsa

import numpy as np
import pandas as pd 
import re, os
import sys
from functools import reduce 
import igraph
from sklearn.preprocessing import scale
from score_functions import *

def process_data(data_net): 
    '''
    @description
        Create data of reverse interactions, Merge directed and reversed data (bi-directed interactions)
        Example: 
                data_net:       ['GLU:22:A', 'ARG:231:A', 'SCSC', '10', 'SB', 'OE1', 'NH1']
                reverse_data:   ['ARG:231:A', 'GLU:22:A', 'SCSC', '10', 'SB', 'NH1', 'OE1']

        For merged data add info about whether acids' atoms are terminal(1) or not terminal(0)
        In merged data for groups of 2 common acids, count total energy (edgeSUM), count chains and number of bonds b/w 2 acids
        Consider data with poistive energies
    @input
        data_net - np.array data of interactions ('net' file)
        data_net must contain info in format 
                [acid1, acid2, Chain types, Energy_value, Bond_type, atom1, atom2]
                ['GLU:22:A', 'ARG:231:A', 'SCSC', '10', 'SB', 'OE1', 'NH1']
    @output
        dict of np.arrays, where
            'out' - grouped data of unique processed acids with counted options
            'influence_all' - merged data with at least one terminal atom
    '''
    # Edges should be listed bi-directionally
    reverse_data = data_net.copy()
    mask1 = np.where(reverse_data[:, 2]=='MCSC')
    mask2 = np.where(reverse_data[:, 2]=='SCMC')
    reverse_data[mask1, 2] = 'SCMC'
    reverse_data[mask2, 2] = 'MCSC'
    if data_net.shape[1] == 7:
        reverse_data[:,[0, 1, 5, 6]] = reverse_data[:, [1, 0, 6, 5]]
    elif data_net.shape[1] == 9:
        reverse_data[:,[0, 1, 5, 6, 7, 8]] = reverse_data[:, [1, 0, 6, 5, 8, 7]] 
    
    data = np.vstack((data_net, reverse_data))
    data = np.unique(data, axis=0)

    # Check if acids' atoms are terminal or not
    if data.shape[1] == 7:
        codes1 = []
        codes2 = []
        for res, atom in data[:, [0, 5]]:
            codes1.append(1 if atom in terminalAtoms[res[:3]] else 0)
        for res, atom in data[:, [1, 6]]:
            codes2.append(1 if atom in terminalAtoms[res[:3]] else 0)
        data = np.insert(data, 7, codes1, axis=1)
        data = np.insert(data, 8, codes2, axis=1)
        out_details = data.copy()
    else:
        out_details = data.copy()
    
    # Take bonds with positive energy and at least with one terminal atom (1)
    influence_all = out_details[(out_details[:,3].astype(float)>0) & np.bitwise_or(out_details[:,7]=='1', out_details[:,8]=='1')]

    # Group data by common 2 acids, calculate WeightSUM, percentMC
    out = []
    processed_acids = []
    for res1, res2 in data[:, [0,1]]:
        if [res1, res2] not in processed_acids:
            processed_acids.append([res1, res2])
            grouped_data = data[(data[:, 0]==res1) & (data[:, 1]==res2)]
            edgeSum = sum(grouped_data[:, 3].astype('float')) # weight_sum
            # chain = np.unique(grouped_data[:, 2], return_counts=True)
            # chain = dict(zip(chain[0], chain[1]))
            # # can we just use value_counts w/o calculating all numbers of sidechain
            # MCMCcount=0; MCSCcount=0; SCSCcount=0; SCMCcount=0
            # for tp in grouped_data[:, 2]:
            #     if tp=='MCMC' or tp=='MCSC':
            #         MCMCcount+=1
            #         MCSCcount+=1
            #     elif tp=='SCMC':
            #         SCMCcount+=1
            #     elif tp=='SCSC':
            #         SCSCcount+=1
            # percentMC = (MCMCcount+MCSCcount) / (MCMCcount+MCSCcount+SCMCcount+SCSCcount)
            # Header: Acid1, Acid2, mixed, Sum of Edges, percent MC, num1 of terminal atoms, num2 of terminal atoms, degree of Acid1
            values = [grouped_data[0,0], grouped_data[0,1], 'mixed', edgeSum, 
                      # percentMC,
                      # sum(grouped_data[:, 7].astype(int)), sum(grouped_data[:, 8].astype(int)), 
                      # len(data[data[:, 0] == grouped_data[0, 0]])
                      ] 
            out.append(values)
    out = np.array(out)
    out = out[out[:,3].astype(float)>0] # Consider only interactions with positive weightSum

    results = {}
    results['out'] = out
    results['influence_all'] = influence_all
    return results
    
def removeRedundancy(data):
    '''
    @description
        In input data first 2 columns are interacted amino acids, so select unique pairs of these acids. Traverse pairs
        Select subdata grouped by 2 acids and calculate sumWeight*0.5 of all interactions
        Example: 
                ['ALA:100:A', 'ALA:99:A', 'mixed', '7.964859491314238', '0.0', '2', '0', '22']
                ['ALA:99:A', 'ALA:100:A', 'mixed', '7.964859491314238', '1.0', '0', '2', '40']
        pairs returned are ['ALA:100:A', 'ALA:99:A']
    @input
        data - np.array, data of interactions 
    @output
        np.array in format [acid1, acid2, SumWeight]
        Example:
                ['ALA:100:A', 'ALA:99:A', '3.982429745657119']
    '''
    pairs = np.sort(data[:,:2], axis=1)
    pairs = np.unique(pairs, axis=0)
    out = []
    for i in range(len(pairs)):
        subdata = data[(data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1])]
        out.append([pairs[i,0], pairs[i,1], 0.5*sum(subdata[:,3].astype(float))])
        # Does this 0.5 matter?

    return np.array(out)

def make_undirected_interactions(data):
    '''
    @description
        Consider bidirected interactions: from input data sort acids by name and select unique pairs,
        Select subdata grouped by unique pair of2 acids and calculate sumWeight
    @input
        data - np.array, data of interactions
    @output
        np.array in format [Acid1, Acid2, SumWeight]
    '''
    pairs = np.sort(data[:,:2], axis=1)
    pairs = np.unique(pairs, axis=0)
    out = []
    for i in range(len(pairs)):
        subdat = data[((data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1])) | ((data[:,1]==pairs[i,0]) & (data[:,0]==pairs[i,1]))]
        out.append([subdat[0,0], subdat[0,1], sum(subdat[:,3].astype(float))])
    return np.array(out)

def make_directed_interactions(data):
    '''
    @description 
        Consider directed interactions, select unique pairs, 
        Select subdata grouped by unique pair of2 acids and calculate sumWeight
    @input
        data - np.array, data of interactions
    @output
        np.array in format [Acid1, Acid2, SumWeight]
    '''
    # pairs = np.sort(data[:,:2], axis=1) # is it important?
    pairs = np.unique(data[:,:2], axis=0)
    out = []
    for i in range(len(pairs)):
        subdat = data[((data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1]))]
        out.append([subdat[0,0], subdat[0,1], sum(subdat[:,3].astype(float))])
    return np.array(out)


if __name__ == '__main__':
    pdb = sys.argv[1] # input pdb file
    try: 
        secStructure = np.loadtxt(pdb.replace('.pdb', '_secondaryStructure2'), dtype=np.str)
        rsa = np.loadtxt(pdb.replace('.pdb', '.rsa'), dtype=np.str)
        data_all = np.loadtxt(pdb.replace('.pdb', '_net'), dtype=np.str)
    except FileNotFoundError:
        print('File not accesible')

    col0 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 0])))
    col1 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 1])))
    data_all[:, 0] = col0[:, 0]
    data_all[:, 1] = col1[:, 0]
    data_all = np.concatenate((data_all, col0[:, [1]]), axis=1)
    data_all = np.concatenate((data_all, col1[:, [1]]), axis=1)
    # data_all must be in format ['LYS:6:A' 'ASP:9:A' 'SCSC' '10' 'SB' 'NZ' 'OD2']
    # where columns 0,1-aminoacids; 2-chain type; 3-energy value; 4-bond type; 5,6-atoms


    # Define nodes (aminoacids) that are interacting in net file
    nodes = np.unique(data_all[:,0:2])
    originalNodes = nodes

    basedata = process_data(data_all)
    basedata_noPP = removeRedundancy(basedata['out'])

    # Consider only sidechain interactions w/o
    # basedata_allsidechain = make_undirected_interactions(data_all[data_all[:,2]!='MCMC'])

    influencedata_all = basedata['influence_all']
    influencedata_collapse = process_data(influencedata_all)
    influencedata = removeRedundancy(influencedata_collapse['out'])


    influencenet = igraph.Graph().TupleList(edges=influencedata[:,:2], directed=False)
    influencenet_W = 1/influencedata[:,2].astype(float) # change 3 to 2

    # Define directed network
    influencedata_directed = make_directed_interactions(influencedata_all[influencedata_all[:,7]=='1'])
    influencenet_directed = igraph.Graph().TupleList(edges=influencedata_directed[:,[1,0]], directed=True)
    influencenet_directed_W = 1/influencedata_directed[:,2].astype(float) 
 
    # Fill secStructure and rsa data with missing nodes
    if any(np.isin(nodes, secStructure[:,0])) == False:
        for node in nodes:
            if (node in secStructure[:,0])==False:
                secStructure = np.append(secStructure, [[node, 'xxx', '360']], axis=0)
        # in R code used rownames(rsa)
        for node in nodes:
            if (node in rsa[:,0]) == False:
                rsa = np.append(rsa, [[node, 0]], axis=0)

    # Define undirected network
    net = igraph.Graph().TupleList(edges=basedata_noPP[:,:2], directed=False)
    net_W = 1/basedata_noPP[:,2].astype(float) 


    # -------------------------WALKTRAP CLUSTERING-------------------------------------------------
    # Use Walktrap algo to define clusters for nodes and write it into net_community_vec dictionary
    net_community = net.community_walktrap(net_W, steps=4)
    cluster = net_community.as_clustering()
    net_community_vec = {name:member for name, member in zip(cluster.graph.vs['name'], cluster.membership)}
    nodes = net.vs['name']

    # Network attributes calculation
    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_W)

    # WALKTRAP community vec
    node_edge_betweenness_wt = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_wt = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_wt, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_wt, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)

    secondOrder_node_intermodular_degree_wt, \
    secondOrder_node_intermodular_degree_sidechain_wt = node_intermodular_degr_2order(node_intermodular_degree_wt, 
                                                                                        node_intermodular_degree_wt, 
                                                                                        node_modules, 
                                                                                        node_modules_sidechain, 
                                                                                        method='energetics')

    firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain = first_second_order_degree(influencenet, influencenet_directed)


    # ------------------------SecondaryStructure CLUSTERING---------------------------------------
    # Use SecondaryStructure file to define clusters for nodes
    factors = sorted(np.unique(secStructure[:,1]))
    factors_dict = {factor:num for num, factor in enumerate(factors, 1)}
    net_community_vec = {node:factors_dict[f] for node, f in zip(secStructure[:,0], secStructure[:,1])}

    # Network attributes calculation
    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_W)

    # SECONDARY STRUCTURE network vec     
    node_edge_betweenness_stride = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_stride = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_stride, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_stride, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)

    secondOrder_node_intermodular_degree_stride, \
    secondOrder_node_intermodular_degree_sidechain_stride = node_intermodular_degr_2order(node_intermodular_degree_stride, 
                                                                                            node_intermodular_degree_stride, 
                                                                                            node_modules, 
                                                                                            node_modules_sidechain,
                                                                                            method='energetics')
    

    nodes = influencenet.vs['name']

    # Create dataframes from dictionaries and merge them
    firstOrderDegree_df = Nanfill_df(firstOrderDegree, nodes)
    firstOrderDegree_sidechain_df = Nanfill_df(firstOrderDegree_sidechain, nodes)

    secondOrderDegree_df = Nanfill_df(secondOrderDegree, nodes)
    secondOrderDegree_sidechain_df = Nanfill_df(secondOrderDegree_sidechain, nodes)

    node_edge_betweenness_stride_df = Nanfill_df(node_edge_betweenness_stride, nodes)
    node_edge_betweenness_sidechain_stride_df = Nanfill_df(node_edge_betweenness_sidechain_stride, nodes)

    node_intermodular_degree_stride_df = Nanfill_df(node_intermodular_degree_stride, nodes)
    node_intermodular_degree_sidechain_stride_df =Nanfill_df(node_intermodular_degree_sidechain_stride, nodes)

    secondOrder_node_intermodular_degree_sidechain_stride_df = Nanfill_df(secondOrder_node_intermodular_degree_sidechain_stride, nodes)
    secondOrder_node_intermodular_degree_stride_df = Nanfill_df(secondOrder_node_intermodular_degree_stride, nodes)

    node_edge_betweenness_wt_df = Nanfill_df(node_edge_betweenness_wt, nodes)
    node_edge_betweenness_sidechain_wt_df = Nanfill_df(node_edge_betweenness_sidechain_wt, nodes)

    node_intermodular_degree_wt_df = Nanfill_df(node_intermodular_degree_wt, nodes)
    node_intermodular_degree_sidechain_wt_df = Nanfill_df(node_intermodular_degree_sidechain_wt, nodes)

    secondOrder_node_intermodular_degree_wt_df = Nanfill_df(secondOrder_node_intermodular_degree_wt, nodes)
    secondOrder_node_intermodular_degree_sidechain_wt_df = Nanfill_df(secondOrder_node_intermodular_degree_sidechain_wt, nodes) #???


    out = [
           firstOrderDegree_df, firstOrderDegree_sidechain_df, 
           secondOrderDegree_df, secondOrderDegree_sidechain_df,
           node_edge_betweenness_stride_df, node_edge_betweenness_sidechain_stride_df,           
           node_intermodular_degree_stride_df, node_intermodular_degree_sidechain_stride_df,
           secondOrder_node_intermodular_degree_stride_df, secondOrder_node_intermodular_degree_sidechain_stride_df,
           node_edge_betweenness_wt_df, node_edge_betweenness_sidechain_wt_df,
           node_intermodular_degree_wt_df, node_intermodular_degree_sidechain_wt_df,
           secondOrder_node_intermodular_degree_wt_df, secondOrder_node_intermodular_degree_sidechain_wt_df
           ]

    # Create final dataset, merging dataframes by common AA column
    df_final = reduce(lambda left,right: pd.merge(left,right,on=0), out)
    colnames = [
                'AA', 
                "Degree","Degree_sidechain",
                "SecondOrderDegree","SecondOrderDegree_sidechain",
                "NodeEdgeBetweennessSTRIDE","NodeEdgeBetweennessSTRIDE_sidechain",
                "IntermodularDegreeSTRIDE",  "IntermodularDegreeSTRIDE_sidechain",
                "SecondOrderIntermodularDegreeSTRIDE",  "SecondOrderIntermodularDegreeSTRIDE_sidechain",
                "NodeEdgeBetweennessWALKTRAP",  "NodeEdgeBetweennessWALKTRAP_sidechain",
                "IntermodularDegreeWALKTRAP",  "IntermodularDegreeWALKTRAP_sidechain",
                "SecondOrderIntermodularDegreeWALKTRAP",  "SecondOrderIntermodularDegreeWALKTRAP_sidechain"
                ]     
    df_final.columns = colnames

    # Create zero dataframe with missing nodes 
    nodes = influencenet.vs['name']
    zeroNodes = []
    for node in originalNodes:
        if node not in nodes:
            zeroNodes.append(node)
    zeroMat = pd.DataFrame(0, index=np.arange(len(zeroNodes)), columns=df_final.columns[0:])
    zeroMat['AA'] =zeroNodes

    # Concat final dataset with zero dataframe, replace Nones to 0, and Add 1 to everything to avoid zeros, except ligand
    out = pd.concat([df_final, zeroMat])
    out1 = out.fillna(value=0)
    out1 = out1[out1.columns[1:]]+1
    out = pd.concat([out[['AA']], out1], axis=1)

    # Check results
    for col in out.columns[1:]:
        print(col, out[col].sum())


    out1 = Scaling(out[out.columns[1:]])
    outZ = pd.concat([out[['AA']], out1], axis=1)

    out.to_csv(re.sub('.pdb', '', pdb)+'_scoresEnergetics', sep='\t', index=False)
    outZ.to_csv(re.sub('.pdb', '', pdb)+'_scoresEnergeticsZ', sep='\t', index=False)