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
        Consider reverse interactions of edges and merge direct and reverse data (bi-directed interactions)
        For merged data add info about whether acids atoms are terminal(1) or not terminal(0)
        In merged data for groups of 2 common acids, count total E (sumWeight), count chains and number of bonds b/w 2 acids
    @input
        data_net - np.array data of interactions ('net' file)
        data_net must contain info in format [acid1, acid2, Chain types, Energy_value, Bond_type, atom1, atom2]
    @output
        dict of np.arrays, where
            'out' - data of unique processed acids with >0 WeightSum
            'influence_all' - data with at least one terminal atom
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
    
    data = np.concatenate((data_net, reverse_data), axis=0, )
    data = np.unique(data.tolist(), axis=0)

    # Check if atom is terminal or not
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
    
    # Consider interactions with terminal atom in one of 2 acids
    influence_all = out_details[(out_details[:,3].astype(float)>0) & np.bitwise_or(out_details[:,7]=='1', out_details[:,8]=='1')]

    # Group data by common 2 acids, calculate WeightSUM, percentMC
    out = []
    processed = []
    for res1, res2 in data[:, [0,1]]:
        if [res1, res2] not in processed:
            processed.append([res1, res2])
            grouped_data = data[(data[:, 0]==res1) & (data[:, 1]==res2)]
            edgeSum = sum(grouped_data[:, 3].astype('float')) # weight_sum
            chain = np.unique(grouped_data[:, 2], return_counts=True)
            chain = dict(zip(chain[0], chain[1]))
            # can we just use value_counts w/o calculating all numbers of sidechain
            MCMCcount=0; MCSCcount=0; SCSCcount=0; SCMCcount=0
            for tp in grouped_data[:, 2]:
                if tp=='MCMC' or tp=='MCSC':
                    MCMCcount+=1
                    MCSCcount+=1
                elif tp=='SCMC':
                    SCMCcount+=1
                elif tp=='SCSC':
                    SCSCcount+=1
            percentMC = (MCMCcount+MCSCcount) / (MCMCcount+MCSCcount+SCMCcount+SCSCcount)
            values = [grouped_data[0,0], grouped_data[0,1], 'mixed', edgeSum, percentMC, sum(grouped_data[:, 7].astype(int)),
                        sum(grouped_data[:, 8].astype(int)), len(data[data[:, 0] == grouped_data[0, 0]])]
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
        For each interaction in data, check if acids names are sorted or not
        If not then replace acids (acid1,acid2 --> acid2, acid1), also MCSC --> SCMC, num terminal atoms
        Select unique pairs, only directed interactions
        Select subdata grouped by 2 acids and calculate sumWeight*0.5 of all interactions (here we consider bonds twice, directed bonds)
    @input
        data - np.array of data of interactions 
    @output
        np.array in format [acid1, acid2, SumWeight]
    '''
    pairs = np.sort(data[:,:2], axis=1)
    pairs = np.unique(pairs, axis=0)
    out = []
    for i in range(len(pairs)):
        subdata = data[(data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1])]
        out.append([pairs[i,0], pairs[i,1], 0.5*sum(subdata[:,3].astype(float))])

    return np.array(out)

def make_undirected_interactions(data):
    '''
    @description
        Consider bidirected interactions: from input data sort acids by name and select unique pairs, (ALA0A GLY-1A)
        Select subdata grouped by 2 acids and calculate sumWeight
    @input
        data - np.array of data of interactions
    @output
        np.array in format [acid1, acid2, SumWeight]
    '''

    # dat = data.copy()
    pairs = np.sort(data[:,:2], axis=1)
    pairs = np.unique(pairs, axis=0)
    out = []
    for i in range(len(pairs)):
        subdat = data[((data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1])) | ((data[:,1]==pairs[i,0]) & (data[:,0]==pairs[i,1]))]
        out.append([subdat[0,0], subdat[0,1], sum(subdat[:,3].astype(float))])
    return np.array(out)

def make_directed_interactions(data):
    '''
    Consider directed interactions 
    same as above function
    '''
    pairs = np.sort(data[:,:2], axis=1)
    pairs = np.unique(data[:,:2], axis=0)
    out = []
    for i in range(len(pairs)):
        subdat = data[((data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1]))]
        out.append([subdat[0,0], subdat[0,1], sum(subdat[:,3].astype(float))])
    return np.array(out)


if __name__ == '__main__':
    # Read all needed files 
    # 1. terminalAtoms 2. SecondaryStructure 3. RSA 4. Net file with all interactions

    pdb = sys.argv[1] # input is pdb file
    tf = sys.argv[2] # time frame, 1btl_0.pdb
    try:        
        terminalAtoms = {}
        with open('../Test_data/terminalAtoms', 'r') as terminalAtomFile:
            for line in terminalAtomFile:
                line = line.strip("\n").split()
                terminalAtoms[line[0]] = [x.strip(',') for x in line[1:]]
                
        secStructure = np.loadtxt(re.sub('.pdb', '_secondaryStructure2', pdb), dtype=np.str)
        rsa = np.loadtxt(re.sub('.pdb', '.rsa', pdb), dtype=np.str)
        data_all = np.loadtxt(re.sub('.pdb', '_net', tf), dtype=np.str)
        # data_all = np.loadtxt(re.sub('.pdb', '_net', pdb), dtype=np.str)
    except FileNotFoundError:
        print('File not accesible')

    col0 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 0])))
    col1 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 1])))
    data_all[:, 0] = col0[:, 0]
    data_all[:, 1] = col1[:, 0]
    data_all = np.concatenate((data_all, col0[:, [1]]), axis=1)
    data_all = np.concatenate((data_all, col1[:, [1]]), axis=1)
    data_all

    # Define nodes that are interacting in net file
    nodes = np.unique(data_all[:,0:2])
    originalNodes = nodes

    # Process net file, group by common 2 interacted acids, calculate SumWeight, add terminal atoms
    basedata = process_data(data_all)
    basedata_noPP = removeRedundancy(basedata['out'])
    # Consider only sidechain interactions w/o MC/MC
    basedata_allsidechain = make_undirected_interactions(data_all[data_all[:,2]!='MCMC'])

    influencedata_all = basedata['influence_all']
    influencedata_collapse = process_data(influencedata_all)
    influencedata = removeRedundancy(influencedata_collapse['out'])


    influencenet = igraph.Graph().TupleList(edges=influencedata[:,:2], directed=False)
    influencenet_W = 1/influencedata[:,2].astype(float) # change 3 to 2

    # Define directed network
    influencedata_directed = make_directed_interactions(influencedata_all[influencedata_all[:,7]=='1'])
    influencenet_directed = igraph.Graph().TupleList(edges=influencedata_directed[:,[1,0]], directed=True)
    influencenet_directed_W = 1/influencedata_directed[:,2].astype(float) # change 3 to 2
 
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
    net_W = 1/basedata_noPP[:,2].astype(float) # change 3 to 2

    # Use Walktrap algo to define clusters for nodes and write it into net_community_vec dictionary
    net_community = net.community_walktrap(net_W, steps=4)
    cluster = net_community.as_clustering()

    net_community_vec = {}
    for name, member in zip(cluster.graph.vs['name'], cluster.membership):
        net_community_vec[name] = member

    net_community_vec_wt = net_community_vec
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

    firstOrderDegree, firstOrderDegree_sidechain, secondOrderDegree, secondOrderDegree_sidechain = first_second_order_degree(influencenet, influencenet_directed, nodes)


    # Use SecondaryStructure file to define clusters for nodes
    factors = sorted(np.unique(secStructure[:,1]))
    factors_dict = {factor:num for num, factor in enumerate(factors, 1)}
    net_community_vec = {node:factors_dict[f] for node, f in zip(secStructure[:,0], secStructure[:,1])}

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


    out.to_csv(re.sub('.pdb', '', tf)+'_scoresEnergetics', sep='\t', index=False)
    outZ.to_csv(re.sub('.pdb', '', tf)+'_scoresEnergeticsZ', sep='\t', index=False)
    # out.to_csv(re.sub('.pdb', '', pdb)+'_scoresEnergetics', sep='\t', index=False)
    # outZ.to_csv(re.sub('.pdb', '', pdb)+'_scoresEnergeticsZ', sep='\t', index=False)
