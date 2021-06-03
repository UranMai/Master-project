# before create frame_1000.pdb and its files
# run >> python centroid_np.py ../Test_data/frame_1000.pdb
# it reads frame_1000_net, frame_1000_secondaryStructure2, frame_1000.rsa, frame_1000_ligands, frame_1000_centroidNetLigand

import numpy as np
import pandas as pd 
import re, os
import sys
from functools import reduce 
import igraph
from sklearn.preprocessing import scale
from functools import reduce # merge dataframes
from score_functions import *

linearcutoff = 1 
sidechainMode = True ##build network of just SCSC, SCMC and MCSC. Still use unique atoms for the _uniqueAtoms directed plots
useDNA=True
removeWaters=True  
ligandCentroidMode = True
ligandmode = False

def add_terminal_atoms(data_net):
    '''
    @description
        Create data of reverse interactions, Merge directed and reversed data (bi-directed interactions)
        Example: 
                data_net:       ['GLU:22:A', 'ARG:231:A', 'SCSC', '10', 'SB', 'OE1', 'NH1']
                reverse_data:   ['ARG:231:A', 'GLU:22:A', 'SCSC', '10', 'SB', 'NH1', 'OE1']

        For merged data add info about whether atoms are terminal(1) or not terminal(0)
    @input
        data_net - np.array data of interactions ('net' file)
        data_net must contain info in format 
                [acid1, acid2, Chain types, Energy_value, Bond_type, atom1, atom2]
                ['GLU:22:A', 'ARG:231:A', 'SCSC', '10', 'SB', 'OE1', 'NH1']
    @output
        data - np.array of processed bidirected data
    '''
    # Edges should be listed bi-directionally
    if data_net.shape[1] == 7:
        reverse_data = data_net.copy()
        mask1 = np.where(reverse_data[:, 2]=='MCSC')
        mask2 = np.where(reverse_data[:, 2]=='SCMC')
        reverse_data[mask1, 2] = 'SCMC'
        reverse_data[mask2, 2] = 'MCSC'
        reverse_data[:,[0, 1, 5, 6]] = reverse_data[:, [1, 0, 6, 5]] 

    # data = np.concatenate((data_net, reverse_data), axis=0, )
    data = np.vstack((data_net, reverse_data))
    data = np.unique(data, axis=0) 

    # Check if acids' atoms are terminal or not
    codes1 = []; codes2 = []
    for bond in data:
        code1 = 0
        code2 = 0 
        if bond[5] in terminalAtoms[bond[0][:3]]:
            code1 = 1
        if bond[6] in terminalAtoms[bond[1][:3]]:
            code2 =1
        if (bond[4] == 'PICATION' or bond[4] == 'PIPI') or (bond[4] == 'DNA' and useDNA==True):
            code1 = 1
            code2 = 1
        codes1.append(code1)
        codes2.append(code2)
    # Add terminal atoms 
    data = np.insert(data, 7, codes1, axis=1)
    data = np.insert(data, 8, codes2, axis=1)
    return data

def make_undirected_interactions(data):
    '''
    @description
        Consider undirected interactions: from input data sort acids by name and select unique pairs,
        Select subdata grouped by unique pair of2 acids and calculate sumWeight
    @input
        data - np.array, data of interactions
    @output
        np.array in format [Acid1, Acid2, SumWeight]
    '''
    pairs = np.sort(data[:,:2], axis=1)
    pairs = np.unique(pairs, axis=0)
    output = []
    for i in range(len(pairs)):
        subdat = data[np.bitwise_and(data[:,0]==pairs[i,0], data[:,1]==pairs[i,1]) | np.bitwise_and(data[:,1]==pairs[i,0], data[:,0]==pairs[i,1])]
        out = []
        for dat in subdat:
            if np.argsort(dat[0:2])[0]==1:
                out.append([dat[1], dat[0], dat[3]])
            else:
                out.append([dat[0], dat[1], dat[3]])
        subdata = np.array(out)
        subdata = np.unique(subdata, axis=0)
        output.append([subdata[0,0], subdata[0,1], sum(subdata[:,2].astype(float))])
    # for i in range(len(pairs)):
    #     subdat = data[((data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1])) | ((data[:,1]==pairs[i,0]) & (data[:,0]==pairs[i,1]))]
    #     # subdat = data[np.bitwise_and(data[:,0]==pairs[i,0], data[:,1]==pairs[i,1]) | np.bitwise_and(data[:,1]==pairs[i,0], data[:,0]==pairs[i,1])]
    #     output.append([pairs[i,0], pairs[i,1], 0.5*sum(subdat[:,3].astype(float))])
    return np.array(output)

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
    pairs = np.unique(data[:,:2], axis=0)
    out = []
    for i in range(len(pairs)):
        subdat = data[(data[:,0]==pairs[i,0]) & (data[:,1]==pairs[i,1])]
        out.append([subdat[0,0], subdat[0,1], sum(subdat[:,3].astype(float))])
    return np.array(out)

def exclude_neirbor_nodes(data):
    # don't include interactions b/w neighboring acids
    out = []
    for row in data:
        condition = np.abs(int(re.sub(r'[^-\d]', '',row[0])) - int(re.sub(r'[^-\d]', '',row[1])))
        if condition > linearcutoff:
            out.append(row)
    return np.array(out)

if __name__ == '__main__':
    pdb = sys.argv[1] # input pdb file
    try:        
        data_all = np.loadtxt(pdb.replace('.pdb', '_net'), dtype=np.str)
        secStructure = np.loadtxt(pdb.replace('.pdb', '_secondaryStructure2'), dtype=np.str)
        rsa = np.loadtxt(pdb.replace('.pdb', '.rsa'), dtype=np.str)
        # bfactor = np.loadtxt(pdb.replace('.pdb','_Bfactor'), dtype=np.str)
        liganddat = np.loadtxt(pdb.replace('.pdb','_ligands'), dtype=np.str)
        if len(liganddat) != 0:
            ligandmode = True
        ligandCentroidFile = np.loadtxt(pdb.replace('.pdb', '_centroidNetLigand'), dtype=np.str)
        
        terminalAtoms = uniqueAtoms
    except FileNotFoundError:
        print('File not found')
    
    col0 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 0])))
    col1 = np.array(list(map(lambda x: x.rsplit(':', 1), data_all[:, 1])))
    data_all[:, 0] = col0[:, 0]
    data_all[:, 1] = col1[:, 0]
    data_all = np.concatenate((data_all, col0[:, [1]]), axis=1)
    data_all = np.concatenate((data_all, col1[:, [1]]), axis=1)
    data_all = data_all[data_all[:,3].astype(float)>0]
    # data_all must be in format ['LYS:6:A' 'ASP:9:A' 'SCSC' '10' 'SB' 'NZ' 'OD2']
    # where columns 0,1-aminoacids; 2-chain type; 3-energy value; 4-bond type; 5,6-atoms

    # Define phi angles from Secondary Structure
    phiAngleDB={acid : angle for acid, angle in secStructure[:,[0,2]]}

    # Fill rsa, secStructure with ligands names; if ligandmode is False set empty list of ligands 
    if ligandmode:
        ligands = np.unique(liganddat[:,1])
        for ligand in ligands:
            rsa = np.append(rsa, [[':'.join(ligand.split(':')), '-1']], axis=0)
            secStructure = np.append(secStructure, [[':'.join(ligand.split(':')), 'ligand', '360']], axis=0)
            phiAngleDB[ligand] = 0
    else:
        ligands = [] # no ligands

    # skip HOH check
    # if removeWaters:
        # delete HOH waters

    # Define nodes (aminoacids) that are interacting in net file
    nodes = np.unique(data_all[:,0:2])
    originalNodes = nodes

    # Process data of interactions
    basedata = add_terminal_atoms(data_all)
    basedata_noPP = make_undirected_interactions(basedata[basedata[:,4]!='PP']) # Consider only undirected interactions excluding PP bonds 

    if sidechainMode: # Consider only sidechain bonds
        influencedata_detailed = basedata[basedata[:,2]!='MCMC']
        influencedata = make_undirected_interactions(influencedata_detailed) 
        influencedata = exclude_neirbor_nodes(influencedata)
        influencedata = influencedata[influencedata[:,2].astype(float)>0] 

        influencedata_directed = make_directed_interactions(influencedata_detailed[influencedata_detailed[:,7]=='1'][:,:4])
        influencedata_directed = exclude_neirbor_nodes(influencedata_directed)

        influencenet = igraph.Graph.TupleList(edges=influencedata[:,:2], directed=False)
        influencenet_weight = 1/influencedata[:,2].astype(float) 

        influencenet_directed = igraph.Graph.TupleList(edges=influencedata_directed[:,:2], directed=True)
        influencenet_directed_weight = 1/influencedata_directed[:,2].astype(float) 

    #Replace missing nodes in the igraph net
    nodesMissing = [node for node in nodes if node not in np.unique(influencedata[:,0:2])]

    if any(np.isin(nodes, secStructure[:,0])) == False:
        for node in nodes:
            if (node in secStructure[:,0])==False:
                secStructure = np.append(secStructure, [[node, 'xxx', '360']], axis=0)
            if  (node in rsa[:,0]) == False:
                rsa = np.append(rsa, [[node, 0]], axis=0)

    # CENTRALITY CALCULATION
    # Compute vertex betweenness for undirected and directed graph with weights
    node_betweenness = influencenet.betweenness(weights=influencenet_weight, directed=False) 
    node_betweenness_unique = influencenet_directed.betweenness(weights=influencenet_directed_weight, directed=False)

    net = igraph.Graph.TupleList(edges=basedata_noPP[:,[0,1]], directed=False)
    net_weight = 1/basedata_noPP[:,2].astype(float)

    # -------------------------WALKTRAP CLUSTERING-------------------------------------------------
    # Use Walktrap algo to define clusters for nodes and write it into net_community_vec dictionary
    net_community = net.community_walktrap(weights=net_weight, steps=4)
    cluster = net_community.as_clustering()
    net_community_vec = {node:cluster for node, cluster in zip(cluster.graph.vs['name'], cluster.membership)}
    nodes = influencenet.vs['name']

    # Compute edge betweenness
    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_weight)

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
                                                                                        method='centroid')
    firstOrderDegree, firstOrderDegree_sidechain, \
    secondOrderDegree, secondOrderDegree_sidechain = first_second_order_degree(influencenet, influencenet_directed)


    # ------------------------SecondaryStructure CLUSTERING---------------------------------------
    # Use SecondaryStructure file to define clusters for nodes
    factors = sorted(np.unique(secStructure[:,1]))
    factors_dict = {factor:num for num, factor in enumerate(factors, 1)}
    net_community_vec = {node:factors_dict[f] for node, f in zip(secStructure[:,0], secStructure[:,1])}

    nodes = influencenet.vs['name']

    node_edge_betweenness_stride = node_edge_btwns(influencedata, edge_betweenness, nodes, net_community_vec)
    node_edge_betweenness_sidechain_stride = node_edge_btwns_SC(influencedata_directed, influencedata, edge_betweenness, nodes, net_community_vec)
    node_intermodular_degree_stride, node_modules = node_intermodular_degree(influencedata, nodes, net_community_vec)
    node_intermodular_degree_sidechain_stride, node_modules_sidechain = node_intermodular_degree_SC(influencedata_directed, nodes, net_community_vec)

    secondOrder_node_intermodular_degree_stride, \
    secondOrder_node_intermodular_degree_sidechain_stride =node_intermodular_degr_2order(node_intermodular_degree_stride, 
                                                                                            node_intermodular_degree_stride, 
                                                                                            node_modules,
                                                                                            node_modules_sidechain, 
                                                                                            method='centroid')
    # Set nodes not including ligands
    nodes = [node for node in nodes if node not in ligands]

    # Create dataframe from dictionaries
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
        #    rsa_df,
            firstOrderDegree_df, firstOrderDegree_sidechain_df, 
            secondOrderDegree_df, secondOrderDegree_sidechain_df,
            node_edge_betweenness_stride_df, node_edge_betweenness_sidechain_stride_df,
            node_intermodular_degree_stride_df, node_intermodular_degree_sidechain_stride_df,
            secondOrder_node_intermodular_degree_stride_df, secondOrder_node_intermodular_degree_sidechain_stride_df,
            node_edge_betweenness_wt_df, node_edge_betweenness_sidechain_wt_df,
            node_intermodular_degree_wt_df, node_intermodular_degree_sidechain_wt_df,
            secondOrder_node_intermodular_degree_wt_df, secondOrder_node_intermodular_degree_sidechain_wt_df
            ]

    final_df = reduce(lambda left,right: pd.merge(left,right,on=0), out)

    colnames = ['AA', 'Degree', 'Degree_uniqueAtoms', 'SecondOrderDegree', 'SecondOrderDegree_uniqueAtoms',
                'NodeEdgeBetweennessSTRIDE', 'NodeEdgeBetweennessSTRIDE_uniqueAtoms', 'IntermodularDegreeSTRIDE',
                'IntermodularDegreeSTRIDE_uniqueAtoms', 'SecondOrderIntermodularDegreeSTRIDE',
                'SecondOrderIntermodularDegreeSTRIDE_uniqueAtoms', 'NodeEdgeBetweennessWALKTRAP',
                'NodeEdgeBetweennessWALKTRAP_uniqueAtoms', 'IntermodularDegreeWALKTRAP',
                'IntermodularDegreeWALKTRAP_uniqueAtoms', 'SecondOrderIntermodularDegreeWALKTRAP',
                'SecondOrderIntermodularDegreeWALKTRAP_uniqueAtoms']
    final_df.columns = colnames

    # consdier only terminal atoms of acids
    # groupby acid name and find ligand which has shortest distance to acid
    if ligandmode==True and ligandCentroidMode==False:
        # liganddat = np.loadtxt(re.sub('.pdb','_ligand',pdb), dtype=np.str)
        # liganddat was already read
        tmp = []
        for row in liganddat:
            acid = row[0].split('-') #return ASP # ASP6A-CA #need approprieate format
            if acid[1] in terminalAtoms[acid[0][:3]]: #if CA in ASP
                tmp.append(row)
        ligandvec = pd.DataFrame(tmp).groupby(0)[2].min().to_dict()
        

    if ligandmode == True and ligandCentroidMode == True:
        liganddat = np.loadtxt(re.sub('.pdb','',pdb)+'_centroidNetLigand', dtype=np.str)
        tmp = []
        for row in liganddat:
            tmp.append(row)
        ligandvec = pd.DataFrame(tmp)
        ligandvec[2] = ligandvec[2].astype(float)
        ligandvec = ligandvec.groupby(0)[2].min().to_dict()

    if ligandmode == False:
        ligandvec = {node:0 for node in nodes}

    # add nodes that were not networked in influencenet graph
    zeroNodes = []
    for node in originalNodes:
        if node not in nodes:
            zeroNodes.append(node)

    zeroMat = pd.DataFrame(0, index=np.arange(len(zeroNodes)), columns=final_df.columns[0:])
    zeroMat['AA'] = zeroNodes

    output_df = pd.concat([final_df, zeroMat])

    # Extract rsa info only for nodes that are in final_df(out['AA'])
    rsa_list = []
    for acid in output_df['AA'].to_list():
        if acid in rsa[:,0]:
            rsa_list.append(rsa[rsa[:,0]==acid][0])
        # add else what if no rsa info for existing node
    rsa_df = pd.DataFrame(rsa_list, columns=['AA', 'RSA'])
    rsa_df['RSA'] = rsa_df['RSA'].astype(float)
    
    ligand_list = []
    for acid in output_df['AA'].to_list():
        if acid in ligandvec:
            ligand_list.append([acid, float(ligandvec[acid])])
        else:
            ligand_list.append([acid, 0])
    ligand_df = pd.DataFrame(ligand_list, columns=['AA', 'Ligand'])

    # Add ligand and rsa df
    output_df = pd.merge(output_df, ligand_df, on='AA')
    output_df = pd.merge(output_df, rsa_df, on='AA')
    #Add 1 to everything to avoid zeros, except ligand
    output_df = output_df.fillna(value=0)
    output_df.loc[output_df['Ligand'] == 0, 'Ligand'] = 150 # set for Nan ligand 150 value

    out_total = output_df[output_df.columns[1:-1]]+1

    final_out = pd.concat([output_df[['AA', 'RSA']], out_total], axis=1)

    # for comparing resutls with R code
    for col in final_out.columns[1:]:
        print(col, final_out[col].sum())

    outZ = Scaling(final_out[final_out.columns[1:]])
    outZ = pd.concat([final_out[['AA']], outZ], axis=1)

    final_out.to_csv(re.sub('.pdb', '', pdb)+'_scoresCentroid', sep='\t', index=False)
    outZ.to_csv(re.sub('.pdb', '', pdb)+'_scoresCentroidZ', sep='\t', index=False)