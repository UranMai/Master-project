import pandas as pd
import numpy as np
import igraph
import sys, re
import time
t0= time.time()
pdb = sys.argv[1]
pdb_file = re.sub(".pdb","",pdb)

terminalAtoms = {}
with open('terminalAtoms', 'r') as terminalAtomFile:
    for line in terminalAtomFile:
        line = line.strip("\n").split()
        terminalAtoms[line[0]] = [x.strip(',') for x in line[1:]]

secStructure = pd.read_csv(pdb_file+'_secondaryStructure', sep='\t', header=None)
rsa = pd.read_csv(pdb_file+'.rsa', sep='\t', header=None)


def COLLAPSE(data):
    '''
    @descirption
    For observing reverse interactions, create data where sidechains are replaced, and concat data
    Add 2 columns to data, whether acids' atom in terminalAtoms (1) or not (0)
    Group data by two aminoacids and iterate through them, counting sidechain types
    and calculate percent of MC (main chains b/w acids)
    Append info about 2 acids, Weightsum, MC percent, amount of terminal atoms, degree of acid1 in whole data->create dataframe 'out_data'
    out_data_noPP w/o polypeptide pairs, and take only with SumWeightEdges > 0
    Create df "influence_all" with terminal atoms in one of the acids
    and create df "out_sub" with non terminal atoms in both acids

    @input
    data - file_net with all interactions in a protein

    @output
    dict of dataframes, where 
        out_data - grouped data by acids with sumWeight, MCpercent, degree of first acid
        out_data_noPP - out_data w/o PP type (polypeptide pairs)
        influence_all - data with one of terminal atoms in acids
        out_subtract - changed data with negative Weight where none terminal atoms in two acids
    '''
    if data.shape[1] == 7:
        replaced_data = data.replace({'MCSC':'SCMC', 'SCMC':'MCSC'})
        replaced_data = replaced_data[['AA2', 'AA1', 'Type1', 'Weight', 'Type2', 'Atom2', 'Atom1']]
        replaced_data = replaced_data.rename(columns={'AA2':'AA1', 'AA1':'AA2', 'Atom2':'Atom1', 'Atom1':'Atom2'})
    else:
        replaced_data = data.replace({'MCSC':'SCMC', 'SCMC':'MCSC'})
        replaced_data = replaced_data[['AA2', 'AA1', 'Type1', 'Weight', 'Type2', 'Atom2', 'Atom1', 'influence2', 'influence1']]
        replaced_data = replaced_data.rename(columns={
            'AA2':'AA1', 'AA1':'AA2', 'Atom2':'Atom1', 'Atom1':'Atom2', 'influence1':'influence2', 'influence2':'influence1'})

    data = pd.concat([data, replaced_data], axis=0, ignore_index=True)
    data = data.drop_duplicates()

    if data.shape[1] == 7:
        codes1 = [];    codes2 = []
        for res, atom in zip(data['AA1'], data['Atom1']):
            codes1.append(1 if atom in terminalAtoms[res[:3]] else 0)
        for res, atom in zip(data['AA2'], data['Atom2']):
            codes2.append(1 if atom in terminalAtoms[res[:3]] else 0)
        data['influence1'] = codes1
        data['influence2'] = codes2
        out_details = data
    else:
        out_details = data

    subdata = data.groupby(by=['AA1', 'AA2'], sort=False)
    out = []; out_noPP = []; completed = []
    for key, item in subdata:
        grouped_df = subdata.get_group(key)
        edgeSum = sum(grouped_df['Weight']) # summ weight for certain 2 aminoacids
        # info = grouped_df.to_numpy()[0] # ddelete
        chain = grouped_df['Type1'].value_counts().to_dict() # calculate sidechain types
        MCMCcount = 0;    MCSCcount = 0
        SCMCcount = 0;    SCSCcount = 0
        # here need to look at counting of MCMC and MCSC, b/c their number will be equal 
        for row in grouped_df.values:
            if 'MCMC' in row[2] or 'MCSC' in row[2]:
                MCMCcount += 1
                MCSCcount += 1
            elif 'SCMC' in row[2]:
                SCMCcount += 1
            elif 'SCSC' in row[2]:
                SCSCcount += 1
        # count Main Chains b/w 2 aminoacids
        percentMC = (MCMCcount+MCSCcount) / (MCMCcount+MCSCcount+SCMCcount+SCSCcount)
        # add to list 
        out.append([info[0], info[1], 'mixed', edgeSum, percentMC, sum(grouped_df['influence1']), sum(grouped_df['influence2']), sum(data['AA1']==''.join(grouped_df['AA1'].unique()))])
        # create other list w/o polypeptide pairs; if all bonds b/w acids are PP dont add
        if all(grouped_df['Type2']!= "PP"):
                out_noPP.append([info[0], info[1], 'mixed', edgeSum, percentMC, sum(grouped_df['influence1']), sum(grouped_df['influence2']), sum(data['AA1']==''.join(grouped_df['AA1'].unique()))])    
        completed.append(grouped_df.index.to_list())

    out_data = pd.DataFrame(out, columns=["AA1", "AA2", "Type", "SumEdges", "percentMC", "influence1", "influence2","degree"])
    out_noPP_data = pd.DataFrame(out_noPP,  columns=["AA1", "AA2", "Type", "SumEdges", "percentMC", "influence1", "influence2","degree"])
    out_data = out_data[out_data['SumEdges'] > 0]
    out_noPP_data = out_noPP_data[out_noPP_data['SumEdges'] > 0]

    # Take only with positive weights and if one of atoms is terminal atom
    influence_all = out_details[(out_details['Weight']>0) & (out_details['influence1']==1) | (out_details['influence2']==1)]
    out_details = out_details[out_details['Weight']>0]
    out_subtract = out_details
    out_sub = out_subtract
    # for non terminal atoms in both acids with influence 0, change Weight value to negative
    out_sub.loc[(out_sub.influence1 == 0) & (out_sub.influence2 == 0), 'Weight'] = -out_sub['Weight']
    out_sub = out_sub[out_sub['Type2']!='PP']
    
    results = {}
    results['out'] = out_data
    results['out_noPP'] = out_noPP_data
    results['influence_all'] = influence_all
    results['out_subtract'] = out_sub
    return results

def removeRedundancy(out): 
    '''
    @input
    list of list or dataframe.value.tolist(); 'out' data

    @desciption
    for each list of acids check whether acids are sorted ny name? if not replace
    Create list of unique aminoacid pairs, and for each of them extract SumWeightEdges
    
    @output
    list data of 2 aminoacids and sumWeight
    '''
    dataTmp = [item[0:7] for item in out]
    df = []
    for row in dataTmp:
        # if argument sorting return 1 at 1st position
        if (np.argsort(row[0:2])[0]==1):
            row[0:2] = row[0:2][::-1]
            row[2] = 'xemi'
            row[5:7] = row[5:7][::-1]
            df.append(row)
        else:
            df.append(row)

    pairs = [] # list of acid pairs
    for row in df:
        if row[0:2] not in pairs:
            pairs.append(row[0:2])

    subData = []
    out = []
    for pair in pairs:
        for i in range(len(df)):
            if df[i][0] == pair[0] and df[i][1] == pair[1]:
                if (pair[0], pair[1], 'mixed', df[i][3]) not in out:
                    out.append([pair[0], pair[1], 'mixed', round(df[i][3],8)])
    unique_out = []
    for x in out: 
        if x not in unique_out:
            unique_out.append(x)
    return unique_out


def collapse_agnostic(data):
    '''
    @input
    dataframe with allsidechain types

    @description
    Sort acids by name and create new_dataframe
    Group new dataframe by acids, and summarize weigths, and return this data
    '''
    new_dataframe = []
    for row in data.values.tolist(): #data_all
        row[0:2] = sorted(row[0:2])
        new_dataframe.append(row)
    data = pd.DataFrame(new_dataframe, columns=['AA1', 'AA2', 'Weight']).groupby(by=['AA1', 'AA2']).sum()
    out = []
    for aminos, value in zip(data.index, data.values):
        out.append([aminos[0], aminos[1], 'mixed', float(value)])
    return out

def collapse_directed(col_dir):
    '''
    @input 
    data with terminal atom in acid1 with weights

    @description
    Group data by acids and summarize weigths for same acids
    and return list of this data
    '''
    x = col_dir.groupby(['AA1', 'AA2'], sort=False).sum()
    out = []
    for aminos, value in zip(x.index, x.values):
        out.append([aminos[0], aminos[1], 'mixed', float(value)])
    return out

print('Define functions', time.time()-t0)


#----------------------------------------------------------------------------------------------------
# read 'pdb_net' file, 
# where AA1, AA2 - aminoacids, 
# Type1 - sidechains type (MCMC); Type2 - bond type
# Weight - energy of interaction, 
# Atom1, Atom2 - involved atoms
data_all = pd.read_csv(pdb_file+'_net', header=None, sep='\t', names=['AA1','AA2','Type1', 'Weight','Type2','Atom1','Atom2']) # 3621x7
nodes = (data_all['AA1'].append(data_all['AA2'])).unique(); sets = nodes
originalNodes = nodes
basedata = COLLAPSE(data_all)
basedata_noPP = removeRedundancy(basedata['out'].values.tolist())

# take only sidechain data with aminoacids and weight
basedata_allsidechain = collapse_agnostic(data_all[data_all['Type1'] != 'MCMC'][['AA1', 'AA2','Weight']])

influencedata_all = basedata['influence_all']
influencedata_collapse = COLLAPSE(influencedata_all)
influencedata = removeRedundancy(influencedata_collapse['out'].values.tolist()) # acids with summ weights


# from data with teminal atoms in one of the acids extract aminoacids and create GRAPH NETWORK NONDIRECTED
df = pd.DataFrame(influencedata, columns=['AA1', 'AA2', 'type', 'Weight'])
influencenet = igraph.Graph().TupleList(edges=df[['AA1', 'AA2']].values.tolist(), directed=False)
influencenet_weight = [1/item[3] for item in influencedata]

# Take data with terminal atom in acid1 and create DIRECTED GRAPH NETWORK
influencedata_directed = collapse_directed(influencedata_all[influencedata_all['influence1']==1][['AA1', 'AA2', 'Type1','Weight']])
influencenet_directed = igraph.Graph().TupleList(edges=pd.DataFrame(influencedata_directed)[[1,0]].values.tolist(), directed=True)
influencenet_directed_weight = [1/item[3] for item in influencedata_directed]

influencedata_directed = pd.DataFrame(influencedata_directed)
influencedata = pd.DataFrame(influencedata)
# pd.DataFrame(influencedata).to_csv(pdb_file+'_allEdges', sep='\t', header=False, index=False)

# Read rsa file 
rsa = pd.read_csv(pdb_file+'.rsa',header=None, sep='\t')
# if node not in nodes of secstructure or rsa, add this node
sec_list = []
for node in nodes:
    if node not in secStructure[0].unique():
        dict1 = {0:node, 1:'xxx', 2:360}
        sec_list.append(dict1)
        # secStructure.append(dict1, ignore_index=True)
secStructure = secStructure.append(pd.DataFrame(sec_list))      

rsa_list = []
for node in nodes:
    if node not in rsa[0].unique():
        dict1 = {0:node, 1:0}
        rsa_list.append(dict1)
        # rsa.append(dict1, ignore_index=True)
rsa = rsa.append(pd.DataFrame(rsa_list))
print('Define all data', time.time()-t0)

# Begin centrality calculations
# Modular calculations
# Walktrap algorithm for clusters detection

net = igraph.Graph().TupleList(edges=pd.DataFrame(basedata_noPP)[[0,1]].values.tolist(), directed=False)
net_weight = [1/item[3] for item in basedata_noPP]

net_community=net.community_walktrap(net_weight, steps=4)
cluster = net_community.as_clustering()

# define dict of node and their cluster number
net_community_vec = {}
for name, member in zip(cluster.graph.vs['name'], cluster.membership):
    net_community_vec[name] = member

net_community_vec_wt = net_community_vec
nodes = net.vs['name']

def calculation(net_community_vec):
    '''
    @input 
    dict of {nodes:belonging cluster number}

    @description
    calculate attributes of network centrality calculations
    '''

    '''
    edge_betweenness values, create dict of nodes and their betweenness values
    iterate through line of influencedata
    if two acids of one connection lie in different clustes, add to created dict edge_betweenness values of node
    '''
    edge_betweenness = influencenet.edge_betweenness(weights=influencenet_weight)
    node_edge_betweenness = {}
    for node in nodes:
        node_edge_betweenness[node] = 0
    for i in range(len(edge_betweenness)):
        if net_community_vec[influencedata[0][i]] != net_community_vec[influencedata[1][i]]:
            node_edge_betweenness[influencedata[0][i]] = node_edge_betweenness[influencedata[0][i]] + edge_betweenness[i]
            node_edge_betweenness[influencedata[1][i]] = node_edge_betweenness[influencedata[1][i]] + edge_betweenness[i]
    '''
    node_intermodular degree - it's number of connecntions for nodes that are connected to nodes from other clusters
    
    iterate through nodes, for each current node select subdata where this current node is involved
    select from this subdata nodes that are connected to current node, and write it to bound
    fot these nodes compare their clusters with current node's cluster
    if they are in diff clusters, define "bound_modules" which is a dict of nodes and their cluster numbers
    add bound_modules{node:cluster number} to node_modules dict and 
    add len(bound_modules), which means degree of node connected to other cluster's nodes
    '''
    node_intermodular_degree = {}
    for node in nodes:
        node_intermodular_degree[node] = 0

    node_modules = dict()
    for node in nodes:
        subdata = influencedata[(influencedata[0] == node) | (influencedata[1] == node)]
        # subdata = pd.DataFrame(influencedata, columns=['A1', 'A2', 'type', 'W'])
        # subdata = subdata[(subdata['A1'] == node) | (subdata['A2'] == node)]
        bound = list(set(subdata[0].unique().tolist() + subdata[1].unique().tolist()))
        bound = [b for b in bound if b!=node]
        bound_modules = {}
        for b in bound:
            if net_community_vec[node] != net_community_vec[b]:
                bound_modules[b] = net_community_vec[b]

        node_modules[node] = bound_modules
        node_intermodular_degree[node] = len(bound_modules)

    '''
    node_edge_betweenness_sidechain dict for nodes

    For two nodes from different clusters 
    from influencedata_directed (which with terminalAtom in acid1) SELECT 2 subdata with right and reverse connections of nodes
    Define weights of these bonds, else define zero
    Add to node_edge_betweennes_sidechain dict weighted edge_betweenness values
    '''
    # influencedata_directed = pd.DataFrame(influencedata_directed)
    # influencedata = pd.DataFrame(influencedata) 
    node_edge_betweenness_sidechain = {}
    for node in nodes:
        node_edge_betweenness_sidechain[node] = 0 

    for i in range(len(edge_betweenness)):
        if net_community_vec[influencedata[0][i]] != net_community_vec[influencedata[1][i]]:
            data = influencedata_directed[(influencedata_directed[0] == influencedata[0][i]) & (influencedata_directed[1] == influencedata[1][i])]
            if len(data) == 0: 
                weight1 = 0 
            else: 
                weight1 = data.values[0][3]

            data = influencedata_directed[(influencedata_directed[0] == influencedata[1][i]) & (influencedata_directed[1] == influencedata[0][i])]
            if len(data) == 0: 
                weight2 = 0 
            else: 
                weight2 = data.values[0][3]

            if weight1 == 0 and weight2 == 0:
                node_edge_betweenness_sidechain[influencedata[0][i]] = node_edge_betweenness_sidechain[influencedata[0][i]] + 0
                node_edge_betweenness_sidechain[influencedata[1][i]] = node_edge_betweenness_sidechain[influencedata[1][i]] + 0
            else:
                node_edge_betweenness_sidechain[influencedata[0][i]] = node_edge_betweenness_sidechain[influencedata[0][i]] + (weight1/(weight1+weight2))*edge_betweenness[i]
                node_edge_betweenness_sidechain[influencedata[1][i]] = node_edge_betweenness_sidechain[influencedata[1][i]] + (weight2/(weight1+weight2))*edge_betweenness[i]
    

    '''
    Same as node_intermodular_degree, but here for another directed data with terminal atoms in acid1
    '''
    node_intermodular_degree_sidechain = {}
    for node in nodes:
        node_intermodular_degree_sidechain[node] = 0

    node_modules_sidechain = {}
    for node in nodes:
        subdata_sidechain = influencedata_directed
        # subdata_sidechain = pd.DataFrame(influencedata_directed.values, columns=['A1', 'A2', 'type', 'W'])
        subdata = subdata_sidechain[subdata_sidechain[0] == node]
        if len(subdata_sidechain) == 0:
            node_intermodular_degree_sidechain[node] = 0
            
        bound = list(set(subdata[0].unique().tolist() + subdata[1].unique().tolist()))
        bound = [b for b in bound if b!=node]
        bound_modules = {}
        for b in bound:
            if net_community_vec[node] != net_community_vec[b]:
                bound_modules[b] = net_community_vec[b]
        node_modules_sidechain[node] = bound_modules
        node_intermodular_degree_sidechain[node] = len(bound_modules)
    
    '''
    For current node that have connections with nodes from other clusters 
    find these nodes' connections with other nodes from other clusters, and
    add to dicts number of secondary connections minus first connection 
    '''
    secondOrder_node_intermodular_degree = {}
    for node, value in node_intermodular_degree.items():
        secondOrder_node_intermodular_degree[node] = 0

    secondOrder_node_intermodular_degree_sidechain = {}
    for node, value in node_intermodular_degree_sidechain.items():
        secondOrder_node_intermodular_degree_sidechain[node] = 0

    for node in node_intermodular_degree:
        if node_intermodular_degree[node] == 0:
            continue
        else: 
            tmp = []
            for nm in node_modules[node]:
                for nm1 in node_modules[nm]:
                    tmp.append(nm1)
            secondOrder_node_intermodular_degree[node] = len(tmp) - 1

            tmp1 = []
            for nm in node_modules_sidechain[node]:
                for nm1 in node_modules_sidechain[nm]:
                    tmp1.append(nm1)
            secondOrder_node_intermodular_degree_sidechain[node] = len(tmp1) - 1

    return [node_edge_betweenness, node_edge_betweenness_sidechain, 
            node_intermodular_degree, node_intermodular_degree_sidechain, 
            secondOrder_node_intermodular_degree, secondOrder_node_intermodular_degree_sidechain]

'''
Define degree of nodes(vertex) as firstOrderDegree
for SecondOrderDegree define node's neighbors and for them find number of connections, excluding node
the len of connections will be SecondOrderDegree
'''
firstOrderDegree = {}
for name, deg in zip(influencenet.vs['name'], influencenet.degree()):
    firstOrderDegree[name] = deg

firstOrderDegree_sidechain = {}
for name, deg in zip(influencenet_directed.vs['name'], influencenet_directed.degree(mode='in')):
    firstOrderDegree_sidechain[name] = deg

secondOrderDegree = {}
secondOrderDegree_sidechain = {}
nodes = influencenet.vs['name']

for num, node in enumerate(nodes):
    firstorder = influencenet.neighbors(node)
    secondorder = []
    for neigh in firstorder:
        secondorder.append(influencenet.neighbors(neigh))
    secondorder = list(set([i for sub in secondorder for i in sub if i!=num]))
    secondOrderDegree[node] = len(secondorder)
    firstorder = influencenet_directed.neighbors(node, mode='in')
    secondorder = []
    for neigh in firstorder:
        secondorder.append(influencenet_directed.neighbors(neigh, mode='in'))
    num1 = influencenet_directed.vs['name'].index(node)
    secondorder = list(set([i for sub in secondorder for i in sub if i!=num1]))
    secondOrderDegree_sidechain[node] = len(secondorder)

result = calculation(net_community_vec)
node_edge_betweenness_wt = result[0]
node_edge_betweenness_sidechain_wt = result[1]
node_intermodular_degree_wt = result[2]
node_intermodular_degree_sidechain_wt = result[3]
secondOrder_node_intermodular_degree_wt = result[4]
secondOrder_node_intermodular_degree_sidechain_wt = result[5]
print('Calculate features for 1st net_community_vec', time.time()-t0)

# Here we define new clusters from secondaryStructure file
# and for new netwrok cluster calculate centrality attributes
secSturcture = pd.read_csv(pdb_file+'_secondaryStructure', sep='\t', header=None)
factors = sorted(secSturcture[1].unique(), key=lambda v: v.upper())
types_dict = {}
for num, factor in enumerate(factors, 1):
    types_dict[factor] = num

net_community_vec = {}
for i, j in zip(secSturcture[0].values, secSturcture[1].values):
    net_community_vec[i] = types_dict[j]

result = calculation(net_community_vec)

node_edge_betweenness_stride = result[0]
node_edge_betweenness_sidechain_stride = result[1]
node_intermodular_degree_stride = result[2]
node_intermodular_degree_sidechain_stride = result[3]
secondOrder_node_intermodular_degree_stride = result[4]
secondOrder_node_intermodular_degree_sidechain_stride = result[5]
print('Calculate features for 2nd net_community_vec', time.time()-t0)


def Nanfill_df(dictionary):
    '''
    if node is not contained in data; add this node with None value
    '''
    df = pd.DataFrame.from_dict(data=dictionary.items())
    for node in nodes:
        if node not in df[0].values:
            dict1 = {0:node, 1:None}
            df = df.append(dict1, ignore_index=True)
    return df

firstOrderDegree_df = Nanfill_df(firstOrderDegree)
firstOrderDegree_sidechain_df = Nanfill_df(firstOrderDegree_sidechain)

secondOrderDegree_df = Nanfill_df(secondOrderDegree)
secondOrderDegree_sidechain_df = Nanfill_df(secondOrderDegree_sidechain)

node_edge_betweenness_stride_df = Nanfill_df(node_edge_betweenness_stride)
node_edge_betweenness_sidechain_stride_df = Nanfill_df(node_edge_betweenness_sidechain_stride)

node_intermodular_degree_stride_df = Nanfill_df(node_intermodular_degree_stride)
node_intermodular_degree_sidechain_stride_df =Nanfill_df(node_intermodular_degree_sidechain_stride)

secondOrder_node_intermodular_degree_sidechain_stride_df = Nanfill_df(secondOrder_node_intermodular_degree_sidechain_stride)
secondOrder_node_intermodular_degree_stride_df = Nanfill_df(secondOrder_node_intermodular_degree_stride)

node_edge_betweenness_wt_df = Nanfill_df(node_edge_betweenness_wt)
node_edge_betweenness_sidechain_wt_df = Nanfill_df(node_edge_betweenness_sidechain_wt)

node_intermodular_degree_wt_df = Nanfill_df(node_intermodular_degree_wt)
node_intermodular_degree_sidechain_wt_df = Nanfill_df(node_intermodular_degree_sidechain_wt)

secondOrder_node_intermodular_degree_wt_df = Nanfill_df(secondOrder_node_intermodular_degree_wt)
secondOrder_node_intermodular_degree_sidechain_wt_df = Nanfill_df(secondOrder_node_intermodular_degree_sidechain_wt) #???

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
print('Merge all data into one', time.time()-t0)



from functools import reduce # concatenate all data

df_final = reduce(lambda left,right: pd.merge(left,right,on=0), out)

colnames = ['AA', "Degree","Degree_sidechain","SecondOrderDegree","SecondOrderDegree_sidechain",
  "NodeEdgeBetweennessSTRIDE","NodeEdgeBetweennessSTRIDE_sidechain","IntermodularDegreeSTRIDE",
  "IntermodularDegreeSTRIDE_sidechain","SecondOrderIntermodularDegreeSTRIDE",
  "SecondOrderIntermodularDegreeSTRIDE_sidechain","NodeEdgeBetweennessWALKTRAP",
  "NodeEdgeBetweennessWALKTRAP_sidechain","IntermodularDegreeWALKTRAP",
  "IntermodularDegreeWALKTRAP_sidechain","SecondOrderIntermodularDegreeWALKTRAP",
  "SecondOrderIntermodularDegreeWALKTRAP_sidechain"]

# Name columns
df_final.columns = colnames

# OriginalNodes - nodes in pdb_net file; nodes - from influencenet graph
# if original nodes not in nodes, add to list of zeronodes
# and create dataframe from these zeronodes with same columns
# make final output (out11)
zeroNodes = []
for i in originalNodes:
    if i not in nodes:
        zeroNodes.append(i) 

zeroMat = pd.DataFrame(0, index=np.arange(len(zeroNodes)), columns=df_final.columns[1:])
zeroMat['AA'] = zeroNodes

out11 = pd.concat([df_final, zeroMat])

# fill None values with zero and add 1 to all rows
out = out11.fillna(value=0)
out3 = out[out.columns[1:]]+1
final_out = pd.concat([out[['AA']], out3], axis=1)

# normalization of data
def Scaling(df):
    cols = df.columns
    for col in cols:
        std = np.std(df[col], ddof=1)
        mean = np.mean(df[col])
        df[col] = (df[col] - mean)/std
    return df

out = final_out[final_out.columns[1:]]
outZ = Scaling(pd.DataFrame(out))
outZ = pd.concat([final_out['AA'], outZ], axis=1)
outZ = outZ.fillna(value=0)

final_out.to_csv(pdb_file+'_scoresEnergetics', sep='\t', index=False)
outZ.to_csv(pdb_file+'_scoresEnergeticsZ', sep='\t', index=False)
print('Write files', time.time()-t0)


