import pandas as pd
import numpy as np
import sys, re
import igraph
import time
t0 = time.time()

linearcutoff = 1
# removeMCMC = False
directed = False
# secondaryStructure2 = False
sidechainMode = True #build network of just SCSC, SCMC and MCSC. Still use unique atoms for the _uniqueAtoms directed plots
uniqueAtomsMode = False #build network of just unique atoms (each edge has to have at least 1). Still use unique atoms for the _uniqueAtoms directed plots
useDNA=True

# uniqueAtomsUnbiased=False
# uniqueAtomsGaurav=False
# uniqueAtomsOLD=False
uniqueAtomsGauravPP = False

weighted = True
removeWaters=True

#if Centroid folder exist
ligand_centroid_mode = True
ligandmode = False
# read file.pdb 
pdb = sys.argv[1] # input file is 4eiy.pdb
pdb_file = re.sub(".pdb","",pdb) # take name of file

# all processed must be in the same folder
# read ligand file, ex. 1BTL_ligand
ligand_file = open(pdb_file+'_ligand', 'r')
ligand_list = []
for line in ligand_file:
    ligand_list.append(line)
if len(ligand_list) != 0:
    liganddat = pd.read_csv(pdb_file+'_ligand', sep='\t', header=None)
    ligandmode = True


terminalAtoms = {}
# In folder must be this file or can create dict here?
with open('uniqueAtoms with tab', 'r') as uniqueAtomFile:
    for line in uniqueAtomFile:
        line= line.strip('\n').split()
        terminalAtoms[line[0]] = [x.strip(',') for x in line[1:]]

# CentroidNetLigand file 
if len(open(pdb_file+'_centroidNetLigand', 'r').readlines()) != 0:
    ligandCentroidfile = pd.read_csv(pdb_file+'_centroidNetLigand', sep='\t', header=None)
else:
    ligandCentroidfile = []

# read secondary structure file and write phi angles
secStructure = pd.read_csv(pdb_file+'_secondaryStructure', sep='\t', header=None) 
phiAngleDB = {}
for acid, angle in secStructure[[0, 2]].values:
    phiAngleDB[acid] = angle

# for file with ligands, add ligands to phi file with 0 angle
# to secStructure with angle of 360
if ligandmode: 
    liganddat = pd.read_csv(pdb_file+'_ligand', sep='\t', header=None)
    ligands = liganddat[1].unique()
    row_list = []
    for ligand in ligands:  
        dict1 = {0:ligand.replace('-', ''), 1:'ligand', 2:360}
        row_list.append(dict1)
        phiAngleDB[ligand] = 0
    secStructure = secStructure.append(pd.DataFrame(row_list))
else:
    liganddat = []
    ligands = []

# read RSA file
rsa = pd.read_csv(pdb_file+'.rsa', header=None, sep='\t', )
# add ligands to rsa with rsa values of -1
row_list = []
for ligand in ligands:
    dict1 = {0:ligand.replace('-', ''), 1:-1}
    row_list.append(dict1)
rsa = rsa.append(pd.DataFrame(row_list))

print('Read files and define dicts', time.time()-t0)

def addTerminalDetails(data): #largedata
    '''
    @descirption
    For observing reverse interactions, create data where sidechains, acids, atoms are replaced, and concat data
    Add 2 columns to data, whether acids' atom in terminalAtoms (1) or not (0)
    !!! largedata input isnot used so remove it 
    '''
    # concat to data taking into account the opposite connection
    replaced_data = data.replace({'MCSC':'SCMC', 'SCMC':'MCSC'})
    replaced_data = replaced_data[['A2', 'A1', 'Type1', 'Weight', 'Type2', 'a2', 'a1']]
    replaced_data = replaced_data.rename(columns={'A2':'A1', 'A1':'A2', 'a2':'a1', 'a1':'a2'})
    data = pd.concat([data, replaced_data], ignore_index=True)
    data = data.drop_duplicates() # delete same bonds

    # check if atoms are terminal atoms and add column of influence
    codes1 = []
    codes2 = []
    for row in data.values:
        acid1 = row[0][:3]
        acid2 = row[1][:3]
        code1 = 0
        code2 = 0
        if (row[5] in terminalAtoms[acid1]) or (row[4] == 'DNA' and useDNA==True):
            code1 = 1
        if (row[6] in terminalAtoms[acid2]) or (row[4] == 'DNA' and useDNA == True): 
            code2 = 1
        if (row[4] == 'PICATION' or row[4] == 'PIPI'): 
            code1 = 1
            code2 = 1
#         #uniqueAtomsGaurav is FALSE, so comment this part
#         if row[4] == 'PP':
#             if uniqueAtomsGaurav:
#                 if acid1 == 'GLY':
#                     code1=1
#                     row[2] = 'SCSC'
#                 if acid2 == 'GLY':
#                     code2 = 1
#                     row[2] = 'SCSC'
        codes1.append(code1)
        codes2.append(code2)
    data['influence1'] = codes1
    data['influence2'] = codes2
    return data

def collapse_agnostic1(data):
    '''
    Sort acids in data and create pairs
    For pairs of acids select subdata, sort acids in subdata
    replace sidechain and summarize weights
    '''
    values = data.values
    tmp = []
    for row in values:
        if (np.argsort(row[0:2])[0] == 1):
            row[0:2] = row[0:2][::-1]
        if row[2] == 'SCMC':
            row[2] = 'MCSC'
        if list(row[0:4]) not in tmp:
            tmp.append(list(row[0:4]))

    out = pd.DataFrame(tmp)
    df = out.groupby(by=[0,1], ).sum()
    out = []
    for aminos, value in zip(df.index, df.values):
        if weighted == False:
            out.append([aminos[0], aminos[1], 'mixed', 1])
        else:
            out.append([aminos[0], aminos[1], 'mixed', float(value)])
    return pd.DataFrame(out)


def collapse_directed1(data):
    '''
    Group data by acids and summarize weights
    '''
    df = data.groupby(by=['A1', 'A2']).sum()
    out = []
    for aminos, value in zip(df.index, df.values):
        if weighted == False:
            out.append([aminos[0], aminos[1], 'mixed', 1])
        else:
            out.append([aminos[0], aminos[1], 'mixed', float(value)])
    return pd.DataFrame(out)

# # if diff bw 2 aminoacid indexes > cutoff
def condition(dataframe):
    out = []
    for row in dataframe.values.tolist():
        cond = np.abs(int(row[0][3:-1]) - int(row[1][3:-1]))
        if cond > linearcutoff:
            out.append(row)
    return pd.DataFrame(out)

#----------------------------------------------------------------------------------------------------
# read 'pdb_net' file, with all protein bonds
# where AA1, AA2 - aminoacids, 
# Type1 - sidechains type (MCMC); Type2 - bond type
# Weight - energy of interaction, 
# Atom1, Atom2 - involved atoms
data_all = pd.read_csv(pdb_file+'_net', sep='\t', header=None, names=['A1', 'A2', 'Type1', 'Weight', 'Type2', 'a1', 'a2'])

# delete waters
if removeWaters:
    data_all = data_all.loc[-data_all['A1'].str.contains('HOH')]
    data_all = data_all.loc[-data_all['A2'].str.contains('HOH')]

nodes = (data_all['A1'].append(data_all['A2'])).unique();
sets = nodes
originalNodes = nodes #all nodes from pdb_net file
# here we can delete names of cols for comforatble work
# account all bonds with reverse connections
# data_all_reverse = data_all[['A2', 'A1', 'Type1', 'Weight', 'Type2', 'a2', 'a1']]
# data_all_reverse = data_all_reverse.rename(columns={'A2':'A1', 'A1':'A2', 'a2':'a1', 'a1':'a2'})
# data_all_original = pd.concat([data_all, data_all_reverse], axis=0, ignore_index=True)
data_all = data_all[data_all['Weight']>0] # take only data with positive weights

basedata = addTerminalDetails(data_all) 
basedata_noPP = collapse_agnostic1(basedata[basedata['Type2']!='PP'])

# Take non MainChain (MC) data
sidechaindata = basedata[basedata['Type1']!='MCMC'].iloc[:,0:4]
sidechaindata_detailed = basedata[basedata['Type1']!='MCMC']

influencedata_detailed = basedata
# take data where atoms are in terminal atoms
if uniqueAtomsMode:
    influencedata_detailed = basedata[(basedata['influence1']==1) | (basedata['influence2']==1)]

print('Define data_net', time.time()-t0)


influencedata = collapse_agnostic1(influencedata_detailed)
# influencedataGLY = influencedata[influencedata[0].str.contains('GLY')].append(influencedata[influencedata[0].str.contains('GLY')])
influencedata = condition(influencedata)


tmp = influencedata_detailed[influencedata_detailed['influence1']==1].iloc[:,0:4]
influencedata_directed = collapse_directed1(tmp)
# influencedata_directedGLY = influencedata_directed[influencedata_directed[0].str.contains('GLY')].append(influencedata_directed[influencedata_directed[1].str.contains('GLY')])
influencedata_directed = condition(influencedata_directed)

# CREATE NETWORK
influencenet = igraph.Graph.TupleList(edges=influencedata[[0,1]].values.tolist(), directed=False)
# WEIGHTS = 1/values
influencenet_weight = [1/item[3] for item in influencedata.values]

# CREATE DIRECTED NETWORK
influencenet_directed = igraph.Graph.TupleList(edges=influencedata_directed[[1,0]].values.tolist(), directed=True)
influencenet_directed_weight = [1/item[3] for item in influencedata_directed.values]


if sidechainMode:
    influencedata_detailed = sidechaindata_detailed
    influencedata = collapse_agnostic1(sidechaindata_detailed)
#     influencedataGLY = influencedata[influencedata[0].str.contains('GLY')].append(influencedata[influencedata[1].str.contains('GLY')])
    influencedata = condition(influencedata)
    
#     if uniqueAtomsGauravPP:
#         influencedata = influencedata.append(influencedataGLY)
    influencedata = influencedata[influencedata[3] > 0]

    influencenet = igraph.Graph.TupleList(edges=influencedata[[0,1]].values.tolist(), directed=False)
    influencenet_weight = [1/item[3] for item in influencedata.values]

    tmp = influencedata_detailed[influencedata_detailed['influence1']==1].iloc[:,0:4]
    influencedata_directed = collapse_directed1(tmp)
#     influencedata_directedGLY = influencedata_directed[influencedata_directed[0].str.contains('GLY')].append(influencedata_directed[influencedata_directed[1].str.contains('GLY')])
    influencedata_directed = condition(influencedata_directed)
#     if uniqueAtomsGauravPP:
#         influencedata_directed = influencedata_directed.append(influencedata_directedGLY)
    
    influencenet_directed = igraph.Graph.TupleList(edges=influencedata_directed[[1,0]].values.tolist(), directed=True)
    influencenet_directed_weight = [1/item[3] for item in influencedata_directed.values]

nodes1 = set(influencedata[0].values.tolist()+ influencedata[1].values.tolist())
nodesMissing = [node for node in nodes if node not in nodes1]


bfactor = pd.read_csv(pdb_file+'_Bfactor', sep='\t', header=None)
#Replace missing nodes in secStrcture matrix
for node in nodes:
    if node not in secStructure[0].values:
        str_dict = {0:node, 1:'xxx', 2:360}
        secStructure = secStructure.append(str_dict, ignore_index=True)
#Replace missing nodes in rsa matrix
#Replace missing nodes in bfactor matrix
for node in nodes:
    if node not in rsa[0].values:
        dict1 = {0:node, 1:0}
        rsa = rsa.append(dict1, ignore_index=True)
    if node not in bfactor[0].values:
        dict1 = {0:node, 1:1}
        bfactor = bfactor.append(dict1, ignore_index=True)


#GIVE DNA attrs
# if any(data_all['Type1'].str.contains('DNA')) == True:
for row in data_all.values:
    if row[4] == 'DNA':
        phiAngleDB[row[0]] = 90
        if row[0] not in rsa[0].values:
            dict1 = {0:row[0], 1:0}
            rsa = rsa.append(dict1, ignore_index=True)
        if row[0] not in secStructure[0].values:
            str_dict = {0:row[0], 1:'DNA', 2:360}
            secStructure = secStructure.append(str_dict, ignore_index=True)
        if row[0] not in bfactor[0].values:
            str_dict = {0:row[0], 1:1}
            bfactor = bfactor.append(str_dict, ignore_index=True)

print('Start calc features for 1st net_community_vec', time.time()-t0)


# Begin centrality calculations
# Modular calculations
# Walktrap algorithm for clusters detection

# betweenness of nodes in network
node_betweenness = influencenet.betweenness(weights=influencenet_weight, directed=directed)
node_betweenness_unique = influencenet_directed.betweenness(weights=influencenet_directed_weight, directed=directed)

net = igraph.Graph.TupleList(edges=basedata_noPP[[0,1]].values.tolist(), directed=False)
net_weight = [1/item[3] for item in basedata_noPP.values]

# Use Walktrap algorithm to select clusters   
net_community = net.community_walktrap(weights=net_weight, steps=4)
cluster = net_community.as_clustering()
memberships = cluster.membership
# Here we create dict of vertex names and their cluster belonging
net_community_vec = {}
for name, member in zip(net.vs['name'], memberships):
    net_community_vec[name] = member

net_community_vec_wt = net_community_vec
nodes = influencenet.vs['name']


def calculate_graph_attrs(net_community_vec):
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
    # if two vertex from the same cluster add betwenness value
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
        if len(subdata) == 0:
            node_intermodular_degree[node] = 0
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
    node_edge_betweenness_sidechain = {}
    for node in nodes:
        node_edge_betweenness_sidechain[node] = 0 
    for i in range(len(edge_betweenness)):
        if net_community_vec[influencedata[0][i]] != net_community_vec[influencedata[1][i]]:
            data = influencedata_directed[(influencedata_directed[0] == influencedata[0][i]) & (influencedata_directed[1] == influencedata[1][i])]
            if len(data) == 0: weight1 = 0 
            else: weight1 = data.values[0][3]
            data = influencedata_directed[(influencedata_directed[0] == influencedata[1][i]) & (influencedata_directed[1] == influencedata[0][i])]
            if len(data) == 0: weight2 = 0 
            else: weight2 = data.values[0][3]
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
            length = 0
            length1 = 0
            for key in node_modules[node]:
                length += len(node_modules[key])
            secondOrder_node_intermodular_degree[node] = length
            
            for key in node_modules_sidechain[node]:
                length1 += len(node_modules_sidechain[key])
            secondOrder_node_intermodular_degree_sidechain[node] = length1

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

    if node in influencenet_directed.vs['name']:
        firstorder_sidechain = influencenet_directed.neighbors(node, mode='in')
        num1 = influencenet_directed.vs['name'].index(node)
    else:
        firstorder_sidechain = []    
    secondorder = []
    for neigh in firstorder_sidechain:
        secondorder.append(influencenet_directed.neighbors(neigh, mode='in'))
    secondorder = list(set([i for sub in secondorder for i in sub if i!=num1]))
    secondOrderDegree_sidechain[node] = len(secondorder)


result = calculate_graph_attrs(net_community_vec)
node_edge_betweenness_wt = result[0]
node_edge_betweenness_sidechain_wt = result[1]
node_intermodular_degree_wt = result[2]
node_intermodular_degree_sidechain_wt = result[3]
secondOrder_node_intermodular_degree_wt = result[4]
secondOrder_node_intermodular_degree_sidechain_wt = result[5]
print('Calculate for 1st net_community_vec', time.time()-t0)
# [sum(node_edge_betweenness_wt.values()), sum(node_edge_betweenness_sidechain_wt.values()), sum(node_intermodular_degree_wt.values()),
#  sum(node_intermodular_degree_sidechain_wt.values()), sum(secondOrder_node_intermodular_degree_wt.values()),
#  sum(secondOrder_node_intermodular_degree_sidechain_wt.values())]

print('Start calc feature for 2nd net_community_vec', time.time()-t0)
# Here we define new clusters from secondaryStructure file
# and for new netwrok cluster calculate centrality attributes
secStructure = secStructure[[0, 1]]
factors = sorted(secStructure[1].unique(), key=lambda v: v.upper())
struct_dict = {}
for i in range(1, len(factors)+1):
    struct_dict[factors[i-1]] = i
net_community_vec = {}
for row in secStructure.values:
    net_community_vec[row[0]] = struct_dict[row[1]]

nodes = influencenet.vs['name']

result = calculate_graph_attrs(net_community_vec)

node_edge_betweenness_stride = result[0]
node_edge_betweenness_sidechain_stride = result[1]
node_intermodular_degree_stride = result[2]
node_intermodular_degree_sidechain_stride = result[3]
secondOrder_node_intermodular_degree_stride = result[4]
secondOrder_node_intermodular_degree_sidechain_stride = result[5]

print('End calc feature for 2nd net_community_vec', time.time()-t0)
# [sum(node_edge_betweenness_stride.values()), sum(node_edge_betweenness_sidechain_stride.values()), 
#  sum(node_intermodular_degree_stride.values()),
#  sum(node_intermodular_degree_sidechain_stride.values()), sum(secondOrder_node_intermodular_degree_stride.values()),
#  sum(secondOrder_node_intermodular_degree_sidechain_stride.values())]
print('Start additional featuring', time.time()-t0)
nodes = [node for node in nodes if node!='DNA1000A']
nodes = [node for node in nodes if node not in ['DA', 'DG', 'DC', 'DT']]
if len(liganddat)!=0:
    ligands = liganddat[1].unique()
    ligands = [l.replace('-', '') for l in ligands]
    nodes = [node for node in nodes if node not in ligands]
    nodes = [node for node in nodes if node in bfactor[0].values]

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

### FOR RSA DATA
df = pd.DataFrame()
out = []
for row in rsa.values:
    if row[0] in nodes and row[0] not in out:
        out.append(row[0])
        dict1 = {0: row[0], 1: row[1]}
        df = df.append(dict1, ignore_index=True)


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

rdsa = df
out = [rdsa,
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
colnames = ["AA","RSA","Degree","Degree_uniqueAtoms","SecondOrderDegree","SecondOrderDegree_uniqueAtoms",
 "NodeEdgeBetweennessSTRIDE","NodeEdgeBetweennessSTRIDE_unqiueAtoms","IntermodularDegreeSTRIDE",
 "IntermodularDegreeSTRIDE_uniqueAtoms","SecondOrderIntermodularDegreeSTRIDE",
 "SecondOrderIntermodularDegreeSTRIDE_uniqueAtoms","NodeEdgeBetweennessWALKTRAP",
 "NodeEdgeBetweennessWALKTRAP_uniqueAtoms","IntermodularDegreeWALKTRAP",
 "IntermodularDegreeWALKTRAP_uniqueAtoms","SecondOrderIntermodularDegreeWALKTRAP",
 "SecondOrderIntermodularDegreeWALKTRAP_uniqueAtoms"]

df_final.columns = colnames

# Account for ligands bonds
if ligandmode == True and ligand_centroid_mode == False:
    liganddat = pd.read_csv(pdb_file+'_ligand', sep='\t', header=None, names=['AA', 'ligand', 'value'])
    ligands = liganddat.ligand.unique()
    tmp = []
    for row in liganddat.values:
        acid = row[0]
        if acid.split('-')[1] in terminalAtoms[acid[:3]]:
            tmp.append([acid.split('-'), row[1], row[2]])
        else:
            continue

if ligandmode == True and ligand_centroid_mode == True:
    liganddat = pd.read_csv(pdb_file+'_centroidNetLigand', sep='\t', header=None, names=['AA', 'ligand', 'value'])
    tmp = []
    for row in liganddat.values:
        tmp.append([row[0], row[1], row[2]])
    out = pd.DataFrame(tmp, columns=['AA', 'ligand', 'value'])
    ligandvec = out.groupby('AA')['value'].min().to_dict()

if ligandmode == False:
    ligandvec = {}
    for node in nodes:
        ligandvec[node] = 0

ligand_df = Nanfill_df(ligandvec)
ligandvec = Nanfill_df(ligandvec)

ligand_df = ligand_df.loc[ligand_df[0].isin(df_final['AA'].tolist())]
ligand_df.columns =['AA', 'Ligand']

out = pd.merge(df_final, ligand_df, on='AA')

# OriginalNodes - nodes in pdb_net file; nodes - from influencenet graph
# if original nodes not in nodes, add to list of zeronodes
# and create dataframe from these zeronodes with same columns
# make final output (out11)
zeroNodes = []
for i in originalNodes:
    if i not in nodes:
        zeroNodes.append(i) 

zeroMat = pd.DataFrame(0, index=np.arange(len(zeroNodes)), columns=out.columns[1:])
zeroMat['AA'] = zeroNodes

out11 = pd.concat([out, zeroMat])
out = out11

for node in zeroNodes:
    # rsa[rsa[0] == node][1]
    out.loc[out['AA']==node, 'RSA'] = rsa[rsa[0] == node][1].values[0]  
    # NEED OR NOT
    out.loc[out['AA']==node, 'Ligand'] = ligandvec[ligandvec[0] == node][1].values[0] 

# fill None values with zero and add 1 to all rows
out = out.fillna(value=0)
out.loc[out['Ligand'] == 0, 'Ligand'] = 150
out3 = out[out.columns[2:]]+1
final_out = pd.concat([out[['AA', 'RSA']], out3], axis=1)
# For checking sum weight with R code
for col in final_out.columns:
    print(col, final_out[col].sum())

final_out1 = final_out.copy()

# normalization of data
def Scaling(df):
    cols = df.columns[1:]
    for col in cols:
        std = np.std(df[col], ddof=1)
        mean = np.mean(df[col])
        df[col] = (df[col] - mean)/std
    return df

out = final_out[final_out.columns[1:]]
# outZ = pd.DataFrame(scale(out), index=out.index, columns=out.columns)
outZ = Scaling(final_out)
# outZ = pd.concat([final_out['AA'], outZ], axis=1)
outZ = outZ.fillna(value=0)

# final_out1.to_csv(pdb_file+'_scoresCentroid', sep='\t', index=False)
# outZ.to_csv(pdb_file+'_scoresCentroidZ', sep='\t', index=False)
print('Write files', time.time()-t0)
