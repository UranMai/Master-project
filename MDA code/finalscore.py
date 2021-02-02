import pandas as pd
import sys, re
import time
# pdb_file = sys.argv[1]
t0 = time.time()
pdb = sys.argv[1]
pdb_file = re.sub(".pdb","",pdb)

# Read Energetics and Centroid files
data1 = pd.read_csv(pdb_file+'_scoresEnergetics', sep='\t', index_col='AA')
data1Z = pd.read_csv(pdb_file + '_scoresEnergeticsZ', sep='\t', index_col='AA')
data1 = data1.round(6)
data1Z = data1Z.round(6)
data1Z = data1Z.fillna(value=0)

data2 = pd.read_csv(pdb_file + '_scoresCentroid', sep='\t', index_col='AA')
data2Z = pd.read_csv(pdb_file + '_scoresCentroidZ', sep='\t', index_col='AA')
data2Z = data2Z.round(6)

data2Z.columns = data2Z.columns + 'CENTROIDSC'

data1Multimer = pd.read_csv(pdb_file + '_scoresEnergetics', sep='\t', index_col='AA')
data1MultimerZ = pd.read_csv(pdb_file + '_scoresEnergeticsZ', sep='\t', index_col='AA')
data1MultimerZ = data1MultimerZ.round(6)

data2Multimer = pd.read_csv(pdb_file + '_scoresCentroid', sep='\t', index_col='AA')
data2MultimerZ =  pd.read_csv(pdb_file + '_scoresCentroidZ', sep='\t', index_col='AA')
data2MultimerZ = data2MultimerZ.round(6)

# Rename columns
data1MultimerZ.columns = data1MultimerZ.columns + 'MULTIMER'
data2MultimerZ.columns = data2MultimerZ.columns + 'MULTIMERCENTROIDSC'

data1 = pd.concat([data2['RSA'], data1],axis=1)
dataZ = pd.concat([data1Z, data2Z], axis=1)

acidsMonomer = []
for i in range(len(dataZ)):
    node = dataZ.index[i]
    acidsMonomer.append(node[:-1])

acidsMultimer = []
for i in range(len(data1MultimerZ)):
    node = data1MultimerZ.index[i]
    acidsMultimer.append(node[:-1])

allAcids = list(set(acidsMonomer+acidsMultimer))

# Calculate Average value of SecondOrderIntermodularDegree
SecondOrderIntermodularDegree_ALL = (data1MultimerZ["SecondOrderIntermodularDegreeSTRIDE_sidechainMULTIMER"]+
                                     data1MultimerZ["SecondOrderIntermodularDegreeSTRIDEMULTIMER"]+
                                     data2MultimerZ['SecondOrderIntermodularDegreeSTRIDEMULTIMERCENTROIDSC']+
                                     data2MultimerZ["SecondOrderIntermodularDegreeWALKTRAPMULTIMERCENTROIDSC"])/4

arr = []
df = SecondOrderIntermodularDegree_ALL
for i, j in zip(df.index, df.values):
    acid = i[:-1] 
    arr.append([acid, j])

new_df = pd.DataFrame(arr)

SecondOrderIntermodularDegree_AVERAGE = {}

for acid in allAcids:
    value = new_df[new_df[0]==acid][1].mean()
    SecondOrderIntermodularDegree_AVERAGE[acid] = value

missing = [acid for acid in allAcids if acid not in acidsMultimer]
# for missing nodes calculate 
for node in missing:
    value = (data1Z[data1Z.index == node]["SecondOrderIntermodularDegreeSTRIDE_sidechain"]+
             data1Z[data1Z.index == node]["SecondOrderIntermodularDegreeSTRIDE"]+
             data2Z[data2Z.index == node]["SecondOrderIntermodularDegreeSTRIDECENTROIDSC"]+
             data2Z[data2Z.index == node]["SecondOrderIntermodularDegreeWALKTRAPCENTROIDSC"])/4
    SecondOrderIntermodularDegree_AVERAGE[node] = value

# Average value of NodeEdgeBetweenness values
NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER = (dataZ["NodeEdgeBetweennessSTRIDE"]+
                                             dataZ["NodeEdgeBetweennessSTRIDE_sidechain"])/2
NodeEdgeBetweennessSTRIDE_AVERAGE_MULTIMER = (data1MultimerZ["NodeEdgeBetweennessSTRIDEMULTIMER"]+
                                              data1MultimerZ["NodeEdgeBetweennessSTRIDE_sidechainMULTIMER"])/2
NodeEdgeBetweennessSTRIDE_sidechain_ALL = NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER.append(NodeEdgeBetweennessSTRIDE_AVERAGE_MULTIMER)

arr = []
df = NodeEdgeBetweennessSTRIDE_sidechain_ALL
for i, j in zip(df.index, df.values):
    acid = i[:-1] 
    arr.append([acid, j])
new_df = pd.DataFrame(arr)

NodeEdgeBetweennessSTRIDE_sidechain_MAX = {}

for acid in acidsMonomer+acidsMultimer:
    value = new_df[new_df[0]==acid][1].max()
    NodeEdgeBetweennessSTRIDE_sidechain_MAX[acid] = value

#  CALCULATIONS OF LIGAND VALUES (MIN VALUE)
LigandMULTIMERCENTROIDSC_ALL = data2MultimerZ["LigandMULTIMERCENTROIDSC"]
LigandMULTIMERCENTROIDSC_ALL.sum()

arr = []
df = LigandMULTIMERCENTROIDSC_ALL
for i, j in zip(df.index, df.values):
    acid = i[:-1] 
    arr.append([acid, j])
new_df = pd.DataFrame(arr)

LigandMULTIMERCENTROIDSC_MIN = {}

for acid in acidsMultimer:
    value = new_df[new_df[0]==acid][1].min()
    LigandMULTIMERCENTROIDSC_MIN[acid] = value

missing = [acid for acid in allAcids if acid not in acidsMultimer]
missing
for node in missing:
    value = dataZ["LigandCENTROIDSC"]
    LigandMULTIMERCENTROIDSC_MIN[node] = value

RSA_ALL = pd.concat([data1['RSA'], data2Multimer['RSA']])
RSA_MIN = {}
# for node in allAcids:
#     value = RSA_ALL[RSA_ALL.index.str.contains(node)].min()
#     RSA_MIN[node] = value

arr = []
df = RSA_ALL
for i, j in zip(df.index, df.values):
    acid = i[:-1] 
    arr.append([acid, j])
new_df = pd.DataFrame(arr)

for acid in acidsMonomer+acidsMultimer:
    value = new_df[new_df[0]==acid][1].min()
    RSA_MIN[acid] = value

# Calculations of SCORE
out = {}
for node in allAcids:
    y = (SecondOrderIntermodularDegree_AVERAGE[node]+NodeEdgeBetweennessSTRIDE_sidechain_MAX[node]-LigandMULTIMERCENTROIDSC_MIN[node])
    out[node] = y


pd.DataFrame(out.items()).sort_values(by=0).to_csv('FinalSum', sep='\t', index=False, header=False)
print('Create FinalSum', time.time()-t0)

# Check results with R.sum(y) (use round())
sum = 0 
for i, j in out.items():
   sum += round(j, 1) # change 1to2,3
print(sum)
