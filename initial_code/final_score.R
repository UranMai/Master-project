## run from within the monomer folder
## there can only be one *_scoresEnergetics and one Centroid/*_scoresCentroid in this folder
## args[1] is the folder for the multimer
## NOTE: results may differ slightly from publication because of missing data in the mutagenesis experiments used in teh paper

#args <- commandArgs(trailingOnly = TRUE)

#Load monomer
data1 = read.table(dir("./","_scoresEnergetics$"),header=T,row.names=1)
#print(data1[,1:2])
data1Z = as.data.frame(apply(data1,2,function(x){return(scale(as.numeric(x)))}))
rownames(data1Z) = rownames(data1)
data1Z[is.na(data1Z)] = 0
data2 = read.table(dir("./","_scoresCentroid$", full.names=TRUE),header=T,row.names=1)
data2Z = as.data.frame(apply(data2,2,function(x){return(scale(as.numeric(x)))}))
rownames(data2Z) = rownames(data2)
data2Z[is.na(data2Z)] = 0
colnames(data2Z) = paste(colnames(data2Z),"CENTROIDSC",sep="")
data2Z = data2Z[rownames(data1Z),]
data2 = data2[rownames(data1),]

#Load multimer
data1Multimer = read.table(dir("./","_scoresEnergetics$"),header=T,row.names=1)
data2Multimer = read.table(dir("./","_scoresCentroid$",full.names=TRUE),header=T,row.names=1)
data2Multimer = data2Multimer[rownames(data1Multimer),]
data1MultimerZ = as.data.frame(apply(data1Multimer,2,function(x){return(scale(as.numeric(x)))}))
rownames(data1MultimerZ) = rownames(data1Multimer)
# data1MultimerZ = round_df(data1MultimerZ, 6)
data2MultimerZ = as.data.frame(apply(data2Multimer,2,function(x){return(scale(as.numeric(x)))}))
rownames(data2MultimerZ) = rownames(data2Multimer)
colnames(data1MultimerZ) = paste(colnames(data1MultimerZ),"MULTIMER",sep="")
colnames(data2MultimerZ) = paste(colnames(data2MultimerZ),"MULTIMERCENTROIDSC",sep="")
data2MultimerZ = data2MultimerZ[rownames(data1MultimerZ),]

#Add RSA
data1 = cbind(data2[,"RSA"],data1)
colnames(data1)[1]="RSA"
dataZ = cbind(data1Z,data2Z)

#Identify max NodeEdgeBetweennessSTRIDE_sidechain
acidsMonomer = c()
for (i in 1:length(rownames(dataZ))) {
  node = rownames(dataZ)[i]
  node = unlist(strsplit(node,split=""))
  node = paste(node[1:(length(node)-1)],collapse="")
  acidsMonomer = c(acidsMonomer,node)
}

acidsMultimer = c()
for (i in 1:length(rownames(data1MultimerZ))) {
  node = rownames(data1MultimerZ)[i]
  node = unlist(strsplit(node,split=""))
  node = paste(node[1:(length(node)-1)],collapse="")
  acidsMultimer = c(acidsMultimer,node)
}

allAcids = unique(c(acidsMonomer,acidsMultimer))

#SecondOrderIntermodularDegree is the average of the multimer only
SecondOrderIntermodularDegree_ALL = (data1MultimerZ[,"SecondOrderIntermodularDegreeSTRIDE_sidechainMULTIMER"]+data1MultimerZ[,"SecondOrderIntermodularDegreeSTRIDEMULTIMER"]+data2MultimerZ[,"SecondOrderIntermodularDegreeSTRIDEMULTIMERCENTROIDSC"]+data2MultimerZ[,"SecondOrderIntermodularDegreeWALKTRAPMULTIMERCENTROIDSC"])/4
SecondOrderIntermodularDegree_AVERAGE = tapply(SecondOrderIntermodularDegree_ALL,as.factor(acidsMultimer),mean)
SecondOrderIntermodularDegree_AVERAGE = SecondOrderIntermodularDegree_AVERAGE[allAcids]
missing = which(allAcids%in%acidsMultimer == FALSE)
SecondOrderIntermodularDegree_AVERAGE[missing] = (data1Z[allAcids[which(allAcids%in%acidsMultimer == FALSE)],"SecondOrderIntermodularDegreeSTRIDE_sidechain"]+data1Z[allAcids[which(allAcids%in%acidsMultimer == FALSE)],"SecondOrderIntermodularDegreeSTRIDE"]+data2Z[allAcids[which(allAcids%in%acidsMultimer == FALSE)],"SecondOrderIntermodularDegreeSTRIDECENTROIDSC"]+data2Z[allAcids[which(allAcids%in%acidsMultimer == FALSE)],"SecondOrderIntermodularDegreeWALKTRAPCENTROIDSC"])/4
names(SecondOrderIntermodularDegree_AVERAGE)[missing] = allAcids[which(allAcids%in%acidsMultimer == FALSE)]

#NodeEdgeBetweenness is the max of the monomer/multimer
NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER = (dataZ[,"NodeEdgeBetweennessSTRIDE"]+dataZ[,"NodeEdgeBetweennessSTRIDE_sidechain"])/2
NodeEdgeBetweennessSTRIDE_AVERAGE_MULTIMER = (data1MultimerZ[,"NodeEdgeBetweennessSTRIDEMULTIMER"]+data1MultimerZ[,"NodeEdgeBetweennessSTRIDE_sidechainMULTIMER"])/2
NodeEdgeBetweennessSTRIDE_sidechain_ALL = c(NodeEdgeBetweennessSTRIDE_AVERAGE_MONOMER,NodeEdgeBetweennessSTRIDE_AVERAGE_MULTIMER)
NodeEdgeBetweennessSTRIDE_sidechain_MAX = tapply(NodeEdgeBetweennessSTRIDE_sidechain_ALL,as.factor(c(acidsMonomer,acidsMultimer)),max,na.rm=TRUE)
NodeEdgeBetweennessSTRIDE_sidechain_MAX = NodeEdgeBetweennessSTRIDE_sidechain_MAX[allAcids]

#Ligand is the min of the multimer only
LigandMULTIMERCENTROIDSC_ALL = data2MultimerZ[,"LigandMULTIMERCENTROIDSC"]
LigandMULTIMERCENTROIDSC_MIN = tapply(LigandMULTIMERCENTROIDSC_ALL,as.factor(acidsMultimer),min)
LigandMULTIMERCENTROIDSC_MIN = LigandMULTIMERCENTROIDSC_MIN[allAcids]

missing = which(allAcids%in%acidsMultimer == FALSE)
LigandMULTIMERCENTROIDSC_MIN[missing] = dataZ[allAcids[which(allAcids%in%acidsMultimer == FALSE)],"LigandCENTROIDSC"]
names(LigandMULTIMERCENTROIDSC_MIN)[missing] = allAcids[which(allAcids%in%acidsMultimer == FALSE)]
LigandMULTIMERCENTROIDSC_MIN[is.na(LigandMULTIMERCENTROIDSC_MIN)] = 0

#RSA
RSA_ALL = c(data1[,"RSA"],data2Multimer[,"RSA"])
RSA_MIN = tapply(RSA_ALL,as.factor(c(acidsMonomer,acidsMultimer)),min,na.rm=TRUE)
RSA_MIN = RSA_MIN[allAcids]

y = SecondOrderIntermodularDegree_AVERAGE+NodeEdgeBetweennessSTRIDE_sidechain_MAX-LigandMULTIMERCENTROIDSC_MIN
#print(sum(SecondOrderIntermodularDegree_AVERAGE))
#print(sum(NodeEdgeBetweennessSTRIDE_sidechain_MAX))
#print(sum(LigandMULTIMERCENTROIDSC_MIN))
print(sum(round(y,1)))
out=cbind(names(y),y)
write.table(out,file="FinalSum",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
