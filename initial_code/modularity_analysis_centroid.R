library(igraph)
library(stats)
## args <- commandArgs(trailingOnly = TRUE)
args = '4eiy' # for running with source in rstudio
##read in data
## arg1: col1 and col2 are edges, col3 is weight of edge
## arg2: weighted, unweighted <-- no longer using
## arg3: forced module definition <-- no longer using

# Set desired parameters
removeMCMC = FALSE
linearcutoff = 1
directed = FALSE
secondaryStructure2 = FALSE
sidechainMode = TRUE ##build network of just SCSC, SCMC and MCSC. Still use unique atoms for the _uniqueAtoms directed plots
uniqueAtomsMode = FALSE ##build network of just unique atoms (each edge has to have at least 1). Still use unique atoms for the _uniqueAtoms directed plots
useDNA=TRUE
uniqueAtomsUnbiased=FALSE
uniqueAtomsGaurav=FALSE
uniqueAtomsOLD=FALSE
atomCorrection = FALSE
uniqueAtomsGauravPP = FALSE
weighted = TRUE
removeWaters=TRUE
# we create Centroid folder and copy all created files to that folder
# so ligandCentroidMode = True
#if (tail(unlist(strsplit(getwd(),"/")),n=1)=="Centroid") {
#  ligandCentroidMode = TRUE
#} else {
#  ligandCentroidMode = FALSE
#}
ligandCentroidMode = TRUE #this influence on out['LIGAND']

# Read in ligand file
ligandfile = scan(paste(args[1],"_ligand",sep=""),what="character")
if (length(ligandfile)==0){
  ligandmode = FALSE
} else {
  ligandmode = TRUE
}

ligandCentroidFile = scan(paste(args[1],"_centroidNetLigand",sep=""),what="character")

# keyword is the base filename
# method is "frag" or "module" or "both"
keyword = args[1]
# regular expression in R: $ means the end of string and here it means that we don't use secStructure2!
args=c(paste(args[1],"_net",sep=""),"weighted",dir("./","_secondaryStructure$"))


# keyword = '1BTL'
# args=c(paste('1BTL',"_net",sep=""),"weighted",dir("./","_secondaryStructure$"))


##Load data
#Unique atoms file
if(uniqueAtomsUnbiased) {
  terminalAtomsFile = as.matrix(read.table("/Users/vjpatel/Dropbox/Gaurav/uniqueAtomsUnbiased",sep="\t"))
  terminalAtoms = list()
  for (i in 1:nrow(terminalAtomsFile)) {
    terminalAtoms[[terminalAtomsFile[i,1]]] = unlist(strsplit(terminalAtomsFile[i,2],split=","))
  }
} else if (uniqueAtomsGaurav) {
  terminalAtomsFile = as.matrix(read.table("/Users/vjpatel/Dropbox/Gaurav/uniqueAtomsGaurav",sep="\t"))
  terminalAtoms = list()
  for (i in 1:nrow(terminalAtomsFile)) {
    terminalAtoms[[terminalAtomsFile[i,1]]] = unlist(strsplit(terminalAtomsFile[i,2],split=","))
  }
} else if (uniqueAtomsOLD) {
  terminalAtomsFile = as.matrix(read.table("../../terminalAtoms",sep="\t"))
  terminalAtoms = list()
  for (i in 1:nrow(terminalAtomsFile)) {
    terminalAtoms[[terminalAtomsFile[i,1]]] = unlist(strsplit(terminalAtomsFile[i,2],split=","))
  }
} else {
  terminalAtomsFile = as.matrix(read.table("uniqueAtoms with tab",sep="\t"))
  terminalAtoms = list()
  otherTA = list()
  for (i in 1:nrow(terminalAtomsFile)) {
    terminalAtoms[[terminalAtomsFile[i,1]]] = unlist(strsplit(terminalAtomsFile[i,2],split=", "))
    #otherTA[[terminalAtomsFile[i,1]]] = unlist(strsplit(terminalAtomsFile[i,4],split=","))
  }
}

#Secondary structure file
#Ligands are all one secondary structure
secStructure = as.matrix(read.table(args[3]))
phiAngleDB = {}
for (i in 1:nrow(secStructure)) {
  phiAngleDB[[secStructure[i,1]]] = as.numeric(secStructure[i,3])
}

if(ligandmode) {
  liganddat = as.matrix(read.table(dir("./","_ligand"),sep="\t"))
  liganddat = unique(liganddat[,2])
  for (ligand in liganddat) {
    secStructure = rbind(secStructure,c(sub("-","",ligand),"ligand",360))
    phiAngleDB[[sub("-","",ligand)]]=0
  }
} else {
  liganddat = c()
}


##Define functions

##Collapse matrix functions
##Will return redundant matrix
##Input needs to be 7 columns
##Second input is different if the main input is a subset of the larger data
addTerminalDetails = function(data,largerdata) {
  #List edges bi-directionally
  x=unlist(lapply(data[,3],function(x){if(x=="MCSC"){return("SCMC")}else if(x=="SCMC"){return("MCSC")} else {return(x)}}))
  data = rbind(data,cbind(data[,2],data[,1],x,data[,4],data[,5],data[,7],data[,6]))
  data = unique(data)
  out = c()
  out_noPP = c()
  out_details = c()
  for (i in 1:nrow(data)) {
    acid1 = paste(unlist(strsplit(data[i,1],split=""))[1:3],collapse="")
    acid2 = paste(unlist(strsplit(data[i,2],split=""))[1:3],collapse="")
    subdata = matrix(data[data[,1]==data[i,1]&data[,2]==data[i,2],],ncol=7)
    allsubdata = matrix(largerdata[largerdata[,1]==data[i,1],],ncol=7)
    code1=0;code2=0
    if (data[i,6] %in% terminalAtoms[[acid1]]) {
      code1 = 1
    }
    if (data[i,5] == "DNA" & useDNA==TRUE) {
      code1 = 1
    }
    allsubdata = matrix(largerdata[largerdata[,2]==data[i,2],],ncol=7)
    if (data[i,7] %in% terminalAtoms[[acid2]]) {
      code2 = 1
    }
    if (data[i,5] == "DNA" & useDNA==TRUE) {
      code2 = 1
    }
    if (data[i,5]=="PICATION"|data[i,5]=="PIPI") {
      code1=1
      code2=1
    }
    if (data[i,5]=="PP") {
      if (uniqueAtomsGaurav) {
        if((acid1) == "GLY") {
          code1 = 1
          data[i,3]=="SCSC"
        }
        if ((acid2) == "GLY") {
          code2 = 1
          data[i,3] == "SCSC"
        }
      }
    }
    out_details = rbind(out_details, c(data[i,],code1,code2))
  }
  return(out_details)
}

##Will return non-redundant matrix
collapse_agnostic = function(data) {
  if(nrow(data)==0) {
    return(matrix(ncol=4,nrow=0))
  } else if (ncol(data)==4 & data[1,3]=="mixed") {
    return(data)    
  } else {
    pairs = unique(t(apply(matrix(data[,1:2],ncol=2),1,function(x){return(sort(x))})))
    out = c()
    for (i in 1:nrow(pairs)) {
      subdat = matrix(data[(data[,1]==pairs[i,1]&data[,2]==pairs[i,2]) | (data[,2]==pairs[i,1]&data[,1]==pairs[i,2]),1:4],ncol=4)
      subdat = unique(t(apply(subdat,1,function(x){a=which(x==sort(x[1:2])[1]);
      if(a==1){return(x)}else{return(c(x[2],x[1],paste(unlist(strsplit(x[3],""))[c(3,4,1,2)],collapse=""),x[4]))}})))
      if (weighted==FALSE) {
        out = rbind(out,c(subdat[1,1],subdat[1,2],"mixed",1))
      } else {
        out = rbind(out,c(subdat[1,1],subdat[1,2],"mixed",sum(as.numeric(subdat[,4]))))
      }
    }
  }
  return(out)
}

##Will return directed collapsed matrix
collapse_directed = function(data) {
  if (length(data)==0) {
    return(matrix(nrow=0,ncol=4))
  } else {
    pairs = matrix(unique(data[,1:2]),ncol=2)
    out = c()
    for (i in 1:nrow(pairs)) {
      subdat = matrix(data[data[,1]==pairs[i,1]&data[,2]==pairs[i,2],1:4],ncol=4)
      if (weighted==FALSE) {
        out = rbind(out,c(subdat[1,1],subdat[1,2],"mixed",1))
      } else {
        out = rbind(out,c(subdat[1,1],subdat[1,2],"mixed",sum(as.numeric(subdat[,4]))))
      }
    }
    return(matrix(out,ncol=4))
  }
}

rsa = as.matrix(read.table(dir("./",".rsa$")))
rownames(rsa) = rsa[,1]
rsa[,1] = rsa[,2]

#add ligand as -1
for (ligand in liganddat) {
  rsa = rbind(rsa,c(-1,-1))
  rownames(rsa)[nrow(rsa)] = sub("-","",ligand)
}

#Function to get amino acid
getAcid = function(residue) {
  if(all((unlist(strsplit(residue,""))[1:2] == c("D","A")) == TRUE) | all((unlist(strsplit(residue,""))[1:2] == c("D","G"))==TRUE) | all((unlist(strsplit(residue,""))[1:2] == c("D","C"))==TRUE) | all((unlist(strsplit(residue,""))[1:2] == c("D","T"))==TRUE)) {return(paste(unlist(strsplit(residue,""))[1:2],collapse=""))
  } else {
    return(paste(unlist(strsplit(residue,""))[1:3],collapse=""))
  }
}

##################
## read in data ##
##################
data_all = as.matrix(read.table(args[1]))[,1:7]

if(removeWaters) {
  y=c(grep("HOH",data_all[,1]),grep("HOH",data_all[,2]))
  if(length(y)>0) {
    data_all = data_all[-y,]
  }
}

nodes = unique(c(data_all[,1:2])); sets = nodes
originalNodes = nodes
data_all_original = rbind(data_all,data_all[,c(2,1,3,4,5,7,6)])
data_all = data_all[as.numeric(data_all[,4])>0,]    


basedata = addTerminalDetails(data_all,data_all_original)
basedata_noPP = collapse_agnostic(basedata[basedata[,5]!="PP",])

sidechaindata = basedata[basedata[,3]!="MCMC",1:4]
sidechaindata_detailed = basedata[basedata[,3]!="MCMC",]

influencedata_detailed = basedata
if (uniqueAtomsMode) {
  influencedata_detailed = basedata[basedata[,8]==1 | basedata[,9]==1,]
}


influencedata = collapse_agnostic(influencedata_detailed)
influencedataGLY = rbind(influencedata[grep("GLY",influencedata[,1]),],influencedata[grep("GLY",influencedata[,1]),])
influencedata = influencedata[abs(as.numeric(gsub("[A-Z]","",influencedata[,1]))-as.numeric(gsub("[A-Z]","",influencedata[,2])))>linearcutoff,]

influencedata_directed = collapse_directed(matrix(influencedata_detailed[influencedata_detailed[,8]==1,1:4],ncol=4))
influencedata_directedGLY = rbind(influencedata_directed[grep("GLY",influencedata_directed[,1]),],influencedata_directed[grep("GLY",influencedata_directed[,2]),])
influencedata_directed = influencedata_directed[abs(as.numeric(gsub("[A-Z]","",influencedata_directed[,1]))-as.numeric(gsub("[A-Z]","",influencedata_directed[,2])))>linearcutoff,]

influencenet = graph.edgelist(influencedata[,1:2],directed=FALSE)
influencenet$weight = 1/as.numeric(influencedata[,4])
influencedata_directed = matrix(influencedata_directed,ncol=4)
influencenet_directed = graph.edgelist(matrix(influencedata_directed[,c(2,1)],ncol=2),directed=TRUE)
influencenet_directed$weight = 1/as.numeric(influencedata_directed[,4])


if(sidechainMode) {
  influencedata_detailed = sidechaindata_detailed
  influencedata = collapse_agnostic(influencedata_detailed)
  influencedataGLY = rbind(influencedata[grep("GLY",influencedata[,1]),],influencedata[grep("GLY",influencedata[,2]),])
  influencedata = influencedata[abs(as.numeric(gsub("[A-Z]","",influencedata[,1]))-as.numeric(gsub("[A-Z]","",influencedata[,2])))>linearcutoff,]
  #Add back for GLY if uniqueAtomsGaurav
  if(uniqueAtomsGauravPP) {
    influencedata = rbind(influencedata,influencedataGLY)
  }
  influencedata = influencedata[as.numeric(influencedata[,4])>0,]
  
  influencenet = graph.edgelist(influencedata[,1:2],directed=FALSE)
  influencenet$weight = 1/as.numeric(influencedata[,4])
  
  influencedata_directed = collapse_directed(matrix(influencedata_detailed[influencedata_detailed[,8]==1,1:4],ncol=4))
  influencedata_directedGLY = rbind(influencedata_directed[grep("GLY",influencedata_directed[,1]),],influencedata_directed[grep("GLY",influencedata_directed[,2]),])
  
  influencedata_directed = influencedata_directed[abs(as.numeric(gsub("[A-Z]","",influencedata_directed[,1]))-as.numeric(gsub("[A-Z]","",influencedata_directed[,2])))>linearcutoff,]
  #Add back for GLY if uniqueAtomsGaurav
  if(uniqueAtomsGauravPP) {
    influencedata_directed = rbind(influencedata_directed,influencedata_directedGLY)
  }
  influencedata_directed = matrix(influencedata_directed,ncol=4)
  influencenet_directed = graph.edgelist(matrix(influencedata_directed[,c(2,1)],ncol=2),directed=TRUE)
  influencenet_directed$weight = 1/as.numeric(influencedata_directed[,4])
  
}

#Replace missing nodes in the igraph net
nodesMissing = nodes[nodes %in% influencedata[,1]==FALSE & nodes %in% influencedata[,2]==FALSE]

##B-factor
bfactor = as.matrix(read.table(dir("./","_Bfactor")))
rownames(bfactor) = bfactor[,1]
#Replace missing nodes in secStructure matrix
if (any(nodes %in% secStructure[,1]==FALSE)) {
  for (node in nodes) {
    if (node %in% secStructure[,1]==FALSE) {
      secStructure = rbind(secStructure,c(node,"xxx",360))
    }
  }
}
for (node in nodes) {
  if (node %in% rownames(rsa) == FALSE) {
    rsa = rbind(rsa,c(0,0))
    rownames(rsa)[nrow(rsa)] = node
  }
  if (node %in% rownames(bfactor) == FALSE) {
    bfactor = rbind(bfactor,1)
    rownames(bfactor)[nrow(bfactor)] = node
  }
}




#Give DNA attributes
if ("DNA" %in% data_all[,5]) {
  for (i in 1:nrow(data_all)) {
    if (data_all[i,5]=="DNA") {
      phiAngleDB[[data_all[i,1]]] = 90
      if (data_all[i,1] %in% rownames(rsa) == FALSE) {
        rsa = rbind(rsa,c(0,0))
        rownames(rsa)[nrow(rsa)] = data_all[i,1]
      }
      if (data_all[i,1] %in% rownames(secStructure) == FALSE) {
        secStructure = rbind(secStructure,c(data_all[i,1],"DNA",360))
      }
      if (data_all[i,1] %in% rownames(bfactor) == FALSE) {
        bfactor = rbind(bfactor,1)
        rownames(bfactor)[nrow(bfactor)] = data_all[i,1]
      }
    }
  }
}


#Betweenness
node_betweenness = betweenness(influencenet,weights=influencenet$weight,directed=directed)
node_betweenness_unique = betweenness(influencenet_directed,weights=influencenet_directed$weight,directed=directed)

#Modular calculations

##WALKTRAP
net = graph.edgelist(basedata_noPP[,1:2],directed=FALSE)
net$weight = 1/as.numeric(basedata_noPP[,4])
net_community=walktrap.community(net,weights=net$weight); 
net_community_vec = net_community$membership; 
names(net_community_vec)=net_community$names

net_community_vec_wt = net_community_vec
#net_community_vec[buriednodes] = max(net_community_vec)+1
nodes = V(influencenet)$name

# edgecolors = rep("grey90",nrow(basedata_noPP))
# nodes.color = c()
# colorPalette=rainbow(max(net_community_vec),s=.5)
# for (i in 1:length(nodes)) {
#     nodes.color = c(nodes.color,colorPalette[net_community_vec[nodes[i]]])
# }
# for (i in 1:nrow(influencedata)) {
#     j = which((basedata_noPP[,1]==influencedata[i,1]&basedata_noPP[,2]==influencedata[i,2])|(basedata_noPP[,1]==influencedata[i,2]&basedata_noPP[,2]==influencedata[i,1]))
#     if (net_community_vec[influencedata[i,1]]!=net_community_vec[influencedata[i,2]]) { edgecolors[j]="grey40" }
# }

##Weighted edge betweenness WALKTRAP
edge_betweenness = edge.betweenness(influencenet,weights=influencenet$weight)
node_edge_betweenness = rep(0,length(nodes))
names(node_edge_betweenness) = nodes

for(i in 1:length(edge_betweenness)) {
  if (net_community_vec[influencedata[i,1]]!=net_community_vec[influencedata[i,2]]) {
    node_edge_betweenness[influencedata[i,1]] = node_edge_betweenness[influencedata[i,1]]+edge_betweenness[i]
    node_edge_betweenness[influencedata[i,2]] = node_edge_betweenness[influencedata[i,2]]+edge_betweenness[i]
  }
}

node_intermodular_degree = rep(0,length(nodes))
names(node_intermodular_degree) = nodes
node_modules = list()
for (node in nodes) {
  subdata = matrix(influencedata[influencedata[,1]==node | influencedata[,2]==node,],ncol=4)
  if (nrow(subdata)==0) {node_intermodular_degree[node] = 0; next}
  bound = unique(c(subdata[,1:2]))
  bound = bound[bound!=node]
  bound_modules = net_community_vec[bound]
  bound_modules = bound_modules[bound_modules!=net_community_vec[node]]
  node_modules[[node]] = bound_modules
  node_intermodular_degree[node] = length(bound_modules)
}

##Weighted edge betweenness WALKTRAP - WEIGHT BY SIDE CHAIN
node_edge_betweenness_sidechain = rep(0,length(nodes))
names(node_edge_betweenness_sidechain) = nodes
for(i in 1:length(edge_betweenness)) {
  if (net_community_vec[influencedata[i,1]]!=net_community_vec[influencedata[i,2]]) {
    weight1 = as.numeric(influencedata_directed[influencedata_directed[,1]==influencedata[i,1]&influencedata_directed[,2]==influencedata[i,2],4]);if(length(weight1)==0){weight1=0}
    weight2 = as.numeric(influencedata_directed[influencedata_directed[,1]==influencedata[i,2]&influencedata_directed[,2]==influencedata[i,1],4]);if(length(weight2)==0){weight2=0}
    if(weight1==0 & weight2==0) {
      node_edge_betweenness_sidechain[influencedata[i,1]] = node_edge_betweenness_sidechain[influencedata[i,1]] + 0
      node_edge_betweenness_sidechain[influencedata[i,2]] = node_edge_betweenness_sidechain[influencedata[i,2]] + 0
    } else {
      node_edge_betweenness_sidechain[influencedata[i,1]] = node_edge_betweenness_sidechain[influencedata[i,1]]+(weight1/(weight1+weight2))*edge_betweenness[i]
      node_edge_betweenness_sidechain[influencedata[i,2]] = node_edge_betweenness_sidechain[influencedata[i,2]]+(weight2/(weight1+weight2))*edge_betweenness[i]
    }
  }
}
node_intermodular_degree_sidechain = rep(0,length(nodes))
names(node_intermodular_degree_sidechain) = nodes
node_modules_sidechain = list()
for (node in nodes) {
  subdata_sidechain = matrix(influencedata_directed[influencedata_directed[,1]==node,],ncol=4)
  if (nrow(subdata_sidechain)==0) {node_intermodular_degree_sidechain[node] = 0; next}
  bound = unique(c(subdata_sidechain[,1:2])); bound = bound[bound!=node]
  bound_modules = net_community_vec[bound]; bound_modules = bound_modules[bound_modules!=net_community_vec[node]]
  node_modules_sidechain[[node]] = bound_modules
  node_intermodular_degree_sidechain[node] = length(bound_modules)
}

##Degree and second order degree for WALKTRAP
firstOrderDegree = degree(influencenet)
firstOrderDegree_sidechain = degree(influencenet_directed,mode=c("in"))

secondOrderDegree = c(); secondOrderDegree_sidechain=c()
nodes = V(influencenet)$name
for (node in nodes) {
  firstorder = neighbors(influencenet,node)
  secondorder = c(); 
  for (neighbor in firstorder){secondorder = c(secondorder,names(neighbors(influencenet,neighbor)))}; 
  secondorder = unique(secondorder); 
  secondorder=secondorder[secondorder!=node]
  secondOrderDegree = c(secondOrderDegree,length(secondorder))
  
  if(node %in% V(influencenet_directed)$name) {
    firstorder_sidechain = neighbors(influencenet_directed,node,mode=c("in"))
  } else {
    firstorder_sidechain = c()
  }
  secondorder = c(); 
  for (neighbor in firstorder_sidechain){
    secondorder = c(secondorder,names(neighbors(influencenet_directed,neighbor,mode=c("in"))))
  }; 
  secondorder = unique(secondorder); 
  secondorder=secondorder[secondorder!=node]
  secondOrderDegree_sidechain = c(secondOrderDegree_sidechain,length(secondorder))
}
names(secondOrderDegree) = nodes; names(secondOrderDegree_sidechain) = nodes


#second order intermodular degree
secondOrder_node_intermodular_degree = rep(0,length(node_intermodular_degree)); 
names(secondOrder_node_intermodular_degree) = names(node_intermodular_degree)

secondOrder_node_intermodular_degree_sidechain = rep(0,length(node_intermodular_degree_sidechain)); 
names(secondOrder_node_intermodular_degree_sidechain) = names(node_intermodular_degree_sidechain)

for (node in names(node_intermodular_degree)) {
  if (node_intermodular_degree[node]==0) {next}
  secondOrder_node_intermodular_degree[node] = length(unlist(node_modules[names(node_modules[[node]])]))
  secondOrder_node_intermodular_degree_sidechain[node] = length(unlist(node_modules_sidechain[names(node_modules_sidechain[[node]])]))
}


#SAVE ALL AS WALKTRAP
node_edge_betweenness_wt = node_edge_betweenness
node_edge_betweenness_sidechain_wt = node_edge_betweenness_sidechain
node_intermodular_degree_wt = node_intermodular_degree
node_intermodular_degree_sidechain_wt = node_intermodular_degree_sidechain
secondOrder_node_intermodular_degree_wt = secondOrder_node_intermodular_degree
secondOrder_node_intermodular_degree_sidechain_wt = secondOrder_node_intermodular_degree_sidechain


##2ARY STRUCTURE
secStructure = secStructure[,1:2]
rownames(secStructure) = secStructure[,1]
net_community_vec = as.numeric(as.factor(secStructure[,2]))
names(net_community_vec) = secStructure[,1]
#net_community_vec[buriednodes] = max(net_community_vec)+1
nodes = V(influencenet)$name

edgecolors = rep("grey90",nrow(basedata_noPP))
colorPalette=rainbow(max(net_community_vec),s=.5)
nodes.color = c()
for (i in 1:length(nodes)) {
  nodes.color = c(nodes.color,colorPalette[net_community_vec[nodes[i]]])
}
for (i in 1:nrow(basedata_noPP)) {
  if (net_community_vec[basedata_noPP[i,1]]!=net_community_vec[basedata_noPP[i,2]]) { edgecolors[i]="grey40" }
}

edge_betweenness = edge.betweenness(influencenet,weights=influencenet$weight)
node_edge_betweenness = rep(0,length(nodes))
names(node_edge_betweenness) = nodes
for(i in 1:length(edge_betweenness)) {
  if (net_community_vec[influencedata[i,1]]!=net_community_vec[influencedata[i,2]]) {
    node_edge_betweenness[influencedata[i,1]] = node_edge_betweenness[influencedata[i,1]]+edge_betweenness[i]
    node_edge_betweenness[influencedata[i,2]] = node_edge_betweenness[influencedata[i,2]]+edge_betweenness[i]
  }
}
node_intermodular_degree = rep(0,length(nodes))
names(node_intermodular_degree) = nodes
node_modules = list()
for (node in nodes) {
  subdata = matrix(influencedata[influencedata[,1]==node | influencedata[,2]==node,],ncol=4)
  if (nrow(subdata)==0) {node_intermodular_degree[node] = 0; next}
  bound = unique(c(subdata[,1:2]))
  bound = bound[bound!=node]
  bound_modules = net_community_vec[bound]
  bound_modules = bound_modules[bound_modules!=net_community_vec[node]]
  node_modules[[node]] = bound_modules
  node_intermodular_degree[node] = length(bound_modules)
}

##Weighted edge betweenness WEIGHT BY SIDE CHAIN
node_edge_betweenness_sidechain = rep(0,length(nodes))
names(node_edge_betweenness_sidechain) = nodes
for(i in 1:length(edge_betweenness)) {
  if (net_community_vec[influencedata[i,1]]!=net_community_vec[influencedata[i,2]]) {
    weight1 = as.numeric(influencedata_directed[influencedata_directed[,1]==influencedata[i,1]&influencedata_directed[,2]==influencedata[i,2],4]);if(length(weight1)==0){weight1=0}
    weight2 = as.numeric(influencedata_directed[influencedata_directed[,1]==influencedata[i,2]&influencedata_directed[,2]==influencedata[i,1],4]);if(length(weight2)==0){weight2=0}
    if(weight1==0 & weight2==0) {
      node_edge_betweenness_sidechain[influencedata[i,1]] = node_edge_betweenness_sidechain[influencedata[i,1]] + 0
      node_edge_betweenness_sidechain[influencedata[i,2]] = node_edge_betweenness_sidechain[influencedata[i,2]] + 0
    } else {
      node_edge_betweenness_sidechain[influencedata[i,1]] = node_edge_betweenness_sidechain[influencedata[i,1]]+(weight1/(weight1+weight2))*edge_betweenness[i]
      node_edge_betweenness_sidechain[influencedata[i,2]] = node_edge_betweenness_sidechain[influencedata[i,2]]+(weight2/(weight1+weight2))*edge_betweenness[i]
    }
  }
}
node_intermodular_degree_sidechain = rep(0,length(nodes))
names(node_intermodular_degree_sidechain) = nodes
node_modules_sidechain = list()
for (node in nodes) {
  subdata_sidechain = matrix(influencedata_directed[influencedata_directed[,1]==node,],ncol=4)
  if (nrow(subdata_sidechain)==0) {node_intermodular_degree_sidechain[node] = 0; next}
  bound = unique(c(subdata_sidechain[,1:2])); bound = bound[bound!=node]
  bound_modules = net_community_vec[bound]; bound_modules = bound_modules[bound_modules!=net_community_vec[node]]
  node_modules_sidechain[[node]] = bound_modules
  node_intermodular_degree_sidechain[node] = length(bound_modules)
}


##Degree and second order degree for 2ARY STRUCTURE

#second order intermodular degree
secondOrder_node_intermodular_degree = rep(0,length(node_intermodular_degree)); 
names(secondOrder_node_intermodular_degree) = names(node_intermodular_degree)
secondOrder_node_intermodular_degree_sidechain = rep(0,length(node_intermodular_degree_sidechain));
names(secondOrder_node_intermodular_degree_sidechain) = names(node_intermodular_degree_sidechain)
for (node in names(node_intermodular_degree)) {
  if (node_intermodular_degree[node]==0) {next}
  #secondOrder_node_intermodular_degree[node] = length(unique(unlist(node_modules[names(node_modules[[node]])])))
  secondOrder_node_intermodular_degree[node] = length(unlist(node_modules[names(node_modules[[node]])]))
  secondOrder_node_intermodular_degree_sidechain[node] = length(unlist(node_modules_sidechain[names(node_modules_sidechain[[node]])]))
}

#SAVE ALL AS STRIDE
node_edge_betweenness_stride = node_edge_betweenness
node_edge_betweenness_sidechain_stride = node_edge_betweenness_sidechain
node_intermodular_degree_stride = node_intermodular_degree
node_intermodular_degree_sidechain_stride = node_intermodular_degree_sidechain
secondOrder_node_intermodular_degree_stride = secondOrder_node_intermodular_degree
secondOrder_node_intermodular_degree_sidechain_stride = secondOrder_node_intermodular_degree_sidechain

### Create final dataset
nodes = nodes[nodes!="DNA1000A"]
nodes = nodes[nodes %in% c("DA","DG","DC","DT") == FALSE]

### Remove ligand amino acids from final dataset
ligands = sub("-","",liganddat)
nodes = nodes[nodes %in% ligands == FALSE]

nodes[which(nodes %in% rownames(bfactor)==FALSE)]

out = cbind(
  rsa[nodes,1],
  firstOrderDegree[nodes],
  firstOrderDegree_sidechain[nodes],
  secondOrderDegree[nodes],
  secondOrderDegree_sidechain[nodes],
  node_edge_betweenness_stride[nodes],
  node_edge_betweenness_sidechain_stride[nodes],
  node_intermodular_degree_stride[nodes],
  node_intermodular_degree_sidechain_stride[nodes],
  
  secondOrder_node_intermodular_degree_stride[nodes],
  secondOrder_node_intermodular_degree_sidechain_stride[nodes],
  node_edge_betweenness_wt[nodes],
  node_edge_betweenness_sidechain_wt[nodes],
  node_intermodular_degree_wt[nodes],
  node_intermodular_degree_sidechain_wt[nodes],
  secondOrder_node_intermodular_degree_wt[nodes],
  secondOrder_node_intermodular_degree_sidechain_wt[nodes]
)

#write column names
colnames(out) = c("RSA","Degree","Degree_uniqueAtoms","SecondOrderDegree","SecondOrderDegree_uniqueAtoms","NodeEdgeBetweennessSTRIDE","NodeEdgeBetweennessSTRIDE_unqiueAtoms","IntermodularDegreeSTRIDE","IntermodularDegreeSTRIDE_uniqueAtoms","SecondOrderIntermodularDegreeSTRIDE","SecondOrderIntermodularDegreeSTRIDE_uniqueAtoms","NodeEdgeBetweennessWALKTRAP","NodeEdgeBetweennessWALKTRAP_uniqueAtoms","IntermodularDegreeWALKTRAP","IntermodularDegreeWALKTRAP_uniqueAtoms","SecondOrderIntermodularDegreeWALKTRAP","SecondOrderIntermodularDegreeWALKTRAP_uniqueAtoms")


## LIGANDS
if(ligandmode==TRUE & ligandCentroidMode==FALSE) {
  liganddat = as.matrix(read.table(dir("./","_ligand"),sep="\t"))
  
  tmp = c()
  i=1
  for (i in 1:nrow(liganddat)) {
    acid=paste(unlist(strsplit(liganddat[i,1],split=""))[1:3],collapse="")
    if (unlist(strsplit(liganddat[i,1],"-"))[2] %in% terminalAtoms[[acid]]) {
      tmp = rbind(tmp,c(unlist(strsplit(liganddat[i,1],"-"))[1],liganddat[i,2],liganddat[i,3]))
    }
  }
  
  #ligandvec = tapply(liganddat[,3],as.factor(liganddat[,1]),min)
  ligandvec = tapply(as.numeric(tmp[,3]),as.factor(tmp[,1]),min)
}

if(ligandmode==TRUE & ligandCentroidMode) {
  liganddat = as.matrix(read.table(dir("./",paste("_centroidNetLigand$",sep="")),sep="\t"))
  
  tmp = c()
  i=1
  for (i in 1:nrow(liganddat)) {
    tmp = rbind(tmp,c(liganddat[i,1],liganddat[i,2],liganddat[i,3]))
  }
  
  ligandvec = tapply(tmp[,3],as.factor(tmp[,1]),min)
}

if (ligandmode==FALSE) {
  ligandvec = rep(0,length(nodes))
}

out = cbind(out,ligandvec[rownames(out)])
colnames(out) = c("RSA","Degree","Degree_uniqueAtoms","SecondOrderDegree","SecondOrderDegree_uniqueAtoms","NodeEdgeBetweennessSTRIDE","NodeEdgeBetweennessSTRIDE_unqiueAtoms","IntermodularDegreeSTRIDE","IntermodularDegreeSTRIDE_uniqueAtoms","SecondOrderIntermodularDegreeSTRIDE","SecondOrderIntermodularDegreeSTRIDE_uniqueAtoms","NodeEdgeBetweennessWALKTRAP","NodeEdgeBetweennessWALKTRAP_uniqueAtoms","IntermodularDegreeWALKTRAP","IntermodularDegreeWALKTRAP_uniqueAtoms","SecondOrderIntermodularDegreeWALKTRAP","SecondOrderIntermodularDegreeWALKTRAP_uniqueAtoms","Ligand")

#add back in nodes that were not networked as zeros
zeroMat = matrix(0,nrow=length(which(originalNodes%in%nodes==FALSE)),ncol=ncol(out));
rownames(zeroMat)=originalNodes[which(originalNodes%in%nodes==FALSE)]
out = rbind(out,zeroMat)

#fix RSA and Bfactor and ligand
out[,"RSA"] = rsa[rownames(out),1]
out[,"Ligand"] = ligandvec[rownames(out)]

#Add 1 to everything to avoid zeros, except ligand
out[is.na(out)]=0
out[out[,"Ligand"]==0,"Ligand"] = 150
out[1:nrow(out),2:ncol(out)] = as.numeric(out[1:nrow(out),2:ncol(out)])+1

##create Z score file
outZ = apply(out,2,function(x){return(scale(as.numeric(x)))})
rownames(outZ) = rownames(out)
colnames(outZ) = colnames(out)
outZ[outZ==NA]=0
outZ[outZ=="NaN"]=0

#Write out scores and z-normed scores files
write.table(out,file=paste(keyword,"_scoresCentroid",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)
write.table(outZ,file=paste(keyword,"_scoresCentroidZ",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

# check with python code 
for (col in colnames(out)) {print(sum(as.numeric(out[,col])))}
