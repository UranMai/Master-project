args <- commandArgs(trailingOnly = TRUE)
# input file, ex. 1BTL.pdb, 4eiy.pdb

##Load libraries
library(igraph)
library(stats)

args = '4eiy' #for running with source in rstudio
##Initialize keywords
# keyword is the base filename
keyword = args[1] # keyword = '4eiy'
# change args[1] to keyword
args=c(paste(keyword,"_net",sep=""),"weighted",paste(keyword,"_secondaryStructure",sep="")) # ['4eiy_net', 'weighted', '4eiy_secondaryStructure']
atomCorrection=FALSE
baseKeyword = unlist(strsplit(keyword,"/"))[1] # 4eiy
runVariableSelection = TRUE

##Load inversal data
#Unique atoms file, some spaces can be problem for search of terminal atoms
terminalAtomsFile = as.matrix(read.table("./terminalAtoms",sep="\t"))
terminalAtoms = {}
for (i in 1:nrow(terminalAtomsFile)) {
  terminalAtoms[[terminalAtomsFile[i,1]]] = unlist(strsplit(terminalAtomsFile[i,2],split=", "))
}

# nrow(secStructure), sum(as.numeric(secstructure[,3]))
#Secondary structure file
secStructure = as.matrix(read.table(args[3]))[,1:2]
rownames(secStructure) = secStructure[,1]

#RSA
rsa = as.matrix(read.table(paste(keyword,".rsa",sep="")))

##Define functions

#Add terminal atom details
addTermDetails = function(data) {
  out_details = c()
  for (i in 1:nrow(data)) {
    acid1 = paste(unlist(strsplit(data[i,1],split=""))[1:3],collapse="")
    acid2 = paste(unlist(strsplit(data[i,2],split=""))[1:3],collapse="")
    if (data[i,6] %in% terminalAtoms[[acid1]]) {
      code1 = 1
    } else {
      code1 = 0
    }
    if (data[i,7] %in% terminalAtoms[[acid2]]) {
      code2 = 1
    } else {
      code2 = 0
    }
    out_details = rbind(out_details, c(data[i,],code1,code2))
  }
  return(out_details)
}

collapse = function(data) {
  if (nrow(data)==0) {
    out = matrix(ncol = 8,nrow=0)
    colnames(out) = c("AA1", "AA2", "Type", "SumEdges", "percentMC", "influence1", "influence2","degree")
    results = c(); results$out = out
    return(results)
    break
  }
  #Edges should be listed bi-directionally
  x=unlist(lapply(data[,3],function(x){if(x=="MCSC"){return("SCMC")}else if(x=="SCMC"){return("MCSC")} else {return(x)}}))
  if (dim(data)[2]==7) {
    data = rbind(data,cbind(data[,2],data[,1],x,data[,4],data[,5],data[,7],data[,6]))
  } else {
    data = rbind(data,cbind(data[,2],data[,1],x,data[,4],data[,5],data[,7],data[,6],data[,9],data[,8]))
  }
  data = unique(data)
  out = c()
  out_noPP = c()
  completed = c()
  if (dim(data)[2]==7) {
    out_details = c()
    for (i in 1:nrow(data)) {
      acid1 = paste(unlist(strsplit(data[i,1],split=""))[1:3],collapse="")
      acid2 = paste(unlist(strsplit(data[i,2],split=""))[1:3],collapse="")
      if (data[i,6] %in% terminalAtoms[[acid1]]) {
        code1 = 1
      } else {
        code1 = 0
      }
      if (data[i,7] %in% terminalAtoms[[acid2]]) {
        code2 = 1
      } else {
        code2 = 0
      }
      out_details = rbind(out_details, c(data[i,],code1,code2))
    }
    data = out_details
  } else {
    out_details = data
  }
  colnames(out_details) = c("AA1","AA2","Type1","Weight","Type2","Atom1","Atom2","influence1","influence2")
  for (i in 1:nrow(data)) {
    if (i %in% completed == FALSE) {
      x=which((data[,1]==data[i,1] & data[,2] == data[i,2]))
      subdata = matrix(data[x,],ncol=9)
      edgeSum = 0
      MCMCcount = 0
      MCSCcount = 0
      SCMCcount = 0
      SCSCcount = 0
      influence1 = sum(as.numeric(subdata[,8]))
      influence2 = sum(as.numeric(subdata[,9]))
      for (j in 1:nrow(subdata)) {
        edgeSum = edgeSum + as.numeric(subdata[j,4])
        if (subdata[j,3]=="MCMC" | subdata[j,3]=="MCSC") {
          MCMCcount = MCMCcount+1
          MCSCcount = MCSCcount+1
        } else if (subdata[j,3]=="SCMC") {
          SCMCcount = SCMCcount+1
        } else if (subdata[j,3]=="SCSC") {
          SCSCcount = SCSCcount+1
        }
      }
      #Header is AA1, AA2, mixed, sum of edges, percent MC, influence1, influence2, degreeAA1
      out = rbind(out,c(data[i,1],data[i,2],"mixed",edgeSum,(MCMCcount+MCSCcount)/(MCMCcount+MCSCcount+SCMCcount+SCSCcount),influence1,influence2,sum(data[,1]==data[i,1])))
      if (any(subdata[,5]=="PP")==FALSE) {
        out_noPP = rbind(out_noPP,c(data[i,1],data[i,2],"mixed",edgeSum,(MCMCcount+MCSCcount)/(MCMCcount+MCSCcount+SCMCcount+SCSCcount),influence1,influence2,sum(data[,1]==data[i,1])))}
      completed = c(completed,x)
    }
  }
  out = out[out[,4]>0,]
  out_noPP = out_noPP[out_noPP[,4]>0,]
  colnames(out) = c("AA1", "AA2", "Type", "SumEdges", "percentMC", "influence1", "influence2","degree")
  results = c()
  results$out = out
  results$out_noPP = out_noPP
  results$influence_all = out_details[as.numeric(out_details[,4])>0 & (out_details[,8]=="1" | out_details[,9]=="1"),]
  out_details = out_details[as.numeric(out_details[,4])>0,]
  out_subtract = out_details; out_subtract[out_subtract[,8]==0 & out_subtract[,9]==0,4] = (-1)*as.numeric(out_subtract[out_subtract[,8]==0 & out_subtract[,9]==0,4])
  out_subtract = out_subtract[out_subtract[,5]!="PP",]
  results$out_subtract = out_subtract
  return(results)
}

removeRedundancy = function(dataTmp) {
  for (i in 1:nrow(dataTmp)) {
    if (order(dataTmp[i,1:2])[1]==2) {
      dataTmp[i,1:2] = c(dataTmp[i,2],dataTmp[i,1])
      dataTmp[i,3] = paste(unlist(strsplit(dataTmp[i,3],split=""))[c(3,4,1,2)],collapse="")
      dataTmp[i,c(6,7)] = c(dataTmp[i,7],dataTmp[i,6])
    }
  }
  dataTmp = unique(dataTmp)
  pairs = unique(dataTmp[,1:2])
  out = c()
  for (j in 1:nrow(pairs)) {
    pair = pairs[j,]
    subData = matrix(dataTmp[dataTmp[,1]==pair[1] & dataTmp[,2]==pair[2],],ncol=7)
    out = rbind(out, c(pair[1],pair[2],"mixed",0.5*sum(as.numeric(subData[,4]))))
  }
  return(out)
}

collapse_directed = function(data) {
  pairs = unique(data[,1:2])
  out = c()
  for (i in 1:nrow(pairs)) {
    subdat = matrix(data[data[,1]==pairs[i,1]&data[,2]==pairs[i,2],],ncol=4)
    out = rbind(out,c(subdat[1,1],subdat[1,2],"mixed",sum(as.numeric(subdat[,4]))))
  }
  return(out)
}

collapse_agnostic = function(data) {
  pairs = unique(t(apply(data[,1:2],1,function(x){return(sort(x))})))
  out = c()
  for (i in 1:nrow(pairs)) {
    subdat = matrix(data[(data[,1]==pairs[i,1]&data[,2]==pairs[i,2]) | (data[,2]==pairs[i,1]&data[,1]==pairs[i,2]),],ncol=4)
    out = rbind(out,c(subdat[1,1],subdat[1,2],"mixed",sum(as.numeric(subdat[,4]))))
  }
  return(out)
}

distributCalc = function(set,nodes,net) {
  #x is assumed to be a set but can also just be 1 node
  distanceNet = distances(net,weights=net$weight)
  distanceTmp = c()
  Dr = c()
  for (node in nodes) {
    distanceTmp = c()
    for (item in set) {
      if (item!=node) {
        if(item%in%rownames(distanceNet)==FALSE){distanceTmp=100;next}
        if(node%in%rownames(distanceNet)==FALSE){distanceTmp=100;next}
        distanceTmp = c(distanceTmp,distanceNet[item,node])
      }
    }
    if (length(distanceTmp)>0) {
      Dr = c(Dr, 1/min(distanceTmp))
    }
  }
  return(sum(Dr)/length(nodes))
}

##Read in user data, read 1BTL_net file 
data_all = as.matrix(read.table(args[1]))
nodes = unique(c(data_all[,1:2])); sets = nodes # all aminoacids define
originalNodes = nodes
basedata = collapse(data_all)
basedata_noPP = removeRedundancy(basedata$out[,1:7])
basedata_allsidechain = collapse_agnostic(data_all[data_all[,3]!="MCMC",c(1,2,3,4)])
influencedata_all = basedata$influence_all
influencedata_collapse = collapse(influencedata_all)
influencedata = removeRedundancy(influencedata_collapse$out[,1:7])
# Create graph network
influencenet = graph.edgelist(influencedata[,1:2],directed=FALSE)
influencenet$weight = 1/as.numeric(influencedata[,4])
influencedata_directed = collapse_directed(influencedata_all[influencedata_all[,8]==1,c(1,2,3,4)])
influencenet_directed = graph.edgelist(influencedata_directed[,c(2,1)],directed=TRUE)
influencenet_directed$weight = 1/as.numeric(influencedata_directed[,4])

nodes = unique(c(data_all[,1:2])); sets = nodes

#Write out edgelist files
#write.table(influencedata,file=paste(keyword,"_allEdges",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
#write.table(influencedata_directed,file=paste(keyword,"_uniqueEdges",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

###Deal with missing nodes
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
}

##Begin centrality calculations

#Modular calculations

##WALKTRAP
net = graph.edgelist(basedata_noPP[,1:2],directed=FALSE)
net$weight = 1/as.numeric(basedata_noPP[,4])
net_community=walktrap.community(net,weights=net$weight); net_community_vec = net_community$membership; names(net_community_vec)=net_community$names
net_community_vec_wt = net_community_vec
nodes = V(net)$name

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
# if(length(nodes)>200) {textsize=.25} else {textsize=.45}
# plot(net,edge.color = edgecolors, vertex.color = nodes.color,vertex.size=4,vertex.label.cex=textsize,main="WALKTRAP")

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
  secondorder = c(); for (neighbor in firstorder){secondorder = c(secondorder,names(neighbors(influencenet,neighbor)))}; secondorder = unique(secondorder); secondorder=secondorder[secondorder!=node]
  secondOrderDegree = c(secondOrderDegree,length(secondorder))
  firstorder = neighbors(influencenet_directed,node,mode=c("in"))
  secondorder = c(); for (neighbor in firstorder){secondorder = c(secondorder,names(neighbors(influencenet_directed,neighbor,mode=c("in"))))}; secondorder = unique(secondorder); secondorder=secondorder[secondorder!=node]
  secondOrderDegree_sidechain = c(secondOrderDegree_sidechain,length(secondorder))
}
names(secondOrderDegree) = nodes; names(secondOrderDegree_sidechain) = nodes
secondOrder_node_intermodular_degree = rep(0,length(node_intermodular_degree)); names(secondOrder_node_intermodular_degree) = names(node_intermodular_degree)
secondOrder_node_intermodular_degree_sidechain = rep(0,length(node_intermodular_degree_sidechain)); names(secondOrder_node_intermodular_degree_sidechain) = names(node_intermodular_degree_sidechain)
for (node in names(node_intermodular_degree)) {
  if (node_intermodular_degree[node]==0) {next}
  secondOrder_node_intermodular_degree[node] = length(unlist(node_modules[names(node_modules[[node]])])) - 1
  secondOrder_node_intermodular_degree_sidechain[node] = length(unlist(node_modules_sidechain[names(node_modules_sidechain[[node]])])) - 1
}

#SAVE ALL AS WALKTRAP
node_edge_betweenness_wt = node_edge_betweenness
node_edge_betweenness_sidechain_wt = node_edge_betweenness_sidechain
node_intermodular_degree_wt = node_intermodular_degree
node_intermodular_degree_sidechain_wt = node_intermodular_degree_sidechain
secondOrder_node_intermodular_degree_wt = secondOrder_node_intermodular_degree
secondOrder_node_intermodular_degree_sidechain_wt = secondOrder_node_intermodular_degree_sidechain

##2ARY STRUCTURE
net_community_vec = as.numeric(as.factor(secStructure[,2]))
names(net_community_vec) = secStructure[,1]
nodes = V(net)$name
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
    node_edge_betweenness_sidechain[influencedata[i,1]] = node_edge_betweenness_sidechain[influencedata[i,1]]+(weight1/(weight1+weight2))*edge_betweenness[i]
    node_edge_betweenness_sidechain[influencedata[i,2]] = node_edge_betweenness_sidechain[influencedata[i,2]]+(weight2/(weight1+weight2))*edge_betweenness[i]
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
nodes = V(influencenet)$name
secondOrder_node_intermodular_degree = rep(0,length(node_intermodular_degree)); names(secondOrder_node_intermodular_degree) = names(node_intermodular_degree)
secondOrder_node_intermodular_degree_sidechain = rep(0,length(node_intermodular_degree_sidechain)); names(secondOrder_node_intermodular_degree_sidechain) = names(node_intermodular_degree_sidechain)
for (node in names(node_intermodular_degree)) {
  if (node_intermodular_degree[node]==0) {next}
  secondOrder_node_intermodular_degree[node] = length(unlist(node_modules[names(node_modules[[node]])])) - 1
  secondOrder_node_intermodular_degree_sidechain[node] = length(unlist(node_modules_sidechain[names(node_modules_sidechain[[node]])])) - 1
}


#SAVE ALL AS STRIDE
node_edge_betweenness_stride = node_edge_betweenness
node_edge_betweenness_sidechain_stride = node_edge_betweenness_sidechain
node_intermodular_degree_stride = node_intermodular_degree
node_intermodular_degree_sidechain_stride = node_intermodular_degree_sidechain
secondOrder_node_intermodular_degree_stride = secondOrder_node_intermodular_degree
secondOrder_node_intermodular_degree_sidechain_stride = secondOrder_node_intermodular_degree_sidechain


### Create final dataset

out = cbind(
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
rownames(out) = nodes
colnames(out) = c("Degree","Degree_sidechain","SecondOrderDegree","SecondOrderDegree_sidechain","NodeEdgeBetweennessSTRIDE","NodeEdgeBetweennessSTRIDE_sidechain","IntermodularDegreeSTRIDE","IntermodularDegreeSTRIDE_sidechain","SecondOrderIntermodularDegreeSTRIDE","SecondOrderIntermodularDegreeSTRIDE_sidechain","NodeEdgeBetweennessWALKTRAP","NodeEdgeBetweennessWALKTRAP_sidechain","IntermodularDegreeWALKTRAP","IntermodularDegreeWALKTRAP_sidechain","SecondOrderIntermodularDegreeWALKTRAP","SecondOrderIntermodularDegreeWALKTRAP_sidechain")

#add back in nodes that were not networked as zeros
zeroMat = matrix(0,nrow=length(which(originalNodes%in%nodes==FALSE)),ncol=ncol(out));rownames(zeroMat)=originalNodes[which(originalNodes%in%nodes==FALSE)]
out = rbind(out,zeroMat)

#Add 1 to everything to avoid zeros, except ligand
out[is.na(out)]=0
out = out+1

#Write out
#write.table(out,file=paste(keyword,"_scoresEnergetics",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

#Standard normalization
outZ = apply(out,2,function(x){return(scale(x))})
rownames(outZ) = rownames(out)
#write.table(outZ,file=paste(keyword,"_scoresEnergeticsZ",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)

for (col in colnames(out)) {print(sum(as.numeric(out[,col])))}
