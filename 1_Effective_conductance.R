library(ggplot2)
library(dominatR)
require(tidyverse)
require(data.table)
library(igraph)
require(dplyr)
require(MASS)
require(entropy)
require(igraph)
# require(dominatR)

eval(parse(file="./CODES/my_function.R"))
eval(parse(file="./CODES/plotTriangle_ver2.R"))

#First download BIOGRID Data
url <- "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.243/BIOGRID-ALL-4.4.243.tab3.zip"
download.file(url, destfile = "BIOGRID-ALL-4.4.243.tab3.zip", mode = "wb")


########################################################################
########################################################################
#########################################Making PPI network from BIOGRID
########################################################################
########################################################################
{


ppi=fread("BIOGRID-ALL-4.4.243.tab3.zip",sep="\t") 
organism=ppi %>% filter(`Experimental System Type`=="physical") %>% group_by(`Organism ID Interactor A`,`Organism Name Interactor A`) %>% summarize(number=n()) %>% ungroup() %>% arrange(desc(number))

ppi2=ppi %>% dplyr::transmute(A=`Official Symbol Interactor A`,B=`Official Symbol Interactor B`,
                              Az=`Synonyms Interactor A`,Bz=`Synonyms Interactor B`,
                              ORGa=`Organism ID Interactor A`,ORGb=`Organism ID Interactor B`,
                              Score,Throughput,Exp=`Experimental System`,Exptype=`Experimental System Type`,
                              Author,pub=`Publication Source`)
ppi2=ppi2 %>% filter(ORGa==9606 & ORGb==9606) %>% dplyr::select(-c(ORGa,ORGb)) #Being Human
# ppi2=ppi2 %>% filter(ORGa==7227 & ORGb==7227) %>% dplyr::select(-c(ORGa,ORGb)) #Being something else

ppi2=ppi2 %>% filter(Exptype == "physical") #Being physical (and not genetical)
table(ppi2$Exptype)

########################################################
#############measuring experiment_number bias
########################################################
{expbias=as.data.frame(table(ppi2$A),responseName = "A")
reality=as.data.frame(table(ppi2$B),responseName = "B")
bias=full_join(expbias,reality) 
bias[is.na(bias)]=0
bias=bias %>% mutate(diff=A-B,Var1=as.character(Var1)) 
bias= bias %>% rowwise() %>% mutate(num_author=length(unique(ppi2$Author[ppi2$A==Var1]))) #might take 5 minute or so to process this line
saveRDS(bias,"bias.Rdata")
bias=readRDS("bias.Rdata")

}##Luka##




ppi2=ppi2 %>% filter(A!=B) # removing the interaction of a protein to itself
table(ppi2$Exp);paste(sort(unique(ppi2$Exp)),collapse=", ")

ppi3= rbind(ppi2 %>% transmute(A,B,C=1),
            ppi2 %>% transmute(C=A,A=B,B=C,C=1)) # reversing the interaction and concatenating it to the original interaction
ppi3=ppi3 %>% group_by(A,B) %>% summarize(num=length(C)) %>% ungroup() # ppi contains 2x interaction, coz there will be A->B and B->A interaction, and both having same score ( unidirectional)
link=data.frame(prot=(unique(c(ppi3$A,ppi3$B)))) %>% mutate(ID=1:nrow(.)) ##labelling proteins
ppi3=left_join(ppi3,link %>% transmute(A=prot,Ai=ID))
ppi3=left_join(ppi3,link %>% transmute(B=prot,Bi=ID))
ppi3=ppi3 %>% filter(Ai!=Bi) # removing the interaction of a protein to itself

matrix= matrix(0,nrow=nrow(link),ncol=nrow(link),dimnames = list(link$prot, link$prot)) #making a null square matrix
#Now filling up the null matrix with interaction (both upper and lower part of the diagonal. since our ppi3 has 2x interaction, 1x of upper and another 1x of lower triangle)
for(i in (1:nrow(ppi3)))
{
  if(i%%1000==0) cat(paste0("\r ",i/1000,"k/",nrow(ppi3)/1000,"k"))
  matrix[ppi3$Ai[i],ppi3$Bi[i]]=ppi3$num[i]
}

########################################################
##############Removing the disjoint networks and nodes if any
########################################################
ifelse(sum(rownames(matrix)!=colnames(matrix))==0," Good, row and columns are aligned", "Stop!! Row and columns are not aligned ")  #making sure the proteins are arranged same way in both row and columns
g <- graph_from_adjacency_matrix(matrix, mode = "undirected", diag = FALSE) # Convert adjacency matrix to an igraph graph object
disconnected_components <- decompose(g) # Find disconnected components
vertices_per_component <- lapply(disconnected_components, function(component) {V(component)$name}) # Extract vertex names from each disconnected component
# Print vertices for each disconnected component
print(paste0("Number of disjoint networks = ",length(vertices_per_component)))
#Removing the smaller disjoint networks if there are any. Manually verify that vertices_per_component contains the biggest network in its first element.
if(length(vertices_per_component)>1)
{
  disconnected_proteins=unlist(sapply(2:length(vertices_per_component),function(x) vertices_per_component[[x]]))
  indx=!rownames(matrix) %in% disconnected_proteins
  matrix=matrix[indx,indx]
}


saveRDS(matrix,"interaction_matrix.Rdata")

}##Luka##


########################################################################
########################################################################
#######################################Calculating Effective Conductance
########################################################################
########################################################################
{

#########################Finding Admittance(A) and degree of graph(D)
A=matrix
diag(A)=0
D=A*0
diag(D)=rowSums(A)
#########################Finding Laplacian
L=D-A  #Computing Normal Laplacian
# D_inv_sqrt <- diag(1 / sqrt(diag(D)))  # Compute D^(-1/2) for normalized laplacian
# L <- D_inv_sqrt %*% L %*% D_inv_sqrt  # Computing normalized laplacian


#########################Finding PseudoInverse
z=ginv(L)
saveRDS(z,"pseudoinversed.Rdata")

#########################Using the formula of Effective conductance
temp=z
for (j in 1:nrow(temp)) 
{ if(j%%100==0) cat(paste0("\r ",j,"/",nrow(temp)))
  for (k in 1:nrow(temp)) {temp[j, k] = z[j, j] + z[k, k] - z[j, k] - z[k, j]}}
EC=1/temp
colnames(EC)=colnames(matrix);rownames(EC)=rownames(matrix)
saveRDS(EC,"EC.Rdata") 
EC=readRDS("EC.Rdata")

}##Luka##


########################################################################
########################################################################
#####################################################Reciprocal analysis
########################################################################
########################################################################
{
#Z-scaling the rows
tem=EC
tem[is.infinite(tem)]=0
for(i in c(1:nrow(tem)))
{
  if(i%%1000==0) cat(paste0("\r ",i,"/",nrow(tem)))
  tem[i,]=scale(tem[i,])
}
saveRDS(tem,"tem.Rdata")
tem=readRDS("tem.Rdata")


##Plot Reciprocal plot
hmm=gimme_EC(prot="POLR3H",threshold=0.01,neighbourhood=T,improved_fit = T,df_improvedfit = 15000,del=tem)


########Calculating Reciprocal plot for all proteins and Creating ECd,ECdz,ECdz_min matrix
del=tem
dim(del)
ECd=del*0
for(i in c(1:nrow(ECd)))
{
  current_protein=rownames(ECd)[i]
  if(i%%10==0) cat(paste0("\r ",i,"/",nrow(ECd)))
  temp=gimme_EC(prot=rownames(ECd)[i],plot = F,improved_fit = T,instructions = F,df_improvedfit = 18000,del = tem)
ECd[rownames(ECd)[i],rownames(temp)]=temp$distances
  }

saveRDS(ECd,"ECd.Rdata")
ECd=readRDS("ECd.Rdata")

####Creating ECdz
ECdz=ECd
for(i in c(1:nrow(ECdz)))
{
  if(i%%100==0) cat(paste0("\r ",i,"/",nrow(ECdz)))
  ECdz[i,]=scale(ECdz[i,])
}
ECdz=abs(ECdz)
####Finding ECdz_avg and ECdz_min
ECdz_avg=(ECdz+t(ECdz))/2
ECdz_min=pmin(ECdz, t(ECdz))
saveRDS(ECdz_min,"ECdz_min.Rdata")


}##Luka##







