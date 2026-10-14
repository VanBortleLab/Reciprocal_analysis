require(dplyr)
require(Matrix)


########################################################################
########################################################################
###################FIGURE1######################################################
##########How zz plot is applied in BIOGRID?? FLOWCHART!!!!!!!
########################################################################
########################################################################
{
  if(!exists("EC")) EC=readRDS('EC.Rdata')
  if(!exists("ECdz_min")) ECdz_min=readRDS('ECdz_min.Rdata')
  matrix=readRDS("interaction_matrix.Rdata")
matrix[1:10,1:10] #taken from Effective_conductance.R after removing the disjoint network
indx=names(sort(rowSums(matrix))) #names sorted as per interaction frequency

bias=readRDS("bias.Rdata") %>% mutate(sum=A+B) %>% arrange(num_author,sum)
indx=bias$Var1;indx=indx[indx %in% rownames(matrix)]
tempC=matrix[indx,indx]
tempS=ECdz_min[indx,indx]
tempEC=EC[indx,indx];diag(tempEC)=0


###############################################################
#########plotting the entire network (not a good idea)
###############################################################
{
if(T){jpeg("IMAGES/fig1_whole_.jpeg",width=1200,height = 1200,quality = 1200)} 
par(mfrow=c(1,3),pty="s")
par(pty="s")
#heatmap(tempEC,Rowv = NA, Colv = NA,col=colorRampPalette(c("white","red","red4"))(100),scale="none",labRow = F,labCol = F) # Use this if you wanna plot entire network
indx=sort(sample(size=5000,c(1:nrow(tempEC))))
Matrix::image(
  x=tempC[indx,indx],
  col=colorRampPalette(c("white","red","red4"))(250),
  main="Conductance (frequency of interaction)",
  useRaster=T,
  scales = list(draw = FALSE)
  )
ew=as.matrix(tempC[indx,indx])[nrow(tempC[indx,indx]):1, ] 
ew <- as(as(Matrix(ew, sparse = TRUE), "generalMatrix"), "TsparseMatrix")
Matrix::image(
  x = ew,
  lwd = 2,
  col.regions= colorRampPalette(c("white", "red", "red4"))(250),
  useRaster = T,
  main = "Conductance (frequency of interactions)",
  scales = list(draw = FALSE)
)

sapply(c(1:4), function(x) axis(x,labels = F,col.ticks = "transparent",lwd=2))
Matrix::image(tempEC[indx,indx],col=colorRampPalette(c("white","red","red4"))(250),axes=F,main="Effective Conductance (closeness)",useRaster=F)
sapply(c(1:4), function(x) axis(x,labels = F,col.ticks = "transparent",lwd=2))
image(tempS[indx,indx],col=colorRampPalette(c("white","red","red4"))(250),axes=F,main="Closeness Score (bias free)")
sapply(c(1:4), function(x) axis(x,labels = F,col.ticks = "transparent",lwd=2))
dev.off()
}##luka##


###############################################################
#########plotting the entire network using dots. (good idea)
###############################################################
{
tempC_sparse<- summary(Matrix(tempC, sparse = TRUE)) %>% dplyr::arrange(x) %>% mutate(alpha= pmin(x/max(x),1),
                                                                 x_scaled=as.numeric(cut(x, breaks = 10000, labels = FALSE)),
                                                                 x_log_scaled=as.numeric(cut(log10(x+100), breaks = 10000, labels = FALSE)))
tempEC_sparse<- summary(Matrix(tempEC, sparse = TRUE)) %>% dplyr::arrange(x)  %>% mutate(alpha= pmin(x/max(x),1),
                                                                   x_scaled=as.numeric(cut(x, breaks = 10000, labels = FALSE)),
                                                                   x_log_scaled=as.numeric(cut(log10(x+100), breaks = 10000, labels = FALSE)))
tempS_sparse<- summary(Matrix(tempS, sparse = TRUE))  %>% dplyr::arrange(x) %>% mutate(alpha= pmin(x/max(x),1),
                                                                 x_scaled=as.numeric(cut(x, breaks = 10000, labels = FALSE)),
                                                                 x_log_scaled=as.numeric(cut(log10(x+100), breaks = 10000, labels = FALSE)))
size=nrow(EC)
#####Plotting colorbars...
#######Same Colorbar is used for all matrix, but the scale is different for different matrix, for example for conductance the end of color bar represents 659, whereas for effective conductacnce the end is some value more than 1000.
# colors=colorRampPalette(c(rep("white",50),rep("orange",100),rep("red",100), rep("red2",500),rep("red3",500),rep("red4",1000)))(10000)
colors=colorRampPalette(c("white",rep("orange",1),rep("red",1),rep("red2",1),rep("red3",1),rep("red4",1)))(10000)
plot(c(1:10000),rep(1,10000),col=colors,pch=19);abline(v=1)
plot(c(1,10000),c(0,5),col="transparent",axes=F,xlab="",ylab="",main="Effective Conductance")
sapply(c(1:10000),function(i) points(c(i,i),c(0.5,2),type="l",lwd=0.4,col=colors[i]))
rect(1,0.5,10000,2,lwd=2)

ew=c(min(tempEC_sparse$x),max(tempEC_sparse$x)) #actual scores at two ends of colorbar
ew=log10(ew+100) #transformed scores at two ends of colorbar
coor=c(1,10000) #x coordinate positions of two end of color bar

gimme_coordinate=function(score,ew,coor) #derived from (y2-y2)/(x2-x1)=(y-y1)/(x-x1) where x is coordinates at color bar, and y is log transformed score,i.e. log10(100+score)
{enquired_coor=coor[1]+diff(coor) * (log10(100+score)-ew[1]) / (ew[2]-ew[1])
  return(enquired_coor)
}
gimme_coordinate(1000,ew,coor)

sapply(c(10,100,1000),function(i) 
  {
  i2=gimme_coordinate(i,ew,coor)
  points(c(i2,i2),c(0.4,0.5),type="l",lwd=2)
  text(i2,0.3,labels=i)
  }
  )
which(tempEC_sparse$x_scaled==1)

# jpeg("IMAGES/fig1_100percent_arrangedby_numauthor_simplelogcolorbar.jpeg",width=2000,height = 1000,quality = 1500)
par(mfrow=c(1,2),pty="s")


temp=tempC_sparse  %>% dplyr::slice_tail(n=as.integer(nrow(.)*0.05)) 
plot(temp$j, temp$i,pch = 20,cex = 1,col = colors[temp$x_log_scaled],  # red, alpha scaled by value
     xlab = "Protein index",ylab = "Protein index",main = "Conductance",
     axes=F,xlim=c(1,size),ylim=c(1,size))
# points(temp$i, temp$j,pch = 20,cex = 1,col = rgb(1, 0, 0,1*temp$alpha),  # red, alpha scaled by value
#        xlab = "Protein index",ylab = "Protein index",main = "ewew")
points(temp$i, temp$j,pch = 20,cex = 1,col = colors[temp$x_log_scaled])  # red, alpha scaled by value)
rect(xleft = 1, ybottom = 1, xright = size, ytop = size, border = "black", lwd = 4)



temp=tempEC_sparse %>% dplyr::slice_tail(n=as.integer(nrow(.)*0.005)) 
plot(temp$j, temp$i,pch = 20,cex = 1,col = colors[temp$x_log_scaled],  # red, alpha scaled by value
     xlab = "Protein index",ylab = "Protein index",main = "Effective conductance",
    axes=F,xlim=c(1,size),ylim=c(1,size))
points(temp$i, temp$j,pch = 20,cex = 1,col = colors[temp$x_log_scaled])  # red, alpha scaled by value)
rect(xleft = 1, ybottom = 1, xright = size, ytop = size, border = "black", lwd = 4)




temp=tempS_sparse %>% dplyr::slice_tail(n=as.integer(nrow(.)*0.05)) 
plot(temp$j, temp$i,pch = 20,cex = 1,col = colors[temp$x_scaled],  # red, alpha scaled by value
     xlab = "Protein index",ylab = "Protein index",main = "Closeness score",
     axes=F,xlim=c(1,size),ylim=c(1,size))
points(temp$i, temp$j,pch = 20,cex = 1,col = colors[temp$x_scaled])  # red, alpha scaled by value)
rect(xleft = 1, ybottom = 1, xright = size, ytop = size, border = "black", lwd = 4)
dev.off()

# plot_dotted_matrix(tempC[indx,indx],"Conductance matrix")
# plot_dotted_matrix(tempEC[indx,indx],"Effective Conductance matrix",sample_frac = 0.0001)
# tempSS=tempS;tempSS[tempSS<2]=0;plot_dotted_matrix(tempSS,"Closeness Score >2, sorted by sum,num_author")

}##luka##


###############################################################
#########plotting the entire network using reduced matrix
#################################obsolete !!!################
reduce_matrix <- function(mat, target_height = 1000, target_width = 1000, summary_func = max) 
  {
  if (target_height > nrow(mat) || target_width > ncol(mat)) {
    stop("Target size must be smaller than the original matrix")
  }
  
  # Divide matrix into chunks
  row_bins <- round(seq(1, nrow(mat) + 1, length.out = target_height + 1))
  col_bins <- round(seq(1, ncol(mat) + 1, length.out = target_width + 1))
  
  # Initialize reduced matrix
  reduced <- matrix(0, nrow = target_height, ncol = target_width)
  
  for (i in 1:target_height) {
    for (j in 1:target_width) {
      block <- mat[row_bins[i]:(row_bins[i+1] - 1), col_bins[j]:(col_bins[j+1] - 1)]
      reduced[i, j] <- summary_func(block, na.rm = TRUE)
    }
  }
  
  return(reduced)
}
tempC_reduced=reduce_matrix(tempC,target_height = 100,target_width = 100,summary_func = max)
image(tempC_reduced,col=colorRampPalette(c("white","red","red4"))(250),axes=F,main="Conductance, using reduced matrix")
sapply(c(1:4), function(x) axis(x,labels = F,col.ticks = "transparent",))

###############################################################
##########Plotting a small slice of a network,, 
##########################################the matrix is to be sorted  num_author and then by sum, #####################
{
e1=which(rownames(tempS) %in% c("POLR3A","BRF1"))
col=colorRampPalette(c("whitesmoke","red","red4"))(10) 
png("IMAGES/slice_of_network.png", width = 4800, height = 1600, res = 300)

par(mfrow = c(1, 3), pty = "s", mar = c(5, 5, 2, 2))  # mar is for better spacing
l=11;r=3 ##l=11 and r=3 from polr3a and brf1 interaction
temp2=tempC[c(-l:r)+e1[1],c(-l:r)+e1[2]]
# heatmap(temp2,Rowv = NA, Colv = NA,col=col,scale="none",cexRow = 0.5,cexCol = 0.5)
image(c(1:(l+r+1)),c(1:(l+r+1)),temp2,col=col,axes=F,xlab="",ylab="")
text(x = 1:(l+r+1), y = par("usr")[3] - 0.1, labels = rownames(temp2), srt = 90, adj = 1, xpd = TRUE)
text(y = 1:(l+r+1), x = par("usr")[1] - 0.1, labels = colnames(temp2), srt = 0, adj =1, xpd = TRUE)
box(lwd=3)
temp2=tempEC[c(-l:r)+e1[1],c(-l:r)+e1[2]]
# heatmap(temp2,Rowv = NA, Colv = NA,col=col,scale="none",cexRow = 0.5,cexCol = 0.5)
image(c(1:(l+r+1)),c(1:(l+r+1)),temp2,col=col,axes=F,xlab="",ylab="")
text(x = 1:(l+r+1), y = par("usr")[3] - 0.1, labels = rownames(temp2), srt = 90, adj = 1, xpd = TRUE)
text(y = 1:(l+r+1), x = par("usr")[1] - 0.1, labels = colnames(temp2), srt = 0, adj =1, xpd = TRUE)
box(lwd=3)
temp2=tempS[c(-l:r)+e1[1],c(-l:r)+e1[2]]
# heatmap(temp2,Rowv = NA, Colv = NA,col=col,scale="none",cexRow = 0.5,cexCol = 0.5)
image(c(1:(l+r+1)),c(1:(l+r+1)),temp2,col=col,axes=F,xlab="",ylab="")
text(x = 1:(l+r+1), y = par("usr")[3] - 0.1, labels = rownames(temp2), srt = 90, adj = 1, xpd = TRUE)
text(y = 1:(l+r+1), x = par("usr")[1] - 0.1, labels = colnames(temp2), srt = 0, adj =1, xpd = TRUE)
box(lwd=3)
dev.off()
}##luka##

###################!!!!!!!!Done till here



#Plotting a sample of zz plot
gimme_EC("POLR3H",improved_fit = T,neighbourhood = F,threshold = 0.1,legend = F,label = F)

}##luka##



########################################################################
############Understanding ZZ Plot#######################################
###################FIGURE1#### How Z-Z plot is calculated in simulated network? FLowchart!!!
#####and supp fig 1#############################################
#################SIMULATION######################################################
########################################################################
{
require(MASS)
require(entropy)
require(igraph)
require(dplyr)


###Important functions
convertC2EC=function(sim)
{A_s=sim;diag(A_s)=0;D_s=A_s*0;diag(D_s)=rowSums(A_s) ##Finding Admittance and degree
L_s=D_s-A_s  #Computing Normal Laplacian
z_s=ginv(L_s) # Finding pseudoinverse
#Making EC
FF_s=z_s
for (j in 1:nrow(FF_s)){for (k in 1:nrow(FF_s)) {FF_s[j, k] = z_s[j, j] + z_s[k, k] - z_s[j, k] - z_s[k, j]}}
EC_s=1/FF_s;colnames(EC_s)=colnames(sim);rownames(EC_s)=rownames(sim)
return(EC_s)
}

plot_matrix_by_dots= function(x,colors,cex=1,pch=19,lwd=10)
{
  ew<- summary(Matrix(x, sparse = TRUE)) %>%
  dplyr::arrange(x) %>%
  mutate(alpha= pmin(x/max(x),1),
         x_scaled=as.numeric(cut(x, breaks = length(colors), labels = FALSE)))

plot(c(1,ncol(x)),c(1,nrow(x)),col="transparent",axes=F,main="",xlab="",ylab="")
points(ew$i,ew$j,pch=pch, cex=cex,col=colors1[ew$x_scaled])
points(ew$j,ew$i,pch=pch, cex=cex,col=colors1[ew$x_scaled])
rect(1,1,nrow(x),ncol(x),lwd=lwd)
# box(lwd=lwd)
}##luka##

# quartz(width=8,height=12)
dev.off()
png("IMAGES/temp33.png", width = 4000, height = 4000, res = 300)

par(mfcol=c(3,4),pty="s")
################################STRATEGY1
##No experiments, but different power of interactions
{
sz=1000
num_exp=c(1:sz) #^2
# num_exp=c(1:250,rep(251,700),951:1000)
# num_exp=ew2
sim=matrix(data=0,nrow=sz,ncol=sz)
for(i in c(1:sz)){sim[i,]=sim[i,]+num_exp[i];sim[,i]=sim[,i]+num_exp[i]}
# if(T){image(sim,col=colors1,xaxt="none",yaxt="none",main="Conductance (C)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
plot_matrix_by_dots(sim,colors=colors1)
###########################Naming rows and columns
rownames(sim)=paste("P",c(1:sz),sep = "")
colnames(sim)=rownames(sim)

#######################
ECs1=convertC2EC(sim);ECs1[is.infinite(ECs1)]=0
# if(T){image(ECs1,col=colors1,xaxt="none",yaxt="none",main="Effective Conductance (EC)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
plot_matrix_by_dots(ECs1,colors = colors1)
tems1=ECs1;tems1[is.infinite(tems1)]=0;for(i in c(1:nrow(tems1))){tems1[i,]=scale(tems1[i,])}
# if(T){image(t(apply(exp(1)^tems1,2,rev)),col=colors,xaxt="none",yaxt="none",main="row z-normalized EC ")} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
# if(T){image(tems1,col=colors2,xaxt="none",yaxt="none",main="row z-normalized EC ");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
plot_matrix_by_dots(2^tems1,colors=colors1)
abline(v=700-10,lty=2);abline(v=700+10,lty=2)
abline(h=(700)-10,lty=2);abline(h=700+10,lty=2)
tems1_c=sim;for(i in c(1:nrow(tems1_c))){tems1_c[i,]=scale(tems1_c[i,])}
###########################Introducing bias
if(T){
  simb=sim
  simb["P700","P300"]=simb["P700","P300"]+200
  simb["P300","P700"]=simb["P300","P700"]+200
  # if(T){image(simb,col=colors1,xaxt="none",yaxt="none",main="Conductance (C)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
plot_matrix_by_dots(simb,colors=colors1)
  
  
  bECs1=convertC2EC(simb);bECs1[is.infinite(bECs1)]=0
  # if(T){image(bECs1,col=colors1,xaxt="none",yaxt="none",main="Effective Conductance (EC)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
plot_matrix_by_dots(bECs1,colors=colors1)
  btems1=bECs1;btems1[is.infinite(btems1)]=0;for(i in c(1:nrow(btems1))){btems1[i,]=scale(btems1[i,])}
  # if(T){image(btems1,col=colors2,xaxt="none",yaxt="none",main="row z-normalized EC ");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
plot_matrix_by_dots(2^btems1,colors=colors1)
  abline(v=700-10,lty=2);abline(v=700+10,lty=2)
  abline(h=700-10,lty=2);abline(h=(700)+10,lty=2)
  btems1_c=simb;btems1_c[is.infinite(btems1_c)]=0;for(i in c(1:nrow(btems1_c))){btems1_c[i,]=scale(btems1_c[i,])}
}
}##luka##

png("IMAGES/temp3.png", width = 4000, height = 4000, res = 300)
par(mfcol=c(3,4),pty="s")
################################STRATEGY2
##differnet number of experiments done in proteins
{
sz=5000
prob= 50/sz # Every protein has 50 interacting protein out of 1000 protein
num_exp=c(1:sz)
sim=matrix(data=0,nrow=sz,ncol=sz)
for(i in c(1:sz))
{
  print(i)
  for(rep in c(1:num_exp[i]))
  {
    indx=sample(c(1:sz),prob*sz)
    sim[i,indx]=sim[i,indx]+1
    sim[indx,i]=sim[indx,i]+1
  }
}
saveRDS(sim,"sim2_5000_50.Rdata")
# if(T){image(sim,col=colors1,xaxt="none",yaxt="none",main="Conductance (C)");box(lwd=3)}
sim=readRDS("sim2_5000_50.Rdata")

###########################Naming rows and columns
rownames(sim)=paste("P",c(1:sz),sep = "")
colnames(sim)=rownames(sim)
#######################
ECs2=convertC2EC(sim);ECs2[is.infinite(ECs2)]=0
# if(T){image(ECs2,col=colors1,xaxt="none",yaxt="none",main="Effective Conductance (EC)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
# plot_matrix_by_dots(ECs2,colors = colors1)
tems2=ECs2;tems2[is.infinite(tems2)]=0;for(i in c(1:nrow(tems2))){tems2[i,]=scale(tems2[i,])}
# if(T){image(tems2,col=colors2,xaxt="none",yaxt="none",main="row z-normalized EC");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
# plot_matrix_by_dots(2^tems2,colors=colors1)
# abline(v=700-10,lty=2);abline(v=700+10,lty=2)
# abline(h=(700)-10,lty=2);abline(h=700+10,lty=2)
tems2_c=sim;for(i in c(1:nrow(tems2_c))){tems2_c[i,]=scale(tems2_c[i,])}
###########################Introducing bias
if(F){
  simb=sim
  simb["P700","P300"]=simb["P700","P300"]+100
  simb["P300","P700"]=simb["P300","P700"]+100
  simb["P5","P300"]=simb["P5","P300"]+100
  simb["P300","P5"]=simb["P300","P5"]+100
  # if(T){image(simb,col=colors1,xaxt="none",yaxt="none",main="Conductance (C)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
 
  # png("IMAGES/temp.png",width = 1000,height = 3000,res=100)
  # par(mfcol=c(3,1)) 
  plot_matrix_by_dots(simb,colors=colors1)
  
  bECs2=convertC2EC(simb);bECs2[is.infinite(bECs2)]=0
  # if(T){image(bECs2,col=colors1,xaxt="none",yaxt="none",main="Effective Conductance (EC)");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
  plot_matrix_by_dots(bECs2,colors=colors1)
  btems2=bECs2;btems2[is.infinite(btems2)]=0;for(i in c(1:nrow(btems2))){btems2[i,]=scale(btems2[i,])}
  # if(T){image(btems2,col=colors2,xaxt="none",yaxt="none",main="row z-normalized EC ");box(lwd=3)} #have only 100 by 100 matrix to show granularity in simulation-2 and smoothness in simulation-1
  plot_matrix_by_dots(2^btems2,colors=colors1)
  abline(v=700-10,lty=2);abline(v=700+10,lty=2)
  abline(h=(700)-10,lty=2);abline(h=700+10,lty=2)
  btems2_c=simb;btems2_c[is.infinite(btems2_c)]=0;for(i in c(1:nrow(btems2_c))){btems2_c[i,]=scale(btems2_c[i,])}

}

###########################Introducing  bias (via 3 intermediates)
# for(numi in c(2,4,6,8,10,12))
{ 
  numi=3 #number of intermediates
  print(numi)
  p1="P140" #protein A
  # p2="P2"
  which(sim[p1,]==0)
  p3="P430" #Manually chosing C so that it doesnt has any interaction with P1
  which(sim[p1,]==0)
  # set.seed(453)
  set.seed(2323544)
  p2_list=sample(setdiff(paste0("P",c(1:sz)),c(p1,p3)),numi) #Choosing B1 B2 B3 proteins.
  p2_list
  
  # for(incrs in c(1,2,3,6,8,10)*1000)
      {
        incrs= 5000
    simb=sim
    for(p2 in p2_list)
    { 
      simb[p1,p2]=simb[p1,p2]+incrs
      simb[p2,p1]=simb[p2,p1]+incrs
      simb[p2,p3]=simb[p2,p3]+incrs
      simb[p3,p2]=simb[p3,p2]+incrs
    }
    
    ##colorbar
    { 
      png("IMAGES/color_bar.png",width=1000,height=1000,res = 300)
      colors1=colorRampPalette(c(rep("white",1),rep("orange",1),rep("red3",1),rep("red4",100)))(5000)
      # colors=colorRampPalette(c("white","red"))(100)
      # plot(c(1:5000),rep(1,5000),col=colors1,pch=19);abline(v=1)
      plot(c(1,310),c(0,5),col="transparent",axes=F,xlab="",ylab="",main="Color bar")
      sapply(c(1:300),function(i) points(c(i,i),c(0.5,2),type="l",lwd=2,col=colors1[i]))
      # rect(1,0.5,par("usr")[2],2,lwd=2) # making a rectange 
      rect(1,0.5,300,2,lwd=2) # making a rectange 
      
      scor=c(min(simb),max(simb)) #actual scores at two ends of colorbar
      # scor=log(1+scor)
      sort(unlist(simb),decreasing = T)[1:20]
      # scor=c(0,5000)
      coor=c(1,5000) #x coordinate positions of two end of color bar
      scors=c(50,100,150,290,500,2500,5000)
      sapply(c(1:length(scors)),function(x) 
      {
        i=(scors[x])
        i2= coor[1] + (i- scor[1]) * (diff(coor)/diff(scor)) # derived from y= y1 + (x-x1) * ((y2-y1)/(x2-x1))
        points(c(i2,i2),c(0.4,0.5),type="l",lwd=2)
        text(i2,0.1,labels=ifelse(scors[x]==290,5000,scors[x])) #writing the labels, but also writing 5000 instead of 290
      }
      )
      dev.off()
      }##Luka##
    png("IMAGES/sim2_biased.png", width = 1000, height = 1000, res = 300)
    plot_matrix_by_dots(simb,colors=colors1)
    dev.off()
    png("IMAGES/sim2_unbiased.png", width = 4000, height = 4000, res = 300)
    {
      ew<- summary(Matrix(sim, sparse = TRUE)) %>%
        dplyr::arrange(x)
      ew$x_scaled=as.numeric(cut(c(ew$x,max(simb)), breaks = length(colors1), labels = FALSE))[1:nrow(ew)]
      plot(c(1,ncol(sim)),c(1,nrow(sim)),col="transparent",axes=F,main="",xlab="",ylab="")
      points(ew$i,ew$j,pch=19, cex=1,col=colors1[ew$x_scaled])
      points(ew$j,ew$i,pch=19, cex=1,col=colors1[ew$x_scaled])
      rect(1,1,nrow(sim),ncol(sim),lwd=10)
    }##luka## plotting the original sim2 matrix using same scalebar as simb2_biased plot.
    dev.off()
    
    bECs2=convertC2EC(simb);bECs2[is.infinite(bECs2)]=0
    btems2=bECs2;btems2[is.infinite(btems2)]=0;for(i in c(1:nrow(btems2))){btems2[i,]=scale(btems2[i,])}
    saveRDS(btems2,"btems2_2000incrs.Rdata")
    btems2=readRDS("btems2_5000incrs.Rdata")
    #Plotting the reciprocal plots of A in simulation2 and biased_simulation2
    {
    quartz(width=4, height = 4, bg="white")
    par(pty="s")
    wew1=gimme_EC(p1,del = tems2,legend = F,threshold = 0.0001,improved_fit = T,main = j,label=F,plot = F)
    plot(wew1$self,wew1$other,cex=1,pch=19,col="skyblue",xaxt="n",yaxt="n",
         main=paste0("Reciprocal plot for protein ", p1), xlab= "Z score",ylab="Reciprocal Z score")
    axis(1, at=c(-2,-1,0,1,2), las=1)
    axis(2, at=c(-1.85,-1.9,-1.95), las=1)
    points(wew1[p2_list,c("self","other")],bg="green",pch=21,cex=1)
    points(wew1[p3,c("self","other")],bg="red",pch=21,cex=1)
    quartz.save("IMAGES/reciprocal_plot_of_A_in_unbiased_sim2.png",dpi=800)
    dev.off()
    quartz(width=4, height = 4, bg="white")
    par(pty="s")
    wew1=gimme_EC(p1,del = btems2,legend = F,threshold = 0.0001,improved_fit = T,main = j,label=F,plot = F)
    plot(wew1$self,wew1$other,cex=1,pch=19,col="skyblue",xaxt="n",yaxt="n",
         main=paste0("Reciprocal plot for protein ", p1), xlab= "Z score",ylab="Reciprocal Z score")
    axis(1, at=c(-2,-1,0,1,2), las=1)
    axis(2, at=c(-1.45,-1.5,-1.55,-1.6), las=1)
    points(wew1[p2_list,c("self","other")],bg="green",pch=21,cex=1)
    points(wew1[p3,c("self","other")],bg="red",pch=21,cex=1)
  quartz.save("IMAGES/reciprocal_plot_of_A_in_biased_sim2.png",dpi=800)
  dev.off()
  
  }##luka##
  
    }
  
}

}##luka#


###############################Plotting colorbars for all the strategies
{
plot(c(1:1000),rep(1,1000),col=colors1,pch=19);abline(v=1)
plot(c(1,1000),c(0,5),col="transparent",axes=F,xlab="",ylab="",main="Effective Conductance")
sapply(c(1:1000),function(i) points(c(i,i),c(0.5,2),type="l",lwd=0.4,col=colors1[i]))
rect(1,0.5,1000,2,lwd=2)

temp=sim
scor=c(min(temp$x),max(temp$x)) #actual scores at two ends of colorbar
coor=c(1,1000) #x coordinate positions of two end of color bar

gimme_coordinate=function(score,ew,coor) #derived from (y2-y2)/(x2-x1)=(y-y1)/(x-x1) where x is coordinates at color bar, and y is log transformed score,i.e. log10(100+score)
{enquired_coor=coor[1]+diff(coor) * (log10(100+score)-ew[1]) / (ew[2]-ew[1])
return(enquired_coor)
}
gimme_coordinate(1000,ew,coor)

sapply(c(10,100,1000),function(i) 
{
  i2=gimme_coordinate(i,ew,coor)
  points(c(i2,i2),c(0.4,0.5),type="l",lwd=2)
  text(i2,0.3,labels=i)
}
)
}##luka##

#################################Plotting
### I TRIED TO MAKE IT IN PHASE BY SYNCING THE TOTAL INTERACTION ,IT DINTO WIRK, EVEN NUM_AUTHOR DOESNT WORK,THE BEST IS TO SYNC THE INDEX, I.E. LINEARLY SYNCING 1000(OF SIMULAITON) WITH 20K(OF BIOGRID)
###, TRY TO SIMULATE(STRATEGY1) WITH  # num_exp=c(1:250,rep(251,750)) AND THEN THE ZZ PLOT DOESNT EVEN START FROM 2*PI. I THINK THERE IS SOMETHIGN VERY FUNCDAMENTAL GOING ON AND WHICH I DONT UNDERSTAND.
####ALSO,COMPARE WITH THE SIMULATION(STRATEGY1) WITH # num_exp=c(1:250,rep(251,700),951:1000)
{
tem=readRDS("tem.Rdata")
bias=readRDS("bias.Rdata") %>% mutate(sum=A+B) %>% arrange(sum)
bias["Pwhat"]=paste0("P",as.integer(MinMaxNorm2(bias$sum,1,4000)))
bias2=bias %>% group_by(Pwhat) %>% summarize(Var1=Var1[1],sum=sum[1]) %>% arrange(sum)


#plotting video of z z plot for siimulation 1 ,2 and biogrid data
par(mfrow=c(3,1))
   for(i in seq(1,100,by=2))
     # for(i in seq(650,800,by=1))
{
  print(i)
     P=paste0("P",i)
     prot=bias$Var1[i*20]

  wew1=data.frame(self=tems1[P,],other=tems1[,P]) %>% filter(rownames(.)!=P)
  wew2=data.frame(self=tems2[P,],other=tems2[,P]) %>% filter(rownames(.)!=P)
  wew3=data.frame(self=tem[prot,],other=tem[,prot]) %>% filter(rownames(.)!=prot)
  
  plot(wew1$self,wew1$other,col="blue")
   plot(wew2$self,wew2$other,col="blue")
  plot(wew3$self,wew3$other,col="blue")
  
  Sys.sleep(0.1)
}
dev.off()


#Crating  animation of z z plot,of simulation 1 and 2 only
for(i in seq(1,100,by=2))
  # for(i in seq(650,800,by=1))
{
  # break
  print(i)
  P=paste0("P",i)
  prot=bias$Var1[i*20]
  col1=colorRampPalette(c("blue","whitesmoke","red"))(1000)
  
  wew1=data.frame(self=tems1[P,],other=tems1[,P]) %>% filter(rownames(.)!=P)
  wew1_c=data.frame(self=tems1_c[P,],other=tems1_c[,P]) %>% filter(rownames(.)!=P)
  
  wew2=data.frame(self=tems2[P,],other=tems2[,P]) %>% filter(rownames(.)!=P)
  wew2_c=data.frame(self=tems2_c[P,],other=tems2_c[,P]) %>% filter(rownames(.)!=P)
  
  jpeg(file=sprintf("delete_fast_%06d.jpg", i),width=1500,height=1000,quality=100)
  par(mfrow=c(3,4),cex.lab=1.5,cex=1.2,cex.axis=1.2,
      mar = c(3, 3, 3, 3),    # margins inside each plot (bottom, left, top, right)
      oma = c(0, 0, 0, 0),    # outer margins (around all plots)
      mgp = c(1.5, 0.5, 0),   # axis title/label spacing
      xpd = NA)               # allow plotting outside plot region if needed
  
  #plotting popularity scale
  plot(1,1,col='transparent',axes=F,xlab="",ylab="")
  plot(c(1,1000),c(2,5),col="transparent",axes=F,xlab="",ylab="")
  sapply(c(1:1000),function(i) points(c(i,i),c(1,2),type="l",lwd=2.5,col=col1[i]))
  rect(1,1,1000,2,lwd=2)
  points(c(i,i),c(1,2),type="l",lwd=3,col="black")
  text(i,0.8,P)
  text(500,2.2,"Popularity",cex=2)
  plot(1,1,col='transparent',axes=F,xlab="",ylab="")
  plot(1,1,col='transparent',axes=F,xlab="",ylab="")
  
  #plotting all z scores
  pad=0.2 # how much padding to have in that dotted rectangle
  plot(wew1_c$self,wew1_c$other,col=col1,xlim=c(-3,3),ylim=c(-3,3),xlab="Zscore",ylab="Reciprocal Z score");box(lwd=4)
  { 
  bottom=min(wew1_c$other)-pad;top=max(wew1_c$other)+pad #coordinates of dotted rectangle
  left=min(wew1_c$self)-pad;right=max(wew1_c$self)+pad   #coordinates of dotted rectangle
  rect(left,bottom,right,top,lty=2,lwd=2)                # creating dotted rectangle
  points(c(right,5.7),c(top,3.15),type="l",lwd=2,lty=2)  #line connecting dotted rectangle to another plot
  points(c(right,5.7),c(bottom,-3.15),type="l",lwd=2,lty=2) #line connecting dotted rectangle  to another plot
  }
  plot(wew1_c$self,wew1_c$other,col=col1,xlab="",ylab="",bty="n");  box(lwd=2,lty=2)
# plot(1,1,pch=21,col="red")
  
  plot(wew2_c$self,wew2_c$other,col=col1,xlim=c(-3,3),ylim=c(-3,3),xlab="Zscore",ylab="Reciprocal Z score");box(lwd=4)
  { 
    bottom=min(wew2_c$other)-pad;top=max(wew2_c$other)+pad
    left=min(wew2_c$self)-pad;right=max(wew2_c$self)+pad
    rect(left,bottom,right,top,lty=2,lwd=2)
    points(c(right,5.7),c(top,3.15),type="l",lwd=2,lty=2)
    points(c(right,5.7),c(bottom,-3.15),type="l",lwd=2,lty=2)
  }
  plot(wew2_c$self,wew2_c$other,col=col1,xlab="",ylab="",bty="n");  box(lwd=2,lty=2)
       
  plot(wew1$self,wew1$other,col=col1,xlim=c(-3,3),ylim=c(-3,3),xlab="Zscore",ylab="Reciprocal Z score");box(lwd=4)
  { 
    bottom=min(wew1$other)-pad;top=max(wew1$other)+pad
    left=min(wew1$self)-pad;right=max(wew1$self)+pad
    rect(left,bottom,right,top,lty=2,lwd=2)
    points(c(right,5.7),c(top,3.15),type="l",lwd=2,lty=2)
    points(c(right,5.7),c(bottom,-3.15),type="l",lwd=2,lty=2)
  }
  plot(wew1$self,wew1$other,col=col1,xlab="",ylab="",bty="n");  box(lwd=2,lty=2)
  
  plot(wew2$self,wew2$other,col=col1,xlim=c(-3,3),ylim=c(-3,3),xlab="Zscore",ylab="Reciprocal Z score");box(lwd=4)
  { 
    bottom=min(wew2$other)-pad;top=max(wew2$other)+pad
    left=min(wew2$self)-pad;right=max(wew2$self)+pad
    rect(left,bottom,right,top,lty=2,lwd=2)
    points(c(right,5.7),c(top,3.15),type="l",lwd=2,lty=2)
    points(c(right,5.7),c(bottom,-3.15),type="l",lwd=2,lty=2)
  }
  plot(wew2$self,wew2$other,col=col1,xlab="",ylab="",bty="n");  box(lwd=2,lty=2)
dev.off()
  #Sys.sleep(0.1)
}
library(magick)
image_files <- list.files("./", pattern = "delete_fast_", full.names = TRUE)
image_list <- lapply(image_files, image_read)
img_all <- image_join(image_list)
animation <- image_animate(img_all, fps = 5)
image_write(animation, path = "ZZ_animation.gif",)
system("rm delete_fast_*")



#Plotting static image of z z  for comparison
{
par(mfcol=c(2,2))
col1=colorRampPalette(c("blue","whitesmoke","red"))(1000)
P="P656"
wew1=data.frame(self=tems1[P,],other=tems1[,P]) %>% filter(rownames(.)!=P) ;  plot(wew1$self,wew1$other,col=col1,xlab="Z score",ylab="Reciprocal Z score",main="Z-Z plot of P65")
points(tems1[P,"P900"],tems1["P900",P],cex=1);text(labels="P900",tems1[P,"900"],tems1["P900",P]+0.1*diff(range(wew1$other)))

wew1=data.frame(self=btems1[P,],other=btems1[,P]) %>% filter(rownames(.)!=P) ;  plot(wew1$self,wew1$other,col=col1,xlab="Z score",ylab="Reciprocal Z score",main=c("Z-Z plot of P656 ","+","an unbiased interaction","with P900 introduced"))
points(btems1[P,"P900"],btems1["P900",P],cex=1);text(labels="P900",btems1[P,"P900"],btems1["P900",P]+0.1*diff(range(wew1$other)))

wew1=data.frame(self=tems2[P,],other=tems2[,P]) %>% filter(rownames(.)!=P) ;  plot(wew1$self,wew1$other,col=col1,xlab="Z score",ylab="Reciprocal Z score",main="Z-Z plot of P656")
points(tems2[P,"P900"],tems2["P900",P],cex=1);text(labels="P900",tems2[P,"P900"],tems2["P900",P]+0.1*diff(range(wew1$other)))

wew1=data.frame(self=btems2[P,],other=btems2[,P]) %>% filter(rownames(.)!=P) ;  plot(wew1$self,wew1$other,col=col1,xlab="Z score",ylab="Reciprocal Z score",main=c("Z-Z plot of P656 ","+","an unbiased interaction","with P900 introduced"))
points(btems2[P,"P900"],btems2["P900",P],cex=1);text(labels="P900",btems2[P,"P900"],btems2["P900",P]+0.1*diff(range(wew1$other)))

dev.off()
i=656
wew1=data.frame(self=tem[bias$Var1[i*20],],other=tem[,bias$Var1[i*20]]) %>% filter(rownames(.)!=bias$Var1[i*20]) ;  plot(wew1$self,wew1$other,col="cyan2",xlab="Z score",ylab="Reciprocal Z score",main="Z-Z Plot of PLEKHM1 (BIOGRID data)")
ew=gimme_EC(bias$Var1[i*20],plot=F) %>% filter(distances>0.01)
points(ew$self,ew$other,col="orange",cex=2,pch=19)
ew=ew %>% filter(distances>0.03)
text(labels=rownames(ew),ew$self,ew$other,cex=0.6)

}##luka##

}##luka##





}##luka##



########################################################################
########################################################################
###################FIGURE supp######################################################
##################ZZ plot for simulation of network with different patterns of distribution of popularity####################
########################################################################
{
  gimme_tem=function(num_exp)
  { 
    P="P250"
    sz=1000
    sim=matrix(data=0,nrow=sz,ncol=sz)
    for(i in c(1:sz)){sim[i,]=sim[i,]+num_exp[i];sim[,i]=sim[,i]+num_exp[i]}
    rownames(sim)=paste("P",c(1:sz),sep = "")
    colnames(sim)=rownames(sim)
    ECs1=convertC2EC(sim)
    diag(ECs1)=0
    tems1=ECs1;for(i in c(1:nrow(tems1))){tems1[i,]=scale(tems1[i,])}
    return(tems1)
  }
  e1=c(1:sz);te1=gimme_tem(e1)
  e2=sort(sample(c(1:10000),1000,replace=T));te2=gimme_tem(e2)
  e3=c(1:250,rep(251,250),251:750);te3=gimme_tem(e3)
  e4=c(1:250,501:1250);te4=gimme_tem(e4)
  e5=c(1:500,rep(501,500));te5=gimme_tem(e5)
  e6=1+as.integer(c( 2^((1:980)/100)/10,5^c((51:70)/15))); te6=gimme_tem(e6)
  
  E= list(list(e1,te1),
          list(e2,te2),
          list(e3,te3),
          list(e4,te4),
          list(e5,te5),
          list(e6,te6))
  
for(i in seq(1,1000,by=20))
{ 
    print(i)
    P=paste0("P",i)
    col1=colorRampPalette(c("blue","whitesmoke","red"))(1000)
    jpeg(file=sprintf("delete_fast_%06d.jpg", i),width=700,height=1500,quality=100)
    par(mfrow=c(length(E),3),cex.lab=1.5,cex=1.2,cex.axis=1.2,
        mar = c(1, 2, 1, 1),    # margins inside each plot (bottom, left, top, right)
        oma = c(0, 0, 0, 0),    # outer margins (around all plots)
        mgp = c(1.5, 0.5, 0))   # axis title/label spacing
    
for(j in c(1:length(E))) 
{ 
   e=E[[j]][[1]]
   te=E[[j]][[2]]
   plot(e,col=col1,xaxt="none",yaxt="none",ylab="");abline(v=i,lty=2)
   plot(cumsum(e),e,col=col1,xaxt="none",yaxt="none",ylab="");abline(v=cumsum(e)[i],lty=2)
   wew1=data.frame(self=te[P,],other=te[,P]) %>% filter(rownames(.)!=P) 
 plot(wew1$self,wew1$other,col=col1,xlab="",ylab="",xaxt="none",yaxt="none",cex=0.4)  
 }
    dev.off()
    #Sys.sleep(0.1)
  }
  library(magick)
  image_files <- list.files("./", pattern = "delete_fast_", full.names = TRUE)
  image_list <- lapply(image_files, image_read)
  img_all <- image_join(image_list)
  animation <- image_animate(img_all, fps = 5)
  image_write(animation, path = "ZZ_distr__animation.gif")
  system("rm delete_fast_*")
}##luka##








