


print("The functions are :")


#######################################
#######################################
####################Rank Order Normalization (RON)
#######################################
print("ron")
ron=function(df)
{
require(dplyr)
sorted= as.data.frame(apply(df,2,sort))
avg=apply(sorted,1,mean)
for(i in (1:ncol(df)))
{k=as.data.frame(cbind(1:nrow(df),df[,i])) %>% arrange(.,across(2))
l=as.data.frame(cbind(k[,1],avg))%>% arrange(.,across(1))
df[,i]=l[,2]}
return(df)
}


#######################################
#######################################
####################MINMAXNORMALIZATION
#######################################
print("MinMaxNorm")
MinMaxNorm=function(x)
{
  if(var(x)==0) return(x*0)
  mn=min(x)
  mx=max(x)
  x2=as.double((x-mn)/(mx-mn))
  return(x2)
} 

#print("use MinMaxNorm2([numberlist],min_limit,max_limit) to set min and max limit")
print("MinMaxNorm2")
MinMaxNorm2=function(x,min,max)
{
  MIN=min(x)
  MAX=max(x)
  y=as.double((x-MIN)*(max-min)/(MAX-MIN)+min)
  return(y)
} 






#######################################
#######################################
#################### REMOVE INTERSECTING GENES
#######################################
print("removeintersecting")
removeintersecting=function(a,b,opposite=F,extra_flags="")
{
  flags=ifelse(opposite==F,"-v","-u")
  print("Writing files to perform bedtools intersect.......")
  write.table(a ,"a.bed",sep="\t",row.names = F,col.names = F,quote=F,append=F)
  write.table(b ,"b.bed",sep="\t",row.names = F,col.names = F,quote=F,append=F)
  flags=paste(flags,extra_flags)
  code=paste("bedtools intersect",flags,"-a a.bed -b b.bed > c.bed")
  print(paste("Running -->",code))
  system(code)
  print("Now reading the intersected file and shortly returning it.....")
  c=fread("c.bed",sep="\t")
  colnames(c)=colnames(a)
  system("rm a.bed b.bed c.bed")
  return(c)
}


#######################################
#######################################
#################### What genes intersect
#######################################
print("whatintersect")
whatintersect=function(a,b,cola=NULL,colb=NULL,f=NULL)
{
  fwrite(a,"a.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)
  fwrite(b,"b.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)
  if(is.null(f)) {system("bedtools intersect -wa -wb -a a.bed -b b.bed > c.bed ")} else
  {system(paste0("bedtools intersect -wa -wb -f ",f," -a a.bed -b b.bed > c.bed "))}
  c=fread("c.bed",sep="\t")
  if(is.null(cola)){cola=colnames(a)}
  if(is.null(colb)){colb=paste0(colnames(b),"_")}
  colnames(c)=c(cola,colb)
  system("rm a.bed b.bed c.bed")
  return(c)
}





#######################################
#######################################
####################MISCELLENEOUS for triangleplotting
#######################################
  ##making colorgradient 
  if(F){
    fun_color_range <- colorRampPalette(c("#1b98e0", "red"))   # Apply colorRampPalette
    my_colors <- fun_color_range(100)                          # Extract 100 color codes
    plot(1:100,pch = 16,col = my_colors)
  }##luka##


  ##function to make an equilateral triangle
print("make triangle")
  maketriangle=function(titl="TITLE",new=T,cex.main=1,cex.vertex = 1)
  {
    if(new==T)
    {plot(c(0,1,0.5),c(0,0,sqrt(3)/2),pch=1,xlim=c(-0.1,1.1),ylim=c(-0.1,1.1),main=titl,cex.main=cex.main,axes=F,xlab="",ylab="",col="transparent")}
    polygon(c(0,1,0.5),c(0,0,sqrt(3)/2),pch=16)#######,xlim=c(-0.2,1),ylim=c(0,0.9))
    if(cex.vertex!=0)
   {text(-0.05,0.1,"Pol2",cex =cex.vertex )
    text(1.05,0.1,"Pol3",cex = cex.vertex)
    text(0.5,0.95,"Pol1",cex = cex.vertex)}
  }


  ###function that takes 6 digit color and transparency in (0,255), then returns 8dgt color
print("gimme_8dgt")
  gimme_8dgt=function(kal,tran)
  {
    #kal="#123456"
    #tran=tran
    k=col2rgb(kal)
    return(rgb(k[1],k[2],k[3],tran,maxColorValue = 255))
  }


##getCOM
print("getCOM")
getCOM=function(data,x=c(0,1,0.5) ,y=c(0,0,sqrt(3)/2))
{return(data.frame(comx=rowSums(as.matrix(data)%*%diag(x)),comy=rowSums(as.matrix(data)%*%diag(y)))/rowSums(data))}

##entropy
print("ent     (for finding entropy)")
ent=function(data)
{
  data=data/rowSums(data)
  eata=log2(data)
  fata=eata %>% dplyr::mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))
  return(-rowSums(data*fata))
  }


##SCORING function
print("SCOORING     (for finding poisson enrichment SCOREs across the bins.just give a bedgraph file containing 3 coordinate columns and raw signal in 4th column )")
SCOORING=function(bed5,binsize=50)
{
  require(data.table)
  require(RcppRoll)
  require(dplyr)
  bed5=as.data.frame(bed5)
  FLAG=0
  colnames(bed5)=c("V1","V2","V3","x")
  bed5$x=round(bed5$x)
  chr=unique(bed5$V1)
  for(i in chr)
  {
    if(i=="" | i=="chrM") next #removing chrM coz its really short and scoring chrM is not possible coz of the extended region to find lambda, exceeds chromose M size
    print(i)
    temp= bed5 %>% filter(V1==i)
    temp= temp %>% arrange(V1,V2)
    #Finding various lambdas
    for(rng in c(5000,10000))
    {print(rng)
      raange=rng/binsize #number of bins in either side to take
      wnd=raange*2+1
      temp=temp %>% mutate(avgg=c(rep(0,(wnd-1)/2),roll_mean(x,n=wnd, align = "center"),rep(0,(wnd-1)/2))) %>%
        mutate(lambda=avgg*wnd/(wnd-1)-(x/(wnd-1)))
      colnames(temp)[which(colnames(temp)=="lambda")]=paste0("lambda_",rng)
    }
    
    ##finding Pvalue and FDR and scores
    gwm=mean(bed5$x)
    temp= temp %>% mutate(lambda=pmax(gwm,lambda_5000,lambda_10000))
    print("max lambda found")
    #temp[temp$xpct<0]=0
    temp= temp %>% mutate(pval=ppois(x,lambda,lower.tail = FALSE,log=F))
    print("pvalue found")
    temp=temp %>% mutate(FDR=p.adjust(pval,method="BH"))
    print("FDR found")
    temp=temp %>% mutate(score=-log10(FDR))
    temp$score[c(1:raange,(nrow(temp)+1-raange):nrow(temp))]=0  #adjusting the scores of the ends
    
    #Removing infinites scores (Assigning them score of 2nd_max_score+50. Second max score coz max score is infinite itself)
    howmanyinfinite=sum(is.infinite(temp$score))
    if(howmanyinfinite>0)
    {
      temp$score[which(is.infinite(temp$score))]=sort(temp$score[!is.infinite(temp$score)],decreasing=T)[1]#for infinite we assign highest non-infinite score  
      temp$score=round(temp$score,digits=2)                  #rounding off score values to 2 digits
    }
    #write.table(temp %>% select(V1,V2,V3,score),"pol3.bedgraph",row.names = F,col.names = F,quote=F,sep="\t",append=F)
    if(FLAG==0){BUDU=temp} else
    {
      BUDU=rbind(BUDU,temp)
    }
    FLAG=FLAG+1
  }##luka##
  
  ##FINDING SCORE2 ,i.e. dividing by log2(freq of that signal+1) 
  # {scr= as.data.frame(BUDU) %>% mutate(score2=as.integer(score)) #assigning closest integer
  #   temp=as.data.frame(table(scr$score2)) #finding frequency
  #   temp= temp %>% transmute(score2=as.integer(as.character(Var1)),freq=Freq)
  #   scr=left_join(scr,temp)
  #   scr = scr %>% mutate(score2=score/log2(freq+1))
  #   scr$score2=round(scr$score2,digits=2)  
  # }##luka##
  
  # Maintaining motonicity of SCOREs when dividing by log2(freq+1)
  # temp=scr %>% arrange(score)
  # FLAG=1
  # while(FLAG!=-9999)
  # {
  #   temp=temp %>% mutate(diff=temp$score2-c(0,temp$score2[1:(nrow(scr)-1)]))
  #   breaks=which(temp$diff<0)
  #   print(paste0("number of monotonic breaks in ",FLAG,"th cycle = ",length(breaks)))
  #   FLAG=FLAG+1
  #   if(length(breaks)>0)
  #   {temp$score2[breaks]=temp$score2[breaks-1]} else {FLAG=-9999}
  # }##luka##
  
  # scr=temp %>% dplyr::select(-c(diff)) %>% arrange(V1,V2)
  # return(scr)
  return(as.data.frame(BUDU))
}##luka##


##give me neighbour from interaction matrix
print("gimme_neighbour(matrix,gene)")
gimme_neighbour=function(matrix,gene)
{
   # gene="BHD2"
   if(!gene %in% rownames(matrix)) {print(paste0(gene, " is not present in the matrix"))}
  tot=c() #making list of neighbours starting with first neighbour
  repeat
  {
    if(length(gene)==1){rsum=unlist(matrix[,gene])} else {rsum=rowSums(matrix[,gene])}
    gene=setdiff(names(which(rsum>=1)),unlist(tot))
    if(length(gene)==0) break
    tot=c(tot,list(gene))
  }
  tot=c(tot,list(c(setdiff(rownames(matrix),unlist(tot)),"Disconnected_Proteins")))
   #putting remaining disconnected genes
  return(tot)
}



##bedtools
print('bedtools(a,b,commands="merge -flag1 -flag2 ..... -a a.bed -b bed ",display.logs=T)')
bedtools=function(a=NULL,b=NULL,commands=NULL,display.logs=T)
{
  if(is.null(commands)) {print("No commands and flags were given"); break}
  if(display.logs==T){print("Writing files in progress ....")}
  if(!is.null(a)){fwrite(a,"a.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)}
  if(!is.null(b)){fwrite(b,"b.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)}
  code=paste0("bedtools ",commands," > c.bed ")
  if(display.logs==T){print("Performing");print(code)}
  system(code)
  if(display.logs==T){print("Now reading and returning the output file ...")}
  c=fread("c.bed",sep="\t")
  system("rm a.bed b.bed c.bed")
  return(c)
}


######
##ENRIZ
#####
print('enriz(tf,threshold=5000,no_permutation=5000)')
enriz=function(tf,threshold=5000,no_permutation=50000)
{
  
  colnames(tf)=c("ID","Ant.class","Antigen","Cell_class","Cell","No.ofPeaks","overlaps_A","overlaps_B","logpval","logqval","FE")
  #removing infinite FE (asuming geometric progression)
  if(is.infinite(max(tf$FE))==T){
    top3=sort(unique(tf$FE),decreasing=T)[1:3]
    tf$FE[is.infinite(tf$FE)]= (top3[2]^2)/top3[3]}
  #removing 0 FE
  if(min(tf$FE)==0){
    bottom3=sort(unique(tf$FE),decreasing=F)[1:3]
    tf$FE[tf$FE==0]= (bottom3[2]^2)/bottom3[3]}
  #removing experiment with less than 5000 peaks
  plot(density(tf$No.ofPeaks));abline(v=threshold)
  tf2=tf %>% dplyr::filter(No.ofPeaks>=threshold)
  tf2=tf2 %>% mutate(overlaps_A=as.double(sapply(strsplit(overlaps_A,"/"),"[[",1)),
                     overlaps_B=as.double(sapply(strsplit(overlaps_B,"/"),"[[",1))) ##grabbing 2 from say 2/56
  
  ##METHOD A
  tfA=tf2 %>% mutate(SCORE=(-logqval)*log10(FE)) %>% group_by(Antigen,Cell_class,Cell) %>% summarize(SCORE=mean(SCORE)) %>% ungroup()
  tfAA=tfA %>% group_by(Antigen) %>% summarize(SCORE=median(SCORE))
  plot(sort(unlist(scale(tfA$SCORE))))
  

  ##METHOD B
  tfB=tf2 %>% group_by(Antigen) %>% summarize(MO=mean(overlaps_A),samples=length(Cell))
  print("Now, MAking Permutation Table ")
  hu <- matrix(rep(0,no_permutation*length(unique(tfB$samples))), nrow = no_permutation) #columns refers to various sample sizes, rows are for vlaues generated during permutation
  colnames(hu)=unique(tfB$samples)
  for(i in c(1:ncol(hu)))
  {print(paste0(i," / ",ncol(hu)));hu[,i]=sapply(1:no_permutation, function(x) mean(sample(tf2$overlaps_A,as.integer(colnames(hu)[i]))))}
  hu2=round(hu,1)
  
  tfBB=tfB %>% rowwise() %>% mutate(perm=list(hu2[,which(colnames(hu2)==samples)])) %>% ungroup()
  tfBB= tfBB %>%rowwise() %>% mutate(pval=mean(perm>MO),exp_MO=mean(perm),z=(MO-exp_MO)/sd(perm)) %>% ungroup() %>% dplyr::select(-c(perm))
  
  TF=left_join(tfAA %>% transmute(Antigen,A=SCORE),tfBB %>% transmute(Antigen,B_z=z,B_pval=pval))
  
  #removing infinite SCORE
  if(is.infinite(max(TF$A))==T){
    top3=sort(unique(TF$A),decreasing=T)[1:3]
    TF$A[is.infinite(TF$A)]= 2*top3[2]-top3[3]}
  if(is.infinite(max(TF$B_z))==T){
    top3=sort(unique(TF$B),decreasing=T)[1:3]
    TF$B[is.infinite(TF$B)]= 2*top3[2]-top3[3]}
  TF["z_means"]=rowMeans(scale(TF[,c("A","B_z")]))
  TF=TF %>% arrange(desc(z_means))
  return(TF)
}##luka##

######
##gimme_EC
#####
print('gimme_EC(prot=\"POLR3G\",threshold=0.01,plot=T,neighbourhood=F,legend=T,label=T,
      del=tem,improved_fit=Finstructions=T,main="Z-Zplot",df_improvedfit=15000),
      xlab="Z score",ylab="Reciprocal Z score"')
gimme_EC=function(prot="POLR3G",threshold=0.01,plot=T,neighbourhood=F,df=1000,spar=NULL,legend=T,label=T,
                  del=tem,col="skyblue",improved_fit=F,instructions=T,main="Z-Z plot",df_improvedfit=15000,
                  xlab="Z score",ylab="Reciprocal Z score")
{
  if(instructions==T){print("Make sure you have that rowwise z-normalized EC matrix in the flag del")}
  
#    prot="POLR3G";threshold=0.01;plot=T;neighbourhood=F;df=1000;spar=NULL;legend=T;label=T;del=tem;col="skyblue";improved_fit=F;instructions=T;main="Z-Z plot";df_improvedfit=15000
  
tub=data.frame(self=del[prot,rownames(del)!=prot],other=del[colnames(del)!=prot,prot]) #%>% filter(other>0.7 & other<0.75)
#fitting
if(improved_fit==F){fit <- smooth.spline(tub$self, tub$other,df=df,spar=spar)
tub$predicted <- predict(fit, x = tub$self)$y
tub= tub %>% mutate(distances=other-predicted) %>% arrange(desc(distances))}
else
{
  df=df_improvedfit
  tub=tub %>% arrange(self)
  fit <- smooth.spline(tub$self, tub$other,df=df,spar=spar)
  tub$predicted <- predict(fit, x = tub$self)$y

  #####################Removing those few outliers in the right side
  hya=tub %>% arrange(desc(self)) %>% dplyr::slice(1:50) %>% mutate(diff=abs(self-lead(self))) #finding consecutive gap length in self
  hya= hya %>% mutate(z_diff=(diff-mean(hya$diff[10:49]))/sd(hya$diff[10:49])) %>% filter(z_diff>50) #Removing those points that deviate signifiacnlty from the gap distribution of points that are in position 10 to 50 form right side
  tub2=tub[!rownames(tub) %in% rownames(hya),]
  tub2=tub2[abs(tub2$other-tub2$predicted)<0.01,] #removing outliers
  fit2 <- smooth.spline(tub2$self, tub2$other,df=df,spar=spar)
  tub$predicted <- predict(fit2, x = tub$self)$y
  tub= tub %>% mutate(distances=other-predicted) # %>% arrange(desc(distances))
  
  }


#plotting
if(plot==T){plot(tub$self, tub$other, pch = 19, cex=1,col = col, xlab = xlab, ylab = ylab, main = main)
tub %>% arrange(self) %>% dplyr::select(self,predicted) %>% points(type="l")
if(legend==T){legend("right", legend = c("Original Points", "Fitted Curve"), col = c(col, "black"), pch = c(19, NA), lwd = c(NA, 1),cex=0.5)}

if(label==T)
  {
  tub2=tub %>% dplyr::filter(!grepl("PO00000LR3",rownames(tub))) %>%
  filter(distances>threshold ) 
# which(rownames(tub) %in% nurd_complex)
points(tub2$self, tub2$other, col = "yellow", lwd = 2,pch=19,)
text(x=tub2$self,y=tub2$other,col="black",labels=rownames(tub2),cex=0.5)}



}
###finding neighbourhood
if(neighbourhood==T)
{
  print("To find neighbourhood, make sure you have interaction matrix in variable -> matrix")
  print("also make sure you have loaded finction -> gimme_neighbourhood()")
  print("neighbourhood doesnt work for those 3 proteins that are disconnected from the ppi")
  nebour=gimme_neighbour(matrix = matrix,gene=prot)
  tub["neighbourhood"]=0
  for(i in c(1:length(nebour)))
  {tub$neighbourhood[rownames(tub) %in% nebour[i][[1]]]=i}
}
return(tub)
}


gimme_ECd_avg=function(prot="POLR3G",neighbourhood=F,score_matrix=NULL)
{
  print("Make sure you have ECd_avg variable loaded else pass ecdzmin or whatever you like in parameter score_matrix(previously oki)")
  if(is.null(score_matrix)) {score_matrix=as.data.frame(ECd_avg[,prot])} else
  {
    score_matrix=as.data.frame(score_matrix[,prot])
  }
  colnames(score_matrix)=prot
  
  if(neighbourhood==T)
  {
    print("To find neighbourhood, make sure you have interaction matrix in variable -> matrix")
    print("also make sure you have loaded finction -> gimme_neighbourhood()")
    print("neighbourhood doesnt work for those 3 proteins that are disconnected from the ppi")
    nebour=gimme_neighbour(matrix = matrix,gene=prot)
    score_matrix["neighbourhood"]=0
    for(i in c(1:length(nebour)))
    {score_matrix$neighbourhood[rownames(score_matrix) %in% nebour[i][[1]]]=i}
  }
  return(score_matrix)
  
}
