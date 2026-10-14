
require(MASS)
require(entropy)
require(igraph)
require(dplyr)
# require(dominatR)

eval(parse(file="./../../CODES/my_function.R"))
eval(parse(file="./../../CODES/plotTriangle_ver2.R"))

##################################################
##################################################
###############Simulation
##################################################
##################################################

################################STRATEGY2
##Every protein have equal P, but unequal number of experiments(poopularity)
sz=1000
prob= 50/sz # Every protein has 20 interacting protein out of 1000 protein
num_exp=as.integer( MinMaxNorm(c(1:sz))*200)+1
num_exp=c(1:1000)
# num_exp=as.integer(MinMaxNorm(2^(c(1:sz)/100))*200+MinMaxNorm(c(1:sz))*20+1)
plot(num_exp)
num_exp     #OR MIGHT TRY THE pattern that purely follows simple power law rather than complex mixture of power law
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

################################STRATEGY1
##No experiments, but different power of interactions
sz=1000
sim=matrix(data=0,nrow=sz,ncol=sz)
# num_exp=as.integer( MinMaxNorm(c(1:sz))*200)+1
num_exp=c(1:sz) #^2
for(i in c(1:sz))
{
 sim[i,]=sim[i,]+num_exp[i]
 sim[,i]=sim[,i]+num_exp[i]
}

###########################Naming rows and columns
rownames(sim)=paste("P",c(1:sz),sep = "")
colnames(sim)=rownames(sim)

###########################Introducing bias
if(F){
sim["P440","P100"]=sim["P440","P100"]+1000
sim["P100","P440"]=sim["P100","P440"]+1000
}



if(F){
  # sim=sim[c(1:500,900:1000),c(1:500,900:1000)]
  # for(i in c(500:900)){ sim[i,]=sim[500,];sim[,i]=sim[,500]}
  indx=sort(sample(c(1:5000),200));sim=sim[indx,indx]
  } #playing around with the data

##################################################
##################################################
###############FINDING EFFECTIVE CONDUCTANCE
##################################################
##################################################

##########################Finding Admittance and degree of graph
A_s=sim
diag(A_s)=0
D_s=A_s*0
diag(D_s)=rowSums(A_s)
#########################Finding Laplacian
L_s=D_s-A_s  #Computing Normal Laplacian
if(F){D_inv_sqrt <- diag(1 / sqrt(diag(D_s)))  # Compute D^(-1/2)
L_s <- D_inv_sqrt %*% L_s %*% D_inv_sqrt  # Computing normalized laplacian
}
#########################Finding PseudoInverse
z_s=ginv(L_s)
# saveRDS(z_s,"z_s_5000.Rdata")
#Final Treatment
FF_s=z_s
for (j in 1:nrow(FF_s)) 
{
  if(j%%100==0){print(j)} ##checking progress of loop
  for (k in 1:nrow(FF_s)) 
  {
    FF_s[j, k] = z_s[j, j] + z_s[k, k] - z_s[j, k] - z_s[k, j]
  }
}
FFF_s=1/FF_s
colnames(FFF_s)=colnames(sim);rownames(FFF_s)=rownames(sim)
FFF_s[is.infinite(FFF_s)]=0
#Z normalization on effective conductance
tem_s=FFF_s
tem_s[is.infinite(tem_s)]=0
for(i in c(1:nrow(tem_s))){tem_s[i,]=scale(tem_s[i,])}

#Z normalization on conductance (OBSOLETE) in simulation data
tem_s2=sim
tem_s2[is.infinite(tem_s2)]=0
for(i in c(1:nrow(tem_s2))){tem_s2[i,]=scale(tem_s2[i,])}


#Z normalization on conductance (OBSOLETE) in real  data
if(T){tem_s2=matrix
tem_s2[is.infinite(tem_s2)]=0
for(i in c(1:nrow(tem_s2))){tem_s2[i,]=scale(tem_s2[i,])}
self=tem_s2["POLR3G",]
other=tem_s2[,"POLR3G"]
plot(self,other,col="#11EE2211",pch=19)
text(x=self,y=other,label=rownames(tem_s2),cex=0.3)
plot(other,pch=19,col="#11EE2211")
text(x=c(1:nrow(tem_s2)),y=other,label=rownames(tem_s2),cex=0.3)
table(other>0)
table(matrix["MYC",]>0)}

#Analyzing sd and variance of each protein across popularity
hq=sim
hq[is.infinite(hq)]=0
par(mfrow=c(3,1))
mean_s=sapply(c(1:1000),function(x) mean(hq[x,]))
var_s=sapply(c(1:1000),function(x) var(hq[x,]))
plot(mean_s)
plot(var_s)
plot(sim[8,]-mean_s)
gimme_EC("P300",del = tem_s)

######################################################################
######################################################################
#######################PLOTTINGGGGGG..................................
######################################################################
######################################################################
#########Plotting P1to P1000
dev.off()
for(i in paste0("P",c(1:(sz/100))*100))
  # for(i in colnames(sim))
{
  print(i)
  if(!i %in% colnames(tem_s)) next
  #Sorted effective conductance for strategy2, normalized laplacian simulation
  # plot(scale(sort(FFF_s[i,])));abline(h=0.501)
  # plot(tem_s[i,],tem_s[,i]);abline(h=0.501)
  
  ###not normalized conductance/effective conductance
   # plot(sim[i,],sim[,i],ylim=c(min(sim),max(sim)),xlim=c(min(sim),max(sim))); Sys.sleep(0.1);next
   ##not normalized effective conductance
   # plot(FFF_s[i,],FFF_s[,i],
   #      ylim=c(min(FFF_s),unique(sort(FFF_s,decreasing = T))[2]),
   #      xlim=c(min(FFF_s),unique(sort(FFF_s,decreasing = T))[2])); Sys.sleep(0.1);next
  ###conductance z normalized
   # wew=gimme_EC_s(i,df=100,label=F,legend=F,del=tem_s2,xlim=NULL,ylim=NULL)
   # wew=gimme_EC_s(i,df=100,label=F,legend=F,del=tem_s2,xlim=NULL,ylim=c(-2,2))
   #Effective conductance z normalized
 if(F){  wew=gimme_EC_s(i,df=49,label=F,legend=F,del=tem_s,xlim=c(-2,2),ylim=NULL)
   abline(v=MinMaxNorm2(c(1,as.integer(substr(i,2,99)),1000),-2,1.5))}
   # wew=gimme_EC_s(i,df=100,label=F,legend=F,del=tem_s,xlim=NULL,ylim=c(-2,2))

    #2 way plots below
  #first
  wew=gimme_EC(i,df=100,label=F,legend=F,del=tem_s,plot = T)
  
  #second  
  # wew=gimme_EC_s(i,del=sim,df=100,label = F,col="#AF2FFA02",legend=F,boxplot="P1000")  ####TRY to make animation of FFF_s, and sim too
   # wew=gimme_EC_s(i,del=FFF_s,df=100,label = F,col="#AF2FFA02",legend=F,xlim=c(200000,800000),ylim=c(200000,800000))  ####TRY to make animation of FFF_s, and sim too
  # hya=sort(wew$distances)
 Sys.sleep(0.2)
}
dev.off()



######################################################################
#######################PLOT Comparison..................................
######################################################################
bias=readRDS("bias.Rdata")
ew=bias %>% mutate(total=A+B,total_n=as.integer(MinMaxNorm2(log(1+total),1,1000))) %>% arrange(total_n)
# par(mfrow=c(3,1))
wew2=gimme_EC_s(as.character(ew$Var1[1*20]),df=100,label=F,legend=F,del=tem,xlim=NULL,ylim=NULL,plot=F) %>% arrange(self)#strategy2
par(mfrow=c(3,1))
Sys.sleep(10)
print("Get> Ready> Record..... in 2 second")
Sys.sleep(2)

for(i in c(1:1000))
{
  print(i)
  if(i%%20!=0) next
  prot=paste0("P",i)
  wew=gimme_EC_s(prot,df=500,label=F,legend=F,del=tem_s_1,xlim=NULL,ylim=NULL,spar=0.8,main="Simulation Strategy1") #strategy1
  wew=gimme_EC_s(prot,df=500,label=F,legend=F,del=tem_s_2,xlim=NULL,ylim=NULL,spar=0.8,main="Simulation Strategy2") #strategy2
  
  plot(wew2$self,wew2$other,pch=19,col="skyblue",main="BIOGRID data",xlab="Self",ylab="Other")
  points(wew2$self,wew2$predicted,type="l")
  i2=which(ew$total_n>=i)[1]
  print(i2)
  wew2=gimme_EC_s(as.character(ew$Var1[i2]),df=10000,label=F,legend=F,del=tem,xlim=NULL,ylim=NULL,plot=F,spar=0.8) %>% arrange(self) #strategy2
  Sys.sleep(0.1)
  # dev.off()
}
dev.off()




################################################################################
######################################## Plotting sine counterpart
################################################################################
#########Plotting sin counterpart
for(p in c(1:1000))
{
  xc=c(-200:200)/100
  xc
  plot(xc,sapply(xc,function(x) 1*sin(x/3-p/10)))
  abline(v=)
  Sys.sleep(0.1)}

########Better way to plot sin counterpart
x=MinMaxNorm2(c(1:1000),0,2*pi) #1000 values from 0 to 2pi 
for(i in c(1:1000))
{
  if(i%%20!=0) next #Making plot faster
  print(i)
  window_center=-x[i]
  window_width=pi/8
  xc= window_center + (c(-500:500)/500)*window_width #Adding width of pi/4 to window center
  
  plot(xc,sapply(xc,function(x) 1*sin(x)),col="skyblue")
  Sys.sleep(0.1)
  }

#########Plotting real data
bias=readRDS("bias.Rdata")
ew=bias %>% mutate(total=A+B) %>% arrange(total)
lim=c(min(tem),max(tem))
for(i in c(1:nrow(ew)))
{
  if(i%%430!=0) next
  noww=as.character(ew$Var1[i])
  print(paste(i,noww))
  print(ew[i,])
  hmm=gimme_EC_s(noww,legend = F,neighbourhood = F,del = tem,threshold = 0.001,xlim=NULL,ylim=NULL,label=F)
  text(0,(max(hmm$other)+min(hmm$other))/1.5,label=noww)
  Sys.sleep(2)
}

gimme_EC("ABCA9",label=F,legend=F,neighbourhood = F,del=tem)

##############################################################################################################################
##############################################################################################################################
##########################       IMPROVING SPLINE FUNCTION           #########################################################
##############################################################################################################################
##############################################################################################################################
# require(mgcv)
for(i in c(1500:nrow(ew)))
{
if(i%%100!=0) next
  # KCNC1 TIMM21 KCNC2 CD177 ENPP2 
  #CTRB2(ctrb1 and this dont have many other interaction) 
  #TP53(right points)
# prot=as.character(ew$Var1[i])  
print(ew[i,])
del=tem
spar=NULL
df=15000
tub=data.frame(self=del[prot,rownames(del)!=prot],other=del[colnames(del)!=prot,prot]) #%>% filter(other>0.7 & other<0.75)

#######################spline fit
if(T){
tub=tub %>% arrange(self)
fit <- smooth.spline(tub$self, tub$other,df=df,spar=spar)
tub$predicted <- predict(fit, x = tub$self)$y

######################Removing outliers 
tub2=tub[abs(tub$other-tub$predicted)<0.01,] #removing outliers
#####################Removing those few outliers in the right side
hya=tub %>% arrange(desc(self)) %>% slice(1:50) %>% mutate(diff=abs(self-lead(self))) #finding consecutive gap length in self
hya= hya %>% mutate(z_diff=(diff-mean(hya$diff[10:49]))/sd(hya$diff[10:49])) %>% filter(z_diff>50) #Removing those points that deviate signifiacnlty from the gap distribution of points that are in position 10 to 50 form right side
tub2=tub2[!rownames(tub2) %in% rownames(hya),]
fit2 <- smooth.spline(tub2$self, tub2$other,df=df,spar=spar)
tub$predicted2 <- predict(fit2, x = tub$self)$y

}

print(matrix[prot,matrix[prot,]>0] )

tub= tub %>% mutate(distances=other-predicted,distances2=other-predicted2) %>% arrange(desc(distances))
plot(tub$self, tub$other, pch = 19, cex=1,col = "skyblue", xlab = "Self",
                 ylab = "Other", main = "Curve Fitting Example",ylim=NULL)
tub %>% arrange(self) %>% dplyr::select(self,predicted) %>% points(type="l")
tub %>% arrange(self) %>% dplyr::select(self,predicted2) %>% points(type="l",col="red",lty=2)
#labels
tub3=tub %>% filter(abs(distances2)>0.01)
if(nrow(tub3)>0){text(x=tub3$self,y=tub3$other,labels=rownames(tub3),cex=0.5)}
text(0,(max(tub$other)+min(tub$other))/1.5,label=prot)
Sys.sleep(3)
dev.off()
}
  print("wait")






################Other fit that I had tried
  #lowess fit
  if(F){fit <- loess(other ~ self, data = hmm, span = 0.4, family = "symmetric") #loess fit
  tub$predicted <- predict(fit, hmm$self)
  }
  # --- Alternative: GAM Fit (More Robust) ---
  if(F) {
    # fit <- gam(other ~ s(self, bs="cs",k=50,fx=T), data = tub)
    # fit <- gam(other ~ s(self, bs="cs") + s(x, bs="cs"), data = tub)
    # fit <- gam(other ~ s(self, bs="cs", k=50, fx=TRUE), data = tub,
    #            family = gaussian(link="identity"), method = "REML")
    # fit <- gam(
    #   other ~ s(self, bs="cs", k=30, fx=TRUE),  # Force boundary conditions to match observed data
    #   data = tub, 
    #   method = "REML",gamma=1.5
    # )
    # fit <- gam(
    #   other ~ s(self, bs = "cs", k = 50),   # Flexible cubic spline with 50 basis functions
    #   data = tub,
    #   weights = tub$weights,  # Use weights to handle outliers
    #   method = "REML", 
    #   family = gaussian(link = "identity")
    # )
    
    
    tub$predicted <- predict(fit, newdata = tub)
  }








require(nortest)
wew=gimme_EC("P440",del=tem_s,threshold = 0.001,label = F,df=100)
hist(wew$distances,breaks=10000,xlim=c(-1,1)*(10^-8))
sd(wew$distances)
mean(wew$distances)
abline(v=0.1)
hya=sort(wew$distances)
# hya=log10(1+sort(wew$distances))
# hya=sample(wew$distances,size=10,replace = F)
ad.test(hya)
ad.test(rnorm(1000000,mean=0,sd=1))
qqnorm(hya)
qqline(hya, col="red")
shapiro.test(hya)

wew=gimme_EC("POLR3D",df=100)
hist(wew$distances,breaks=1000,xlim=c(-0.01,0.01))


hist(FFF_s["P800",],breaks=1000)



gimme_EC_s=function(prot="POLR3G",threshold=0.01,plot=T,neighbourhood=F,EC,df=1000,spar=NULL,
                    legend=T,label=T,del=tem,col="skyblue",
                    xlim=c(-5,5),ylim=c(-2,2),boxplot=NULL,main="Z-Z plot")
{
  print("Make sure you have that rowwise z-normalized EC matrix in the flag del")
  tub=data.frame(self=del[prot,rownames(del)!=prot],other=del[colnames(del)!=prot,prot]) #%>% filter(other>0.7 & other<0.75)
  fit <- smooth.spline(tub$self, tub$other,df=df,spar=spar)
  tub$predicted <- predict(fit, x = tub$self)$y
  tub= tub %>% mutate(distances=other-predicted) %>% arrange(desc(distances))
  #plotting
  if(plot==T){plot(tub$self, tub$other, pch = 19, cex=1,col = col, xlab = "Self",
                   ylab = "Other", main = main,
                   xlim=xlim,ylim=ylim)
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
  
  if(!is.null(boxplot))
  {
    points(rep(tem_s[prot,boxplot],1000),unlist(tem_s[boxplot,]),col="green",pch=19)
    points(rep(tem_s[prot,"P1"],1000),unlist(tem_s["P1",]),col="green",pch=19)
    points(rep(tem_s[prot,"P250"],1000),unlist(tem_s["P250",]),col="green",pch=19)
    points(rep(tem_s[prot,"P500"],1000),unlist(tem_s["P500",]),col="green",pch=19)
    points(rep(tem_s[prot,"P750"],1000),unlist(tem_s["P750",]),col="green",pch=19)
    
  }
  return(tub)
}
gimme_EC_s("P10",del=del_s,plot = F,df=100)






####Finding eqn
# Load required libraries
# Load required libraries
{
  library(nls.multstart)
library(dplyr)

# Define the curve model
curve_model <- function(x, a, b, c, k) {
  x + a * x^2 + b * x^3 + c * sin(k * x)
}

# Initialize storage for results
fitted_params <- data.frame(
  popularity = numeric(),
  a = numeric(),
  b = numeric(),
  c = numeric(),
  k = numeric()
)

# Loop through popularity levels (P1 to P1000)
for (pop in 1:1000) {
  # Generate the data for the current popularity level
  print(pop)
  popularity_label <- paste0("P", pop)
  wew <- gimme_EC_s(popularity_label, del = del_s, plot = FALSE, df = 100)
  
  # Try fitting the model, skip if it fails
  tryCatch({
    fit <- nls_multstart(
      other ~ curve_model(self, a, b, c, k),
      data = wew,
      start_lower = c(a = -1, b = -1, c = -1, k = 0.1),
      start_upper = c(a = 1, b = 1, c = 1, k = 10),
      iter = 500
    )
    
    # Extract parameters
    params <- coef(fit)
    fitted_params <- rbind(fitted_params, 
                           data.frame(popularity = pop, 
                                      a = params["a"], 
                                      b = params["b"], 
                                      c = params["c"], 
                                      k = params["k"]))
  }, error = function(e) {
    message(paste("Failed to fit model for popularity level:", pop))
  })
}

# Save fitted parameters for inspection
write.csv(fitted_params, "fitted_params.csv", row.names = FALSE)

# Function to reconstruct the equation for a given popularity
generate_equation <- function(a, b, c, k) {
  paste0("y = x + (", round(a, 4), ") * x^2 + (", round(b, 4), 
         ") * x^3 + (", round(c, 4), ") * sin(", round(k, 4), " * x)")
}

# Add equations to the fitted parameters
fitted_params <- fitted_params %>%
  rowwise() %>%
  mutate(equation = generate_equation(a, b, c, k))

# Save equations to a file for inspection
write.csv(fitted_params, "equations_with_popularity.csv", row.names = FALSE)

# Print the first few equations
print(head(fitted_params))



ggplot(fitted_params, aes(x = popularity)) +
  geom_line(aes(y = a, color = "a")) +
  geom_line(aes(y = b, color = "b")) +
  geom_line(aes(y = c, color = "c")) +
  geom_line(aes(y = k, color = "k")) +
  labs(title = "Parameter Trends Across Popularity Levels",
       x = "Popularity",
       y = "Parameter Value") +
  scale_color_manual(name = "Parameter", values = c("a" = "blue", "b" = "red", "c" = "green", "k" = "purple")) +
  theme_minimal()
}##luka

####Finding eqn in sin(Ax+C) curve
{
  # Initialize storage for fitted parameters
  fitted_params_sin <- data.frame(
    popularity = numeric(),
    B=numeric(),
    C = numeric()
  )
  
  # Fit the model for each popularity level
  for (pop in 1:1000) {
    # Generate data
    popularity_label <- paste0("P", pop)
    wew <- gimme_EC_s(popularity_label, del = del_s, plot = FALSE, df = 100)
    
    # Try fitting the sinusoidal model
    tryCatch({
      fit <- nls(
        other ~ B*sin(2*pi*self + C),
        data = wew,
        start = list(B = 1, C=pi),  # Initial guesses
        algorithm = "port"
      )
      
      # Extract parameters
      params <- coef(fit)
      fitted_params_sin <- rbind(fitted_params_sin, 
                                 data.frame(popularity = pop, 
                                            B = params["B"],
                                            C = params["C"]))
    }, error = function(e) {
      message(paste("Failed to fit model for popularity level:", pop))
    })
  }
  
  # Save results
  write.csv(fitted_params_sin, "fitted_params_sin.csv", row.names = FALSE)
  
  # View first few results
  head(fitted_params_sin)
  
  ###VISUALIZING
  
  ggplot(fitted_params_sin, aes(x = popularity)) +
    geom_line(aes(y = C, color = "C")) +
    geom_line(aes(y = B, color = "B")) +
    labs(title = "Trends in A and C Across Popularity Levels",
         x = "Popularity",
         y = "Parameter Value") +
    scale_color_manual(name = "Parameter",
                       values = c( "C" = "red","B"="green")) +
    theme_minimal()+geom_vline(xintercept=813)+
    geom_vline(xintercept=194)+
    geom_line(aes(data=NULL,x=c(1:1000),
                  y=sapply(c(1:1000),function(x) 1.137^(-(x+10)/100)-0.55)),lty=3)
  
  plot((fitted_params_sin$C))
  plot((fitted_params_sin$B))
  plot(fitted_params_sin$C/fitted_params_sin$B)
  plot(fitted_params_sin$B/fitted_params_sin$C,
       type="l",col="red");points(x=c(1:1000),y=sapply(c(1:1000),
                                                       function(x) 1.14^(-(x-460)/600)-0.99),type="l")
    # plot(1,1,ylim=c(-0.1,0.1)*10,xlim=c(0,1000)*10)
    
  
  
  ten=fitted_params_sin$C/fitted_params_sin$B
  
  plot((c(ten,0)-c(0,ten))[2:999])
  
  plot(c(-5000:5000)/10000,-sapply(c(-5000:5000)/10000,function(x) cos(x)/sin(x)))
  plot(c(-5000:5000)/10000,-sapply(c(-5000:5000)/10000,function(x) 1/x))
  
  
}##luka
