

#' Vrtx - creates regular polygon
#' @description Returns coordinates of a regular polygon with n vertices. It plots the polygon with the flag "plot=T"
#'
#' @param n an integer that indicates the number of vertex of the polygon
#' @param plot a boolean T/F that indicates wheteher to plot the polygon or not
#'
#' @return t
#' Returns the coordinates of the polygon
#' @export
#' @examples
vrtx=function(n,plot=F,label=NULL,push_label=1)
{
  par(pty="s", mar = 0.5*c(5, 5, 5, 5))
  # n=6
  t=data.frame(x=sapply(0:(n-1),  function(x) sinpi((360/n)*x/180)),
               y=sapply(0:(n-1),  function(x) cospi((360/n)*x/180)))
  if(plot==T){ 
    plot(c(1.2,-1.2),c(1.2,-1.2),col="transparent",xlab = "",ylab="")
    polygon(t)
    if(is.null(label)){label=paste0("V",1:n)}
    for(i in c(1:n)){text(t[i,]*(1+push_label/10),label=label[i],cex = 1)}
    }
  return(t)
}

#' Center of Mass
#' @description Returns coordinates of center of mass from the data provided.
#'
#' @param data a table that provides mass of different objects
#' @param x a vector that provides x coordinates of  objects
#' @param y a vector that provides y coordinates of  objects
#' @return
#' Returns a table with 2 columns(i.e. comx and comy) containing x coordinate and y coordinate of the center of mass
#' @export
#'
#' @examples
getCOM=function(data,x=c(0,1,0.5) ,y=c(0,0,sqrt(3)/2))
{return(data.frame(comx=rowSums(as.matrix(data)%*%diag(x)),comy=rowSums(as.matrix(data)%*%diag(y)))/rowSums(data))}






#' Plotting Triangle with center of mass
#'
#' @param data a table with three columns
#' @param entropyrange a vector with 2 elements in the range [0,1] i.e. c(minimum entropy , maximum entropy)
#' only the entries that has their entropy within this range are highlighted
#' @param magnituderange a vector with 2 elements in the range [0,1] i.e. c(minimum magnitude , maximum magnitude)
#' only the entries that has their magnitude within this range are highlighted
#' @param col a vector of 3 elements that contains color names
#' @param output_table a boolean value of T or F indicating whether to return the output table or not
#' @param plotAll a boolean value
#' If T, all the points that are out of entropyrange or magnitude range are plotted in whitesmoke color, else they are not plotted at all.
#' @param cex a numeric value indicating the size of the points
#' @param pch a integer value indicating different shapes of the points
#'
#' @return
#'  Returns a table with the coordinates of center of mass entropy(in the range of 0 and 1) and magnitude(in the range of 0 and 1)
#' @export
#' @import dplyr
#'
#' @examples
plotTriangle=function(data,
                      entropyrange=c(0,Inf),normalized_entropyrange=c(0,1),
                      magnituderange=c(0,Inf),normalized_magnituderange=c(0,1),
                      maxvaluerange=c(0,Inf),
                      col=c("red","green","blue"),
                      output_table=T,plotAll=T,cex=1,pch=16,
                      title="",push_text=1.2,
                      margin=5)
{
   #require(dplyr)
    n=ncol(data)
    px=vrtx(n)$x
    py=vrtx(n)$y
   #data=log2(data)
    colu=colnames(data)
    colnames(data)=c("a","b","c")
    final=cbind(data,getCOM(data,px,py))
  
   #Finding entropy and magnitude
    final=final %>% dplyr::mutate(entropy=ent(data),
                                entropy_norm=MinMaxNorm2(entropy,0,1))
    final=final %>% dplyr::mutate(magnitude=sqrt(a^2+b^2+c^2),
                                magnitude_norm=MinMaxNorm2(magnitude,0,1))
    
   #Finding max value and assigning color on the basis of max value
    final=final %>% dplyr::mutate(max=pmax(a,b,c),max_norm=MinMaxNorm2(max,0,1)) %>% 
      dplyr::mutate(color=ifelse(a==max,col[1],ifelse(b==max,col[2],col[3])))
    
   #Filtering on the basis of entropy
    final$color[!(final$entropy>=entropyrange[1] & 
                    final$entropy<=entropyrange[2])]="whitesmoke"
    
    final$color[!(final$entropy_norm>=normalized_entropyrange[1] & 
                  final$entropy_norm<=normalized_entropyrange[2] )]="whitesmoke"
    
   #Filtering on the basis of magnitude
    final$color[!(final$magnitude_norm>=normalized_magnituderange[1] & 
                    final$magnitude_norm<=normalized_magnituderange[2])]="whitesmoke"
    final$color[!(final$magnitude>=magnituderange[1] & 
                  final$magnitude<=magnituderange[2])]="whitesmoke"
    
   #Filtering on the basis of maxvalue
    final$color[!(final$max>=maxvaluerange[1] & 
                    final$max<=maxvaluerange[2])]="whitesmoke"
  
   #Plotting
  {
    par(pty="s", mar = 0.5*c(margin, margin, margin, margin))
    plot(c(1.2,-1.2),c(1.2,-1.2),pch=16,xlab = "",ylab="",frame.plot = F,ann = T,axes = F,col="transparent",main=title)
    if(plotAll==T){
      badpoints=final %>% dplyr::filter(color=="whitesmoke")
      points(badpoints$comx,badpoints$comy,col=badpoints$color,pch=pch,cex=cex)}

    goodpoints=final %>% dplyr::filter(color!="whitesmoke") 
   
    points(goodpoints$comx,goodpoints$comy,col=unlist(ifelse(pch!=21,list(goodpoints$color),"black")),
           bg=unlist(ifelse(pch!=21,"black",list(goodpoints$color))),pch=pch,cex=cex)
    polygon(px,py,pch=16)
    text(colu,x=push_text*px,y=push_text*py)
  }##Luka##

    #Returning table
  colnames(final)[1:3]=colu
  if(output_table==T) return(final)
}



plotRope=function(data,
                  push_text=1,rope_width=1,rope_color="#CCCCCCCC",rope_border=T,
                  col=c("red","blue"),col_bg="whitesmoke",pch = c(21,21),pch_bg=19,cex=1,
                  entropyrange=c(0,Inf),maxvaluerange=c(0,Inf),
                  plotAll=T,label=T,title="title",
                  output_table=T
                  )
{
  require(dplyr)
  rope_width=rope_width*0.25
  # data=as.data.frame(abs(matrix(rnorm(nu*2,100,sd=100),nrow=nu,ncol=2)))
  n=ncol(data)
  px=c(-1,1);py=c(0,0)
  colu=colnames(data)
  colnames(data)=c("a","b")
  final=cbind(data,getCOM(data,px,py))
  final$comy=rope_width*sample(c(-500:500),nrow(final),replace = T)/500
  final=final %>% dplyr::mutate(color=ifelse(a>b,col[1],col[2]))
  
  ##Filtering on the basis of maxvalue and entropy
  final=final %>% dplyr::mutate(maxvalue=pmax(a,b),entropy=ent(data))
  final$color[!(final$entropy>=entropyrange[1] & 
                  final$entropy<=entropyrange[2])]=col_bg
  final$color[!(final$max>=maxvaluerange[1] & 
                  final$max<=maxvaluerange[2])]=col_bg
  
  #plotting background
  # par(pty="s", mar = 0.5*c(5, 5, 5, 5))
  plot(0,0,xlim=px*push_text*1.5,pch=16,xlab = "",ylab="",frame.plot = F,ann = T,axes = F,col="transparent",main=title)
  polygon(x=c(-1,-1,1,1),y=rope_width*c(-1,1,1,-1),col=rope_color,border = rope_border)
  if(plotAll==T){
    badpoints=final %>% dplyr::filter(color==col_bg)
    points(badpoints$comx,badpoints$comy,pch=pch_bg,cex=cex,
           col=unlist(ifelse(pch_bg %in% c(21:25),"black",list(badpoints$color))),
           bg=unlist(ifelse(pch_bg %in% c(21:25),list(badpoints$color),"black")))
    }
  
  goodpoints=final %>% dplyr::filter(color!=col_bg) 
  
  #plotting highlighted points
    type1=  goodpoints %>% filter(color==col[1])
    points(type1$comx,type1$comy,pch=pch[1],cex=cex,
           col=unlist(ifelse(pch[1] %in% c(21:25),"black",list(type1$color))),
           bg=unlist(ifelse(pch[1] %in% c(21:25),list(type1$color),"black")))
    type2=  goodpoints %>% filter(color==col[2])
    points(type2$comx,type2$comy,pch=pch[2],cex=cex,
           col=unlist(ifelse(pch[2] %in% c(21:25),"black",list(type2$color))),
           bg=unlist(ifelse(pch[2] %in% c(21:25),list(type2$color),"black")))


        
  if(label==T) text(colu,x=px*push_text*1.2,y=py)
  if(output_table==T) return(final)
  
}




plotTetrahedron=function(data,
                         entropyrange=c(0,Inf),normalized_entropyrange=c(0,1),
                         magnituderange=c(0,Inf),normalized_magnituderange=c(0,1),
                         maxvaluerange=c(0,Inf),
                         col=c("red","green","blue","orange"),
                         labels=NULL,
                         output_table=F,remove_background=T,plot=T,cex=1,pch=16,
                         title="",push_text=1.2,
                         margin=5,
                         edge_width=10,edge_color="black")
{
  #require(dplyr)
  # data=y %>% select("POLR3A","POLR3GL","POLR3C","POLR3G") %>% slice_sample(n=500)
  # col=c("red","green","blue","orange")
  # normalized_magnituderange=c(0.01,1)
  # data=hmm;remove_background = T;normalized_magnituderange = c(0.001,1);output_table = T
  
  #Creating a Tetrahedron and finding center of mass
  {  vertices <- matrix(c( 1,  1,  1,
                           1, -1, -1,
                           -1,  1, -1,
                           -1, -1,  1), ncol = 3, byrow = TRUE)
  com_=data.frame(comx=rowSums(as.matrix(data)%*%diag(vertices[,1])),
                  comy=rowSums(as.matrix(data)%*%diag(vertices[,2])),
                  comz=rowSums(as.matrix(data)%*%diag(vertices[,3])))/rowSums(data)
  vertices11 = as.data.frame(vertices) %>% mutate(labels = colnames(data)) #for plotting vertex
  edges11=as.data.frame(vertices[c(1,2,3,1,4,2,3,4),]) #for plotting edges
  final=cbind(data,com_)
  }##luka##
  
  #Finding entropy,magnitude,max value and color
  {  
    #Finding entropy and magnitude
    final=final %>% dplyr::mutate(entropy=ent(data),
                                  entropy_norm=MinMaxNorm2(entropy,0,1))
    final=final %>% dplyr::mutate(magnitude=sqrt(rowSums(data^2)),
                                  magnitude_norm=MinMaxNorm2(magnitude,0,1))
    
    #Finding max value and assigning color on the basis of max value
    final=final %>% dplyr::mutate(max=apply(data,1,max),max_norm=MinMaxNorm2(max,0,1)) %>% 
      dplyr::mutate(color=ifelse(.[[1]]==max,col[1],ifelse(.[[2]]==max,col[2],ifelse(.[[3]]==max,col[3],col[4]))))
  }##Luka##
  
  #Filtering on the basis of entropy, magnitude and max value
  {
    #Filtering on the basis of entropy
    final$color[!(final$entropy>=entropyrange[1] & 
                    final$entropy<=entropyrange[2])]="whitesmoke"
    
    final$color[!(final$entropy_norm>=normalized_entropyrange[1] & 
                    final$entropy_norm<=normalized_entropyrange[2] )]="whitesmoke"
    
    #Filtering on the basis of magnitude
    final$color[!(final$magnitude_norm>=normalized_magnituderange[1] & 
                    final$magnitude_norm<=normalized_magnituderange[2])]="whitesmoke"
    final$color[!(final$magnitude>=magnituderange[1] & 
                    final$magnitude<=magnituderange[2])]="whitesmoke"
    
    #Filtering on the basis of maxvalue
    final$color[!(final$max>=maxvaluerange[1] & 
                    final$max<=maxvaluerange[2])]="whitesmoke"
    #Adding labels and colors
    if(is.null(labels)) 
    {print("No labels provided. Using rownames as labels");labels=rownames(final)}
    final=final %>% mutate(color=factor(color),labels=labels)  #factorizing color coz thats needed in plotly
  }##Luka##
  
  #Plotting
  {
    
    final_good= final %>% filter(color!="whitesmoke")
    if(remove_background==T) {final_bad=final %>% dplyr::slice()} else
    {   final_bad= final %>% filter(color=="whitesmoke")}
    

    
    if(plot==T){
    now_plotting= plot_ly(
      data = vertices11,
      x = ~V1,y = ~V2,z = ~V3,
      type = 'scatter3d',mode = 'text',
      text = ~labels,textposition = 'top center',
      textfont = list(size = 20, color = col)
    ) %>%    ##Now adding background
      add_trace(
        data = final_bad ,
        x = ~comx,y = ~comy,z = ~comz,
        type = 'scatter3d',mode = 'markers+text',
        text = ~labels,textposition = 'top center',
        textfont = list(size = 10, color=~color,colors=col),
        marker = list(size = 10,color=~color,colors=col )
      )  %>%   ##Now adding frontground
      add_trace(
        data = final_good ,
        x = ~comx,y = ~comy,z = ~comz,
        type = 'scatter3d',mode = 'markers+text',
        text = ~labels,textposition = 'top center',
        textfont = list(size = 10, color=~color,colors=col),
        marker = list(size = 10,color=~color,colors=col )
      ) %>%   ##Now drawing the edges
      add_trace(
        data = edges11,
        x = ~V1,y = ~V2,z = ~V3,
        type = 'scatter3d',mode = 'lines+text',
        line = list(color = edge_color, width = edge_width),
        text=NULL
      )
print(now_plotting)
}
    
  }##Luka##
  
  #Returning table
  if(output_table==T) return(final)
}##luka##







#' MinMax Normalization
#'
#' @param x data that needed to be normalized
#' @param min the minimum value of the final normalized data
#' @param max the maximum value of the final normalized data
#' @return a vector with normalized data
#' @export
#'
#' @examples
MinMaxNorm2=function(x,min,max)
{
  MIN=min(x)
  MAX=max(x)
  y=as.double((x-MIN)*(max-min)/(MAX-MIN)+min)
  return(y)
}



# finding entropy
ent=function(data)
{
  data=data/rowSums(data)
  eata=log2(data)
  fata=eata %>% dplyr::mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))
  return(-rowSums(data*fata))
}



#' Plotting a dummy triangle
#'
#' @export
#' @import dplyr
#' @examples
plotTriangle_dummy=function(n=NULL)
{
  # par(mar=c(1,0.5,1,0.5))
  #require(dplyr)
 
  repeat{
    if(is.null(n)){
    ok=as.integer(readline(prompt=paste0("Enter the number of random points you want to plot or enter 0 to end = ","\n")))
    }else {ok=n}
    
    if(ok==0) {dev.off();break}
    
    if(ok>50000 & is.null(n)){
      ok2=readline(prompt=paste0('Plotting more than 50 thousands points may hang your system. Enter "YES" if you wanna proceed = ',"\n"))
      if(ok2!="YES") next}
    
    if(ok>50000 & (!is.null(n))){
      print("Can't plot more than 50000 points in vignette")
       break}

    data=data.frame(a=runif(ok,0,1),b=runif(ok,0,1),c=runif(ok,0,1))
    data=data %>% dplyr::mutate(color=rgb(a,b,c,maxColorValue = 1))
    data=cbind(data,getCOM(data[,c(1,2,3)]))
    data=data %>% dplyr::arrange(a^2+b^2+c^2)

    #dev.off()
    par(pty="s", mar = 0.5*c(5, 5, 5, 5))
    plot(c(0,1,0.5),c(0,0,sqrt(3)/2),pch=16,xlab = "",ylab="",frame.plot = F,ann = F,axes = F,col="transparent")
    points(data$comx,data$comy,col=data$color,pch=16,cex=0.5)
    polygon(c(0,1,0.5),c(0,0,sqrt(3)/2),pch=16)
    text(0.2,0.85,"Plotting",cex=0.5)
    text(0.2,0.8,ok,cex=0.7)
    text(0.2,0.75,"random points",cex=0.5)
    
    if(!is.null(n)) break
  }
  
  
  
}





