##
source("./CODES/plotTriangle_ver2.R")
library("tibble")
library(plotly)
#####


#############################################
#############################################
############### Grabbing complex from CORUM
##########################################
########https://mips.helmholtz-muenchen.de/corum/

corum=fread("CORUM_5.3.txt")
ct =corum %>% transmute(`#Complex ac`=complex_id, `Recommended name`= complex_name, unique_prots=subunits_gene_name)
ct= ct %>% rowwise() %>% mutate(unique_prots = str_split(unique_prots,";"))
ct=ct %>% rowwise() %>% mutate(num_unique_prots=length(unlist(unique_prots))) %>% ungroup()
ct=ct[sapply(c(1:nrow(ct)), function(x) sum(!ct$unique_prots[[x]] %in% rownames(ECdz_min)))==0,] # Getting rid of some complexes whose one or more subunits are not in my ECd_min because I removed many proteins that didnot have physical but genetic interaction while making ECd_min.
#Now remobing some actual duplicates, then one having duplicate name but distinct subunits.
ct = ct %>%
  distinct(`Recommended name`, unique_prots, num_unique_prots, .keep_all = TRUE) %>%
  arrange(`Recommended name`, desc(num_unique_prots)) %>%  # Highest first
  group_by(`Recommended name`) %>%
  mutate(
    n_versions = n(),
    `Recommended name` = if_else(
      n_versions > 1,
      paste0(`Recommended name`, " ver", row_number()),
      `Recommended name`
    )
  ) %>%
  select(-n_versions) %>%
  ungroup()

saveRDS(ct,"complexinfo_CORUM_human.Rdata")
###!!!!!!Try to transform corum table as close to the previous complexinfo table, sothat you can use same code for downstream processing.





#############################################
######################takes some time like 4 hours or so, but this is one time job####################
############### finding CIM from uncollapsed complex
#############################################
#############################################
{
  ct=readRDS("complexinfo_CORUM_human.Rdata")
  head(ct)
  cmplx_sizes=sort(unique(ct$num_unique_prots))
  backg_mean= matrix(rep(0,length(cmplx_sizes)*length(cmplx_sizes)),nrow=length(cmplx_sizes),ncol=length(cmplx_sizes))
  colnames(backg_mean)=cmplx_sizes
  rownames(backg_mean)=cmplx_sizes
  backg_sd= backg_mean
  
  ECdz_min=readRDS("ECdz_min.Rdata")
  proteins_in_complexes=table(unlist(ct$unique_prots)) %>% as.data.frame() %>% mutate(Var1=as.character(Var1))
  indx=rownames(ECdz_min) %in% proteins_in_complexes$Var1
  temp=ECdz_min[indx,indx]
  diag(temp)=NA
  
  nperms=10000
  for(i9 in c(1:nrow(backg_mean)))
  {
    i=cmplx_sizes[i9]
    cat(paste0("\r ",i9,"/",nrow(backg_mean)))
    for(j9 in c(i9:nrow(backg_mean)))
    {
      j=cmplx_sizes[j9]
      perms= sapply(c(1:nperms),function(x) 
        median(temp[sample(1:nrow(temp),i,replace = F),sample(1:nrow(temp),j,replace = F)],na.rm = T)
      )
      backg_mean[i9,j9]=mean(perms)
      backg_mean[j9,i9]=backg_mean[i9,j9]
      backg_sd[i9,j9]=sd(perms)
      backg_sd[j9,i9]=backg_sd[i9,j9]
      
    }
  }
  saveRDS(backg_mean,"./backg_mean29n.Rdata")
  saveRDS(backg_sd,"./backg_sd29n.Rdata")
  
  
  
  
  
  #######Making median interaction matrix of complexes #Takes about 4 hours or so
  temp= matrix(rep(0,nrow(ct)^2),nrow=nrow(ct))
  rownames(temp)=ct$`#Complex ac`
  colnames(temp)=ct$`#Complex ac`
  
  for(i in c(1:nrow(ct)))
  {
    
    for(j in c(i:nrow(ct)))
    {
      cat(paste0("\r         complexes done:      ",i,"/",nrow(ct),"                 ",j))
      # which(ct$`#Complex ac`==colnames(temp)[i])
      temp2=ECdz_min[ct$unique_prots[[i]],ct$unique_prots[[j]],drop=F]
      if(0 %in% dim(temp2)) next
      temp[i,j]=(median(temp2,na.rm = T)-backg_mean[as.character(nrow(temp2)),as.character(ncol(temp2))])/backg_sd[as.character(nrow(temp2)),as.character(ncol(temp2))]
      temp[j,i]=temp[i,j]
    }
  }
  
  CIM=temp
  rownames(CIM)=ct$`Recommended name`
  colnames(CIM)=ct$`Recommended name`
  saveRDS(CIM,"CIM_corum.Rdata")
  
}##luka##

#############################################
######################Adding custom made complexes ####################
############### in complexinfo table(ct) and CIM
#############################################
#############################################
{
  CIM=readRDS("CIM_corum.Rdata")
  backg_mean=readRDS("./backg_mean29n.Rdata")
  backg_sd= readRDS("./backg_sd29n.Rdata")
  
  # Write your custom made complexes here.
  rna_pol_complexes=list(
  RNA_Pol1= c("POLR1A","POLR1B","POLR1E","CD3EAP","TWISTNB","ZNRD1"),
  RNA_Pol2= c("POLR2A","POLR2B","POLR2C","POLR2D","POLR2G","POLR2I","POLR2J","POLR2M"),
  RNA_Pol3= c("POLR3A","POLR3B","POLR3C","POLR3D","POLR3E","POLR3F","POLR3G","POLR3H","POLR3K","CRCP"),
  RNA_Pol3_GL= c("POLR3A","POLR3B","POLR3C","POLR3D","POLR3E","POLR3F","POLR3G","POLR3H","POLR3K","CRCP")
  )
  

  
  custom_complex=rna_pol_complexes
  
  # Checking if permutaton is already done for the size of your custom made complexes.
  print(ifelse(sum(!sapply(custom_complex, length) %in% as.integer(rownames(backg_mean)))==0,"Great ! All your custom made complex size are in permutation table.", "Stop !! One or more of your custom made complex has a size which is missing in the permutation table"))
  
  
#Adding your complex in complexinfo table
  temp_ct=ct
  temp_CIM=CIM
  
  for(i in seq_along(custom_complex))
  {
    
    print(i)
    head(temp_ct)
    temp_ct=rbind(temp_ct,
          list(which(is.na(match(c(1:9999),temp_ct$`#Complex ac`)))[1],
               names(custom_complex)[i],
               unname(custom_complex[i]),
               length(unlist(custom_complex[i])))) #adding a row in complexinfo table
    
    
    
    temp_CIM=rbind(temp_CIM,0);rownames(temp_CIM)[nrow(temp_CIM)]=names(custom_complex)[i]
    temp_CIM=cbind(temp_CIM,0);colnames(temp_CIM)[ncol(temp_CIM)]=names(custom_complex)[i]
    
    
      i2=nrow(temp_ct)
      for(j2 in c(1:nrow(temp_ct)))
      {
        cat(paste0("\r         complexes done:      ",i2,"/",nrow(temp_ct),"                 ",j2))
        temp2=ECdz_min[temp_ct$unique_prots[[i2]],temp_ct$unique_prots[[j2]],drop=F]
        if(0 %in% dim(temp2)) next
        temp_CIM[i2,j2]=(median(temp2,na.rm = T)-backg_mean[as.character(nrow(temp2)),as.character(ncol(temp2))])/backg_sd[as.character(nrow(temp2)),as.character(ncol(temp2))]
        temp_CIM[j2,i2]=temp_CIM[i2,j2]
      }

    
  }
  
 
  
  
  
  
  }##Luka##  

saveRDS(temp_CIM,"CIM_CORUM+custommadepolymerase.Rdata")





