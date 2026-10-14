library(httr)
library(jsonlite)
library(dplyr)
library(IRanges)
library(data.table)
library(GO.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)


#############################################
#############################################
############### Grabbing strucutral info of all proteins
#################### from uniprot REST API website
#############################################
#############################################
h1=bias %>% mutate(change=gained-lost)
#Analayzing the domains in the protein
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
id_map <- getBM(attributes = c("hgnc_symbol", "uniprotswissprot"),filters = "hgnc_symbol",values = h1$Var1,mart = mart)
# Remove rows where uniprotswissprot is blank or NA
id_map <- id_map[!is.na(id_map$uniprotswissprot) & id_map$uniprotswissprot != "", ]
h1 <- merge(h1, id_map, by.x = "Var1", by.y = "hgnc_symbol", all.x = TRUE)
# Initialize with total_length, idr_length, idr_num
h1 <- h1 %>% mutate(total_length = NA)
h1_freq=h1
h1_length=h1
metdat= data.frame(type="",description="")[0,]
  
# Loop over all proteins
for (i in seq_len(nrow(h1))) {
  
   # i=which(h1$Var1=="TP53")
    
  cat(paste0("\r ", i, "/", nrow(h1),"    ||   Number of features = ",ncol(h1_freq)-12,"   ||   metdata so far = ",nrow(metdat)))
  
  #Does it have uniprot id?
  uniprot_id <- h1$uniprotswissprot[i]
  if (is.na(uniprot_id) || uniprot_id == "") next
  
  # Does it have data in uniprot REST API website?
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".json")
  res <- try(GET(url, timeout(60)), silent = TRUE)
  if (inherits(res, "try-error") || status_code(res) != 200) next
  dat <- try(fromJSON(rawToChar(res$content)), silent = TRUE)
  if (inherits(dat, "try-error") || is.null(dat)) next
  
  # How long protein is ??
  seq_len <- dat$sequence$length
  if (is.null(seq_len) || seq_len == 0) next
  h1_freq$total_length[i] <- seq_len
  h1_length$total_length[i] <- seq_len
  
  # What are its features.
  feats <- dat$features
  if (is.null(feats) || nrow(feats) == 0) next
  metdat=rbind(metdat,feats %>% dplyr::group_by(type,description) %>% summarize(.groups = "drop")) #Adding  features to metadata table
  
  # Lets remove unnecessary features
  feats <- feats %>%
    filter(!type %in% c("Sequence conflict","Natural variant","Chain","Site",
                        "Mutagenesis","Alternative sequence"))

  #Copying ligand position to the description
  indx=(feats$type=="Binding site")
  feats$description[indx]=feats$ligand$name[indx]
  
  # For Modified residue: take only before ";"
  feats$description <- sapply(strsplit(feats$description, ";"), `[`, 1)
  
  #Making FEATS,starts and end column
  feats$FEATS= ifelse(is.na(feats$description),feats$type,paste(feats$type,feats$description,sep="__"))
  feats$start=feats$location$start$value;  feats$end=feats$location$end$value

  feats_clean <- feats %>%
    filter(!is.na(start) & !is.na(end)) %>%
    dplyr::select(FEATS, start, end) %>%
    group_by(FEATS) %>%
    summarise( freq=n(),
      merged_length = {
        merged <- IRanges::reduce(IRanges(start, end))
        sum(width(merged))
      },
      .groups = "drop"
    )
  

  
  # Frequency of all remaining features for this protein

  for (j in seq_len(nrow(feats_clean)))
    {
    feat_name <- feats_clean$FEATS[j]
    freq_val <- feats_clean$freq[j]
    merged_length_val <- feats_clean$merged_length[j]
    
    # Create column if it doesn't exist
    if (!(feat_name %in% colnames(h1_freq))) {
      h1_freq[[feat_name]] <- 0
      h1_length[[feat_name]] <- 0
    }
    
    # Assign count for this protein
    h1_freq[i, feat_name] <- freq_val
    h1_length[i, feat_name] <- merged_length_val
    
  }
}
saveRDS(h1_freq,"h1_freq.Rdata") #table that has info about the frequency of the strucutural features.
saveRDS(h1_length,"h1_length.Rdata") #table that has info about the total length of the structural features.
saveRDS(metdat,"metdat.Rdata") #Metdata that contains all the types and description of the structural features of proteins.




#############################################
#############################################
############### Find which features are enriched in the hub proteins
#############################################
#############################################
##Doing permutation analysis to find which featuers are enriched in hub proteins.
{
h1_freq=readRDS("h1_freq.Rdata")  %>% mutate(final_degree=initial_degree-lost) ##Here you decide to take which degree !!
h2=h1_freq[,c(8,13:ncol(h1_freq))]
h2 <- h2[(!is.na(h2$final_degree)) & (!is.na(h2$total_length)), ]
bginfo=data.frame(sample_size=sort(unique(colSums(h2[,3:ncol(h2)]>0))),
                  sd=0,mean=0)
for(i in c(1:nrow(bginfo)))
{
  print(i)
  perms=sapply(c(1:100000),function(x) mean(sample(h2$final_degree,size = bginfo$sample_size[i])))
  bginfo$sd[i]=sd(perms)
  bginfo$mean[i]=mean(perms)
}

h3=data.frame(size=colSums(h2[3:ncol(h2)]>0),obs=0,z=0)
for(i in c(1:nrow(h3)))
{
  print(i)
  h3$obs[i]=mean(h2$final_degree[h2[,rownames(h3)[i]]>0])
  j=which(bginfo$sample_size==h3$size[i])
  h3$z[i]=(h3$obs[i]-bginfo$mean[j])/bginfo$sd[j]
}

 # write.table(h3,"features_enrichment_init-lost_100kperms.tsv",sep="\t")

ew=h3 %>% arrange(z) %>% mutate(indx=c(1:n()))
par(mar = c(6, 6, 4, 6),xpd=T)
plot(ew$indx,ew$z,pch=19,col="whitesmoke",type="p",main="All proteins")
eww=ew %>% filter(abs(z)>=5 & size>=5)
text(eww$indx,eww$z,labels = rownames(eww),cex=0.4)

head(ew,n=30)
tail(ew,n=30) %>% arrange(-z)

# h3=read.table("features_enrichment.tsv")
ew=h3 %>% filter(abs(z)>8)  %>% mutate(feature=rownames(.))
ew$feature <- factor(ew$feature, levels = ew$feature[order(ew$z)])

ggplot(ew, aes(x = feature, y = z, fill = z > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  coord_flip() +   # flip to put categories on y-axis
  labs(x = "", y = "Z score") +
  theme_minimal() +
  theme(legend.position = "none")

}##luka##


####Doing permutation in sub groups of proteins
{
metdat=readRDS("metdat.Rdata")
h1_freq=readRDS("h1_freq.Rdata")  %>% mutate(final_degree=initial_degree-lost) ##Here you decide to take which degree !!
metdat2= metdat[grepl("DNA",metdat$type,ignore.case = T),] %>% distinct()
go2gene <- AnnotationDbi::select(org.Hs.eg.db,keys = keys(org.Hs.eg.db, keytype = "GO"),columns = c("GO", "SYMBOL"),keytype = "GO")

#############TRanscription Factors subgroup
url <- "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
tfs_symbols=read.table(url, header = FALSE, quote = "", comment.char = "",stringsAsFactors = FALSE)[,1]
tfs=h1_freq %>% filter(Var1 %in% tfs_symbols) %>% filter(!is.na(final_degree)) %>% filter(!is.na(total_length))



############Subcellular components group
#s#locations
subcell=fread("subcellular_locations.tsv")
subcell=subcell %>% rowwise() %>% mutate(median_degree= median( h1_freq$final_degree[h1_freq$Var1 %in% (go2gene %>% filter(GO==`GO term`) %>% dplyr::select(SYMBOL) %>% unlist(use.names = F))],na.rm = T))
shortlisted_symbols=go2gene %>% filter(GO==subcell$`GO term`[subcell$`Subcellular location`=="Mitochondria"]) %>%
  dplyr::select(SYMBOL) %>% unlist(use.names = F) %>% unique()
tfs=h1_freq %>% filter(Var1 %in% shortlisted_symbols) %>% filter(!is.na(final_degree)) %>% filter(!is.na(total_length))



############RNA binding protein
shortlisted_symbols=go2gene %>%
  filter(GO %in% c("GO:0003723","GO:0003729","GO:0019843","GO:0003727")) %>% #RNA binding, mrna binding, rRNA binding,snRNA binding
  dplyr::select(SYMBOL) %>% unlist(use.names = F) %>% unique()
# shortlisted_symbols=fread("~/Downloads/Table_HS_RBP.txt") %>% filter(RBP2GO_Score>=3) %>% dplyr::select(Gene_Name) %>% unlist(use.names = F) %>% unique() #https://rbp2go.dkfz.de/
rbp_symbols=fread("~/Downloads/High_confidence_human_RBPs.txt") %>% dplyr::select(Uniprot_ID) %>% unlist(use.names = F) %>% unique() # https://rbp2go-2-beta.dkfz.de/
mapped <- AnnotationDbi::select(org.Hs.eg.db,keys = rbp_symbols,keytype = "UNIPROT",columns = c("SYMBOL"))
rbp_symbols=unique(mapped$SYMBOL)
tfs=h1_freq %>% filter(Var1 %in% rbp_symbols) %>% filter(!is.na(final_degree)) %>% filter(!is.na(total_length))
# shortlisted_symbols= fread("~/Downloads/RBPDB_v1.3.1_human_2012-11-21_CSV/RBPDB_v1.3.1_proteins_human_2012-11-21.csv") %>%
#   dplyr::select(V5) %>% unlist(use.names = F) %>% unique()


drbp_symbols=intersect(tfs_symbols,rbp_symbols)
###############translational factors
tls_symbols=go2gene %>%
  filter(GO %in% c("GO:0003743","GO:0003746","GO:0003747")) %>% #initiaition , elongation and release factor
  dplyr::select(SYMBOL) %>% unlist(use.names = F) %>% unique()
tls_symbols= fread("~/Downloads/tlf_genes_geneid.txt") %>% pull(V1) # https://compbio.uth.edu/TranslFac/download.html
tfs=h1_freq %>% filter(Var1 %in% tls_symbols) %>% filter(!is.na(final_degree)) %>% filter(!is.na(total_length))


#Doing permutation
{
ew=tfs[,c(8,12+unname(which(colSums(tfs[,c(13:ncol(tfs))])!=0)))]
bginfo=data.frame(sample_size=sort(unique(colSums(ew[,3:ncol(ew)]>0))),
                  sd=0,mean=0)
for(i in c(1:nrow(bginfo)))
{
  cat(paste0("\r",i," / ",nrow(bginfo)))
  perms=sapply(c(1:100000),function(x) mean(sample(ew$final_degree,size = bginfo$sample_size[i])))
  bginfo$sd[i]=sd(perms)
  bginfo$mean[i]=mean(perms)
}

eww=data.frame(size=colSums(ew[3:ncol(ew)]>0),obs=0,exp=0,sd=0)
for(i in c(1:nrow(eww)))
{
  print(i)
  eww$obs[i]=mean(ew$final_degree[ew[,rownames(eww)[i]]>0])
  j=which(bginfo$sample_size==eww$size[i])
  eww$exp[i]=bginfo$mean[j]
  eww$sd[i]=bginfo$sd[j]
  # eww$z[i]=(eww$obs[i]-eww$exp[i])/eww$sd[i]
}

err=eww %>% mutate(z=(obs-exp)/sd) %>% arrange(z)

}##luka##
# write.table(err,"features_enrichment_tf_init-lost_100kperms.tsv",sep="\t")
write.table(err,"features_enrichment_rbp_init-lost_100kperms.tsv",sep="\t")

er=err %>% filter(size>=5)
par(mar = c(6, 6, 4, 6),xpd=T)
plot(er$z,pch=19,col="whitesmoke",main="Features affecting degree of Nucleoplasm ",cex=er$size/1000,type="p")
text(c(1:nrow(er)),er$z,labels = rownames(er),cex=0.4)
View(er)
plot(er$z,er$obs-er$exp,pch=19,col="#EE44EE22",cex=er$size/300)
abline(h=c(-10,10),lty=2,col="red")
text(er$z,er$obs-er$exp,labels = rownames(er),cex=0.4)


# err=read.table("features_enrichment_tf_init-lost_100kperms.tsv")
ew=err %>% filter(abs(z)>2 & size>5)  %>% mutate(feature=rownames(.))
ew$feature <- factor(ew$feature, levels = ew$feature[order(ew$z)])

ggplot(ew, aes(x = feature, y = z, fill = z > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  coord_flip() +   # flip to put categories on y-axis
  labs(x = "", y = "Z score") +
  theme_minimal() +
  theme(legend.position = "none")

#comparing the degree of few GO terms
ddf<- bind_rows(
  h1_freq %>% filter(Var1 %in% tfs_symbols) %>% transmute(group = "TF", initial_degree,gained,lost),
  h1_freq %>% filter(Var1 %in% tls_symbols) %>% transmute(group = "TL", initial_degree,gained,lost),
  h1_freq %>% filter(Var1 %in% rbp_symbols) %>% transmute(group = "RBP", initial_degree,gained,lost)
)
# write.table(ddf,"TF_TL_RBP_degree.tsv",sep="\t")

ddf<- bind_rows(
  h1_freq %>% filter(Var1 %in% tfs_symbols) %>% filter(!(Var1 %in% drbp_symbols)) %>% transmute(group = "DBP",  initial_degree,gained,lost),
  h1_freq %>% filter(Var1 %in% rbp_symbols) %>% filter(!(Var1 %in% drbp_symbols)) %>% transmute(group = "RBP",  initial_degree,gained,lost),
  h1_freq %>% filter(Var1 %in% drbp_symbols) %>% transmute(group = "DRBP", initial_degree,gained,lost)
)
# write.table(ddf,"DBP_RBP_DRBP_degree.tsv",sep="\t")

ggplot(ddf %>% filter(group!="RB00P"), aes(x = group, y = initial_degree-lost)) + theme_classic()+
  geom_jitter(col="#EE442244",width = 0.2)+geom_violin(fill="transparent")+
  geom_boxplot(outlier.shape = NA,fill="transparent",width=0.1)
  # geom_hline(yintercept=ddf %>% filter(group=="TL") %>% pull( initial_degree-lost) %>% median(),lty=2,col="blue")

# stat_compare_means(method = "wilcox.test") 


#comparing GO terms of all GO terms
h1_freq=readRDS("h1_freq.Rdata") %>% mutate(corrected_degree=initial_degree-lost)
subcell=fread("subcellular_locations.tsv")
subcell=subcell %>% rowwise() %>% mutate(median_degree= median( h1_freq$corrected_degree[h1_freq$Var1 %in% (go2gene %>% filter(GO==`GO term`) %>% dplyr::select(SYMBOL) %>% unlist(use.names = F))],na.rm = T))



ddf=data.frame()
for(i in c(1:nrow(subcell)))
{
  print(i)
  temptab=h1_freq %>% filter(Var1 %in% (go2gene %>% filter(GO==subcell$`GO term`[i]) %>%
                               dplyr::select(SYMBOL) %>% unlist(use.names = F))) %>%
    transmute(group = subcell$`Subcellular location`[i], initial_degree,lost,gained,size=n())
ddf=rbind(ddf,temptab)
}


write.table(ddf,"subcellular_components_degree.tsv",sep="\t")


library(ggpubr)
vv=subcell %>% arrange(median_degree) %>% dplyr::select(`Subcellular location`) %>% unlist(use.names = F)
ggplot(ddf %>% 
         filter(size>=10) %>% 
         mutate(group = factor(group, levels = vv)),
       aes(x = group, y = initial_degree-lost)) + 
  theme_classic() +
  geom_jitter(col = "#EE442244", width = 0.2) +
  geom_violin(fill = "transparent") +
  geom_boxplot(outlier.shape = NA, fill = "transparent", width = 0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # geom_hline(yintercept=79,lty=2,col="blue")




# Comparing differnet groups of TFs
humantfs=fread("~/Downloads/DatabaseExtract_v_1.01.csv")##https://humantfs.ccbr.utoronto.ca/download.php
ddf<- bind_rows(
  h1_freq %>% filter(Var1 %in% (humantfs %>% filter(`Is TF?`=="Yes" & `Is C2H2 ZF(KRAB)?`==TRUE) %>% pull(`HGNC symbol`))) %>% transmute(group="TF with c2h2 zf (KRAB)",final_degree,num_author,initial_degree,change),
  h1_freq %>% filter(Var1 %in% (humantfs %>% filter(`Is TF?`=="Yes" & `Is C2H2 ZF(KRAB)?`==FALSE) %>% pull(`HGNC symbol`))) %>% transmute(group="TF without c2h2 zf(KRAB)",final_degree,num_author,initial_degree,change),
  # h1_freq %>% filter(h1_freq %>% dplyr::select(contains("C2H2")) %>% rowSums() >= 1) %>% filter(Var1 %in% (humantfs %>% filter(`Is TF?`=="Yes") %>% pull(`HGNC symbol`))) %>% transmute(group="TF with c2h2",final_degree,num_author,initial_degree,change)
   )
ddf<- bind_rows(
  h1_freq %>% filter(h1_freq %>% dplyr::select(starts_with("Zinc finger"),-contains("C2H2")) %>% rowSums() >= 1)  %>% transmute(group="Prot: ZF - c2h2 domain",final_degree,num_author,initial_degree,change),
 h1_freq %>% filter(h1_freq %>% dplyr::select(contains("C2H2")) %>% rowSums() >= 1)  %>% transmute(group="Prot: ZF+ c2h2 domain",final_degree,num_author,initial_degree,change),
 h1_freq %>% filter(h1_freq %>% dplyr::select(starts_with("Zinc finger"),-contains("C2H2")) %>% rowSums() >= 1) %>% filter(Var1 %in% (humantfs %>% filter(`Is TF?`=="Yes") %>% pull(`HGNC symbol`))) %>% transmute(group="TF: ZF - c2h2 domain",final_degree,num_author,initial_degree,change),
 h1_freq %>% filter(h1_freq %>% dplyr::select(contains("C2H2")) %>% rowSums() >= 1)  %>% filter(Var1 %in% (humantfs %>% filter(`Is TF?`=="Yes") %>% pull(`HGNC symbol`))) %>% transmute(group="TF: ZF+ c2h2 domain",final_degree,num_author,initial_degree,change)
)
ddf$initial_n=ddf$initial_degree/(ddf$num_author+1)
ddf$final_n=ddf$final_degree/(ddf$num_author+1)


plot1=function(plotwhat){
  # plotwhat="initial_degree"
  return(ggplot(ddf %>% filter(group!="RB00P"), aes(x = group, y = !!sym(plotwhat))) + theme_classic()+
  geom_jitter(col="#EE442211",width = 0.2,cex=3)+geom_violin(fill="transparent")+
  geom_boxplot(outlier.shape = NA,fill="transparent",width=0.1)+
  geom_hline(yintercept=ddf %>% group_by(group) %>% summarize(median=median(!!sym(plotwhat))) %>% pull(median),lty=2,col="blue")
       )
  }


p1= plot1("initial_degree")+coord_cartesian(ylim=c(0,200)) 
p2= plot1("final_degree")+coord_cartesian(ylim=c(0,200)) 
p3= plot1("num_author")+coord_cartesian(ylim=c(0,50))
p4= plot1("change")+coord_cartesian(ylim=c(-30,30)) + ylab("Change in interaction by our method")
(p1 | p2)/ (p3 | p4)

#tissue specificity of the tfs
temp1=fread("~/Downloads/rna_tissue_consensus.tsv")
temp1 <- temp1 %>% 
  group_by(Gene, `Gene name`) %>%
  summarise(max_tissue = Tissue[which.max(nTPM)],
            max_expr = max(nTPM),
            mean_expr = mean(nTPM),
            detected_tissues = sum(nTPM > 1),
            total_tissues = n(),
            .groups = "drop") %>%
  mutate(group = case_when(
    detected_tissues == 0 ~ "Not detected",
    max_expr / mean_expr >= 5 ~ "Tissue enriched/enhanced",
    detected_tissues > 0.75 * total_tissues ~ "Low tissue specificity",
    TRUE ~ "Intermediate"
  )) 
temp1= temp1 %>% dplyr::rename(Var1=`Gene name`) %>% left_join(. , h1_freq %>% dplyr::select(Var1,final_degree,initial_degree,num_author,change))
temp1=temp1 %>% filter(Var1 %in% (humantfs %>% filter(`Is TF?`=="Yes") %>% pull(`HGNC symbol`)))
# temp1=temp1 %>% filter(Var1 %in% (h1_freq %>% filter(`Zinc finger__C2H2-type 2`>=1) %>% pull(Var1)))
temp1 <- temp1 %>%
  mutate(group = factor(group, 
                        levels = c("Not detected", 
                                   "Intermediate", 
                                   "Tissue enriched/enhanced", 
                                   "Low tissue specificity")))

plot1=function(plotwhat,ylim){
return(ggplot(temp1 %>% filter(group!="RB00P"), aes(x = group, y = !!sym(plotwhat))) + theme_classic()+
  geom_jitter(col="#EE442244",width = 0.2)+geom_violin(fill="transparent")+
  geom_boxplot(outlier.shape = NA,fill="transparent",width=0.1)+
  geom_hline(yintercept=temp1 %>% group_by(group) %>% summarize(median=median(!!sym(plotwhat),na.rm = T)) %>% pull(median),lty=2,col="blue")+
  coord_cartesian(ylim = ylim)   # <-- zooms in on y-axis
)
}

(plot1("initial_degree",c(0,100)) | plot1("final_degree",c(0,100))) /
  (plot1("num_author",c(0,20)) | plot1("change",c(-15,15)))


}##luka##





