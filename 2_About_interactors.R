

matrix=readRDS("interaction_matrix.Rdata")
tem=readRDS("tem.Rdata")
ECdz_min=readRDS("ECdz_min.Rdata")

table(rownames(matrix)==rownames(ECdz_min))
threshold=2
hi=data.frame(Var1=rownames(matrix),
              initial_degree=rowSums(matrix>=1),
              final_degree=rowSums(ECdz_min>=threshold),
              gained=rowSums(matrix==0 & ECdz_min>=threshold),
              lost=rowSums(matrix>=1 & ECdz_min<threshold))
bias=readRDS("bias.Rdata")
bias=bias %>% mutate(total=A+B)
bias=left_join(bias,hi,by=c("Var1"))
bias=bias %>% mutate(change=gained-lost)
saveRDS(bias,"bias2.Rdata")



#############################################
#############################################
############### Whats the GO of the GAINED interaction?
##############################Doing GO enrichment
#############################################
{
#Curating GO terms
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(dplyr)
library(tidyr)
library(tibble)
  
go2gene <- AnnotationDbi::select(org.Hs.eg.db,keys = keys(org.Hs.eg.db, keytype = "GO"),columns = c("GO", "SYMBOL","ENTREZID"),keytype = "GO")
ew=go2gene %>% group_by(ONTOLOGY,SYMBOL) %>% summarize(GO=list(unique(GO)))

table(rownames(matrix)==rownames(ECdz_min))
gain=(matrix==0 & ECdz_min>=2)
lost=(matrix>0 & ECdz_min<2)

go= data.frame(SYMBOL=rownames(matrix))
go=left_join(go,ew)
go= go %>% group_by(SYMBOL) %>% mutate(gain=list(names(which(gain[SYMBOL[1],]==T))),
                                       lost=list(names(which(lost[SYMBOL[1],]==T))),
                                       orig=list(names(which(matrix[SYMBOL[1],]>=1))))  %>% ungroup()

# STEP 1: Collapse GO terms into GO_BP, GO_CC, GO_MF
go_collapsed <- go %>%
  group_by(SYMBOL) %>%
  summarise(
    GO_BP = list(GO[ONTOLOGY == "BP"]),
    GO_CC = list(GO[ONTOLOGY == "CC"]),
    GO_MF = list(GO[ONTOLOGY == "MF"]),
    gain = gain[1]
  ) %>% ungroup()

# STEP 2: Prepare universe
all_genes <- unique(unlist(go_collapsed$gain))
entrez_universe <- bitr(all_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
N_BP=go_collapsed %>% dplyr::select(GO_BP) %>% unlist() %>% unique() %>% length() #number of BP GO in universe
N_CC=go_collapsed %>% dplyr::select(GO_CC) %>% unlist() %>% unique() %>% length()
N_MF=go_collapsed %>% dplyr::select(GO_MF) %>% unlist() %>% unique() %>% length()
N_ALL= N_MF+N_CC+N_BP
  
# STEP 3/4: DOing GO enrichment for all gained interactors, and do hypergeometric test too at same time.
for (i in seq_len(nrow(go_collapsed))) {
    # for (i in c(7216:nrow(go_collapsed))) {
  #i=26
  sym <- go_collapsed$SYMBOL[i]
  genes <- go_collapsed$gain[[i]]
  go_bp <- unique(unlist(go_collapsed$GO_BP[[i]]))
  go_cc <- unique(unlist(go_collapsed$GO_CC[[i]]))
  go_mf <- unique(unlist(go_collapsed$GO_MF[[i]]))

  cat("Processing:", sym, "(", i, "/", nrow(go_collapsed), ")\n")
  
  if (length(genes) < 10) {
    final_df <- bind_rows(final_df, tibble(SYMBOL=sym,
                                           enr_BP=NA,enr_BP_vals=NA, n_BP=NA, k_BP=NA, K_BP=NA, p_BP=NA,
                                           enr_CC=NA,enr_CC_vals=NA, n_CC=NA, k_CC=NA, K_CC=NA, p_CC=NA,
                                           enr_MF=NA,enr_MF_vals=NA, n_MF=NA, k_MF=NA, K_MF=NA, p_MF=NA))
    next
  }
  
  # Try mapping SYMBOL to ENTREZID safely
  entrez_genes <- tryCatch({
    suppressMessages(
      bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )$ENTREZID
  }, error = function(e) {
    # If error occurs, return NULL instead of stopping
    NULL
  })
  
  # If mapping failed (entrez_genes is NULL) or too few genes
  if (is.null(entrez_genes) || length(entrez_genes) < 10) {
    final_df <- bind_rows(final_df, tibble(SYMBOL = sym,
                                           enr_BP = NA, enr_BP_vals = NA, n_BP = NA, k_BP = NA, K_BP = NA, p_BP = NA,
                                           enr_CC = NA, enr_CC_vals = NA, n_CC = NA, k_CC = NA, K_CC = NA, p_CC = NA,
                                           enr_MF = NA, enr_MF_vals = NA, n_MF = NA, k_MF = NA, K_MF = NA, p_MF = NA))
    next
  }
  enr <- tryCatch({
    enrichGO(gene = entrez_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "ALL", universe = entrez_universe,
             pAdjustMethod = "BH", qvalueCutoff = 0.05,pvalueCutoff = 1)
  }, error = function(e) NULL)
  
  if (is.null(enr) || nrow(enr) == 0) {
    final_df <- bind_rows(final_df, tibble(SYMBOL=sym,
                                           enr_BP=NA,enr_BP_vals=NA, n_BP=0, k_BP=0, K_BP=length(go_bp), p_BP=1,
                                           enr_CC=NA,enr_CC_vals=NA, n_CC=0, k_CC=0, K_CC=length(go_cc), p_CC=1,
                                           enr_MF=NA,enr_MF_vals=NA, n_MF=0, k_MF=0, K_MF=length(go_mf), p_MF=1))
    next
  }
  
  enr_df <- enr@result
  enr_df_BP=enr_df %>% filter(ONTOLOGY=="BP") 
  enr_df_MF=enr_df %>% filter(ONTOLOGY=="MF")
  enr_df_CC=enr_df %>% filter(ONTOLOGY=="CC")
  
  get_stats <- function(ontology, original_go) {
    enr_go <- enr_df %>% filter(ONTOLOGY == ontology) %>% pull(ID) %>% unique()
    N <- eval(parse(text=paste0("N_",ontology)))  # universe size in this run
    K <- length(original_go)   # original GO terms of the symbol
    n <- length(enr_go)                # enriched GO terms
    k <- length(intersect(enr_go, original_go))  # overlap
    p <- if (K > 0 && N > 0 && n > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA
    return(c(n, k, K, p))
  }
  
  vals_bp <- get_stats("BP", go_bp)
  vals_cc <- get_stats("CC", go_cc)
  vals_mf <- get_stats("MF", go_mf)
  
  final_df <- bind_rows(final_df, tibble(
    SYMBOL = sym,
    enr_BP=list(enr_df_BP$ID),enr_BP_vals=list(enr_df_BP$qvalue),n_BP = vals_bp[1], k_BP = vals_bp[2], K_BP = vals_bp[3], p_BP = vals_bp[4],
    enr_CC=list(enr_df_CC$ID),enr_CC_vals=list(enr_df_CC$qvalue),n_CC = vals_cc[1], k_CC = vals_cc[2], K_CC = vals_cc[3], p_CC = vals_cc[4],
    enr_MF=list(enr_df_MF$ID),enr_MF_vals=list(enr_df_MF$qvalue),n_MF = vals_mf[1], k_MF = vals_mf[2], K_MF = vals_mf[3], p_MF = vals_mf[4]
  )
  )
}
goo=left_join(go_collapsed,final_df)
goo$gain_n= goo %>% rowwise() %>% transmute(length(unlist(gain))) %>% unlist(use.names = F)
goo=left_join(goo, bias %>% transmute(SYMBOL=Var1,num_author,initial_degree,initial_interaction=total))
saveRDS(goo,"goo.Rdata")


#STEP 5: Finding hypergeometric test on ALL GO terms to find p_ALL
goo=goo %>% rowwise() %>%  mutate(n_ALL= ifelse(gain_n>=10,1,NA),
                                  k_ALL= ifelse(gain_n>=10,1,NA),
                                  K_ALL= ifelse(gain_n>=10,1,NA),
                                  p_ALL= ifelse(gain_n>=10,1,NA))
for(i in c(1:nrow(goo)))
{
  #i=1
  cat(paste0("\r ",i))
  if(goo$gain_n[i]<10) next
  N <- N_ALL # universe size in this run
  original_go=c(unlist(goo$GO_BP[i]),unlist(goo$GO_CC[i]),unlist(goo$GO_MF[i]))
  enr_go= unique(unlist(goo[i,c("enr_BP","enr_CC","enr_MF")]))
  K <- length(original_go)  # original GO terms of the symbol
  n <- length(enr_go)                # enriched GO terms
  k <- length(intersect(enr_go, original_go))  # overlap
  p <- if (K > 0 && N > 0 && n > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA
  goo$n_ALL[i]= n
  goo$k_ALL[i]= k
  goo$K_ALL[i]= K
  goo$p_ALL[i]= p
}
saveRDS(goo,"goo.Rdata")




# STEP 6: Analyzing hypergeometric p value patterns
goo=readRDS("goo.Rdata")
{plot(goo$gain_n,-log10(goo$p_ALL))
plot(goo$num_author,-log10(goo$p_ALL))

table(goo$p_BP<=0.05 | goo$p_CC<=0.05 | goo$p_MF<=0.05)
table(goo$p_BP<=0.05 & goo$p_CC<=0.05 & goo$p_MF<=0.05)
table(goo$p_ALL<=0.05)

dev.off()
hya=goo %>% filter(K_ALL>0 ) %>% filter(!is.na(p_ALL)) #some p values are NA, because the small universe is empty
plot(log2(hya$K_ALL+1),-log10(hya$p_ALL),
     xlab="log2(number of GO terms related to  poi)",ylab="how significant the overlap?")
plot(sort(-log10(goo$p_ALL)))

#plotting boxplots
w <- 10;m <- 50
brk <- c(seq(0, m, by = w), Inf)
lbl <- c(paste(seq(1, m - w + 1, by = w), seq(w, m, by = w), sep = "-"), paste0(m + 1, "+"))
hya=goo %>% filter(K_ALL>0) %>% filter(!is.na(p_ALL)) %>% filter(gain_n>10) %>%
  mutate(g = cut(K_ALL, breaks = brk, labels = lbl, include.lowest = TRUE, right = TRUE),
         lp = -log10(p_ALL))

ggplot(hya, aes(x = g, y = lp)) +
  geom_jitter(color = "#11FF1155") +
  geom_boxplot(fill = "transparent") +
  geom_hline(yintercept = 1.3, color = "red", linetype = "dashed") +
  xlab("Number of GO terms of a protein") +
  ylab("Significance of overlap between GO terms of a protein and it's new interactors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  print("hi")



hya=goo %>% filter(K_CC!=0) %>% filter(!is.na(p_CC)) #some p values are NA, because the small universe is empty
table(hya$p_CC<=0.05)
plot(hya$gain_n,-log10(hya$p_CC))


hya=goo %>% filter(K_MF!=0) %>% filter(!is.na(p_MF)) #some p values are NA, because the small universe is empty
table(hya$p_MF<=0.05)
plot(hya$gain_n,-log10(hya$p_MF))

hya=-log10(goo %>% dplyr::select(starts_with("p_")))
rownames(hya)=goo$SYMBOL
plot(hya)
plot(hya$p_BP,hya$p_CC,pch=19,col="whitesmoke");text(hya$p_BP,hya$p_CC,labels = final_df$SYMBOL,cex=0.4)
plot(hya$p_BP,hya$p_MF,pch=19,col="whitesmoke");text(hya$p_BP,hya$p_MF,labels = final_df$SYMBOL,cex=0.4)
plot(hya$p_MF,hya$p_CC,pch=19,col="whitesmoke");text(hya$p_MF,hya$p_CC,labels = final_df$SYMBOL,cex=0.4)

hya2=plotTriangle(hya,plotAll = T,col=c("#22EE22AA","#AA222222","#1199FF88"),
                  maxvaluerange = c(2,Inf)) %>% filter(color!="whitesmoke")
text(hya2$comx,hya2$comy,labels=rownames(hya2),cex=0.4,)

}##luka##



##Using semantic similairty
{
library(GOSemSim)
go1= unlist(hya$GO_BP[12])
go2= unlist(hya$enr_BP[12])
semdata_BP <- godata('org.Hs.eg.db', ont = "BP")
semdata_CC <- godata('org.Hs.eg.db', ont = "CC")
semdata_MF <- godata('org.Hs.eg.db', ont = "MF")

simil <- mgoSim(go1, go2, semData = semdata, measure="Wang", combine="BMA")

goo2=goo %>% rowwise() %>% mutate(semsim_BP=mgoSim(unlist(GO_BP), unlist(enr_BP), semData = semdata, measure="Wang", combine="BMA"))
goo2=goo2 %>% rowwise() %>% mutate(semsim_CC=mgoSim(unlist(GO_CC), unlist(enr_CC), semData = semdata, measure="Wang", combine="BMA"))
goo2=goo %>% rowwise() %>% mutate(semsim_BP=mgoSim(unlist(GO_BP), unlist(enr_BP), semData = semdata, measure="Wang", combine="BMA"))

goo2$semsim_BP=NA;goo2$semsim_CC=NA;goo2$semsim_MF=NA

for(i in c(1:nrow(goo2)))
{
  if(goo2$gain_n[i]>=10)
  mgoSim(unlist(goo$GO_BP[1]), unlist(goo$enr_BP[1]), semData = semdata, measure="Wang", combine="BMA")
}

#plotting boxplots
w <- 10;m <- 50
brk <- c(seq(0, m, by = w), Inf)
lbl <- c(paste(seq(1, m - w + 1, by = w), seq(w, m, by = w), sep = "-"), paste0(m + 1, "+"))
hya=goo2 %>% filter(K_BP>0) %>% filter(!is.na(semsim_BP)) %>% filter(gain_n>=10) %>%
  mutate(g = cut(K_BP, breaks = brk, labels = lbl, include.lowest = TRUE, right = TRUE),
         lp = semsim_BP)
# Count per group
counts <- hya %>% group_by(g) %>% summarise(n = n(), .groups = "drop")

# Plot
ggplot(hya, aes(x = g, y = lp)) +
  geom_jitter(color = "#11FF1155") +
  geom_boxplot(fill = "transparent") +
  geom_hline(yintercept = 0.5, color = "blue", linetype = "dotted") +
  geom_text(data = counts, aes(x = g, y = 1.05, label = n), inherit.aes = FALSE, size = 3.5) +
  xlab("Number of GO terms of a protein") +
  ylab("GO semantic similarity between GO terms of a protein and its new interactors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mean(hya$semsim_BP>=0.5)


}##luka##



##Doing overall analysis on hypergeometric p value
{
  goo4=goo %>% filter(gain_n>=10) %>%  rowwise() %>% filter(length(enr_BP_vals)>0 & length(GO_BP)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  goo4=goo4[!sapply(goo4$GO_BP, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  table(goo4$p_BP<=0.05)
  obs=sum(goo4$p_BP<=0.05)
perms=sapply(c(1:1000), function(x)
  {goo5=goo4
  idx = sample(nrow(goo5))  # Random permutation of row indices
  goo5$GO_BP <- goo5$GO_BP[idx]
  goo5$K_BP  <- goo5$K_BP[idx]
 goo5 %>% rowwise() %>% transmute(p_suffl= phyper(sum(unlist(enr_BP) %in% unlist(GO_BP)) -1,K_BP,N_BP - K_BP,
                                                     n_BP,lower.tail = F) <= 0.05) %>% unlist(use.names = F) %>% sum()}
)

hist(perms,breaks = 10000,xlim=c(0,2000),main="Number of proteins that have overlap significance <= 0.05")
abline(v=obs,lty=2)

  }

}##luka##


#############################################
#############################################
############### Whats the GO of the LOST interaction? 
###############################Doing GO enrichment
#############################################
{
  #Curating GO terms
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GO.db)
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  go2gene <- AnnotationDbi::select(org.Hs.eg.db,keys = keys(org.Hs.eg.db, keytype = "GO"),columns = c("GO", "SYMBOL","ENTREZID"),keytype = "GO")
  ew=go2gene %>% group_by(ONTOLOGY,SYMBOL) %>% summarize(GO=list(unique(GO)))
  
  table(rownames(matrix)==rownames(ECdz_min))
  gain=(matrix==0 & ECdz_min>=2)
  lost=(matrix>0 & ECdz_min<2)
  
  go= data.frame(SYMBOL=rownames(matrix))
  go=left_join(go,ew)
  go= go %>% group_by(SYMBOL) %>% mutate(gain=list(names(which(gain[SYMBOL[1],]==T))),
                                         lost=list(names(which(lost[SYMBOL[1],]==T))),
                                         orig=list(names(which(matrix[SYMBOL[1],]>=1))))  %>% ungroup()
  
  # STEP 1: Collapse GO terms into GO_BP, GO_CC, GO_MF
  go_collapsed <- go %>%
    group_by(SYMBOL) %>%
    summarise(
      GO_BP = list(GO[ONTOLOGY == "BP"]),
      GO_CC = list(GO[ONTOLOGY == "CC"]),
      GO_MF = list(GO[ONTOLOGY == "MF"]),
      lost = lost[1]
    ) %>% ungroup()
  
  # STEP 2: Prepare universe
  all_genes <- unique(unlist(go_collapsed$lost))
  entrez_universe <- bitr(all_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  N_BP=go_collapsed %>% dplyr::select(GO_BP) %>% unlist() %>% unique() %>% length() #number of BP GO in universe
  N_CC=go_collapsed %>% dplyr::select(GO_CC) %>% unlist() %>% unique() %>% length()
  N_MF=go_collapsed %>% dplyr::select(GO_MF) %>% unlist() %>% unique() %>% length()
  N_ALL= N_MF+N_CC+N_BP
  
  # STEP 3/4: DOing enrichment for all lost interactors, and do hypergeometric test too at same time.
  for (i in seq_len(nrow(go_collapsed))) {
    # for (i in c(7216:nrow(go_collapsed))) {
    #i=26
    sym <- go_collapsed$SYMBOL[i]
    genes <- go_collapsed$lost[[i]]
    go_bp <- unique(unlist(go_collapsed$GO_BP[[i]]))
    go_cc <- unique(unlist(go_collapsed$GO_CC[[i]]))
    go_mf <- unique(unlist(go_collapsed$GO_MF[[i]]))
    
    cat("Processing:", sym, "(", i, "/", nrow(go_collapsed), ")\n")
    
    if (length(genes) < 10) {
      final_df <- bind_rows(final_df, tibble(SYMBOL=sym,
                                             enr_BP=NA,enr_BP_vals=NA, n_BP=NA, k_BP=NA, K_BP=NA, p_BP=NA,
                                             enr_CC=NA,enr_CC_vals=NA, n_CC=NA, k_CC=NA, K_CC=NA, p_CC=NA,
                                             enr_MF=NA,enr_MF_vals=NA, n_MF=NA, k_MF=NA, K_MF=NA, p_MF=NA))
      next
    }
    
    # Try mapping SYMBOL to ENTREZID safely
    entrez_genes <- tryCatch({
      suppressMessages(
        bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      )$ENTREZID
    }, error = function(e) {
      # If error occurs, return NULL instead of stopping
      NULL
    })
    
    # If mapping failed (entrez_genes is NULL) or too few genes
    if (is.null(entrez_genes) || length(entrez_genes) < 10) {
      final_df <- bind_rows(final_df, tibble(SYMBOL = sym,
                                             enr_BP = NA, enr_BP_vals = NA, n_BP = NA, k_BP = NA, K_BP = NA, p_BP = NA,
                                             enr_CC = NA, enr_CC_vals = NA, n_CC = NA, k_CC = NA, K_CC = NA, p_CC = NA,
                                             enr_MF = NA, enr_MF_vals = NA, n_MF = NA, k_MF = NA, K_MF = NA, p_MF = NA))
      next
    }
    enr <- tryCatch({
      enrichGO(gene = entrez_genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
               ont = "ALL", universe = entrez_universe,
               pAdjustMethod = "BH", qvalueCutoff = 0.05,pvalueCutoff = 1)
    }, error = function(e) NULL)
    
    if (is.null(enr) || nrow(enr) == 0) {
      final_df <- bind_rows(final_df, tibble(SYMBOL=sym,
                                             enr_BP=NA,enr_BP_vals=NA, n_BP=0, k_BP=0, K_BP=length(go_bp), p_BP=1,
                                             enr_CC=NA,enr_CC_vals=NA, n_CC=0, k_CC=0, K_CC=length(go_cc), p_CC=1,
                                             enr_MF=NA,enr_MF_vals=NA, n_MF=0, k_MF=0, K_MF=length(go_mf), p_MF=1))
      next
    }
    
    enr_df <- enr@result
    enr_df_BP=enr_df %>% filter(ONTOLOGY=="BP") 
    enr_df_MF=enr_df %>% filter(ONTOLOGY=="MF")
    enr_df_CC=enr_df %>% filter(ONTOLOGY=="CC")
    
    get_stats <- function(ontology, original_go) {
      enr_go <- enr_df %>% filter(ONTOLOGY == ontology) %>% pull(ID) %>% unique()
      N <- eval(parse(text=paste0("N_",ontology)))  # universe size in this run
      K <- length(original_go)   # original GO terms of the symbol
      n <- length(enr_go)                # enriched GO terms
      k <- length(intersect(enr_go, original_go))  # overlap
      p <- if (K > 0 && N > 0 && n > 0) phyper(k - 1, K, N - K, n, lower.tail = FALSE) else NA
      return(c(n, k, K, p))
    }
    
    vals_bp <- get_stats("BP", go_bp)
    vals_cc <- get_stats("CC", go_cc)
    vals_mf <- get_stats("MF", go_mf)
    
    final_df <- bind_rows(final_df, tibble(
      SYMBOL = sym,
      enr_BP=list(enr_df_BP$ID),enr_BP_vals=list(enr_df_BP$qvalue),n_BP = vals_bp[1], k_BP = vals_bp[2], K_BP = vals_bp[3], p_BP = vals_bp[4],
      enr_CC=list(enr_df_CC$ID),enr_CC_vals=list(enr_df_CC$qvalue),n_CC = vals_cc[1], k_CC = vals_cc[2], K_CC = vals_cc[3], p_CC = vals_cc[4],
      enr_MF=list(enr_df_MF$ID),enr_MF_vals=list(enr_df_MF$qvalue),n_MF = vals_mf[1], k_MF = vals_mf[2], K_MF = vals_mf[3], p_MF = vals_mf[4]
    )
    )
  }
  loo=left_join(go_collapsed,final_df)
  loo$lost_n= goo %>% rowwise() %>% transmute(length(unlist(lost))) %>% unlist(use.names = F)
  loo=left_join(loo, bias %>% transmute(SYMBOL=Var1,num_author,initial_degree,initial_interaction=total))
  saveRDS(loo,"loo.Rdata")
  



}##luka##



#############################################
#############################################
############### Comparing GO of gained and lost interactors with GO of original protein
#################### by random shuffling (permutation test)
#############################################
#############################################
{
  
  goo=readRDS("goo.Rdata")
  loo=readRDS("loo.Rdata")
  
  #Function for shuffling permutation
  permit=function(df,posn=1) {
    colnames(df)=c("GO1","GO2")
    df$GO1=sample(df$GO1)
    perm=df %>% rowwise() %>% transmute(isitthere=unlist(GO2)[posn] %in% unlist(GO1)) %>% unlist(use.names = F) %>% sum()
    return(perm)
  }
  
  #Shuffle analysis for gained interactors
  {
  sflzg=data.frame(posn=c(1:5),z_BP=0,z_CC=0,z_MF=0)
  #####Doing overall analysis on BP
  hya=goo %>% filter(gain_n>=10) %>%  rowwise() %>% filter(length(enr_BP_vals)>0 & length(GO_BP)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  hya=hya[!sapply(hya$GO_BP, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  for(i in c(1:nrow(sflzg)))
  { print(i)
    perms=sapply(c(1:10000),function(x) permit(hya %>% dplyr::select(GO_BP,enr_BP),i))
    obs=hya %>% rowwise() %>% transmute(isitthere=unlist(enr_BP)[i] %in% unlist(GO_BP)) %>% dplyr::select(isitthere) %>% unlist(use.names = F) %>% sum()
    sflzg$z_BP[i]=(obs-mean(perms))/sd(perms)
  }
  #####Doing overall analysis on CC
  hya=goo %>% filter(gain_n>=10) %>%  rowwise() %>% filter(length(enr_CC_vals)>0 & length(GO_CC)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  hya=hya[!sapply(hya$GO_CC, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  for(i in c(1:nrow(sflzg)))
  { print(i)
    perms=sapply(c(1:10000),function(x) permit(hya %>% dplyr::select(GO_CC,enr_CC),i))
    obs=hya %>% rowwise() %>% transmute(isitthere=unlist(enr_CC)[i] %in% unlist(GO_CC)) %>% dplyr::select(isitthere) %>% unlist(use.names = F) %>% sum()
    sflzg$z_CC[i]=(obs-mean(perms))/sd(perms)
  }
  #####Doing overall analysis on MF
  hya=goo %>% filter(gain_n>=10) %>%  rowwise() %>% filter(length(enr_MF_vals)>0 & length(GO_MF)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  hya=hya[!sapply(hya$GO_MF, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  for(i in c(1:nrow(sflzg)))
  { print(i)
    perms=sapply(c(1:10000),function(x) permit(hya %>% dplyr::select(GO_MF,enr_MF),i))
    obs=hya %>% rowwise() %>% transmute(isitthere=unlist(enr_MF)[i] %in% unlist(GO_MF)) %>% dplyr::select(isitthere) %>% unlist(use.names = F) %>% sum()
    sflzg$z_MF[i]=(obs-mean(perms))/sd(perms)
  }
  
  write.table(sflzg,"sflzg.tsv",sep = "\t")  
}##luka##
  
  #Shuffle analysis for lost interactors
  { 
  sflzl=data.frame(posn=c(1:5),z_BP=0,z_CC=0,z_MF=0)
  #####Doing overall analysis on BP
  hya=loo %>% filter(lost_n>=10) %>%  rowwise() %>% filter(length(enr_BP_vals)>0 & length(GO_BP)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  hya=hya[!sapply(hya$GO_BP, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  for(i in c(1:nrow(sflzl)))
  { print(i)
    perms=sapply(c(1:10000),function(x) permit(hya %>% dplyr::select(GO_BP,enr_BP),i))
    obs=hya %>% rowwise() %>% transmute(isitthere=unlist(enr_BP)[i] %in% unlist(GO_BP)) %>% dplyr::select(isitthere) %>% unlist(use.names = F) %>% sum()
    sflzl$z_BP[i]=(obs-mean(perms))/sd(perms)
  }
  #####Doing overall analysis on CC
  hya=loo %>% filter(lost_n>=10) %>%  rowwise() %>% filter(length(enr_CC_vals)>0 & length(GO_CC)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  hya=hya[!sapply(hya$GO_CC, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  for(i in c(1:nrow(sflzl)))
  { print(i)
    perms=sapply(c(1:10000),function(x) permit(hya %>% dplyr::select(GO_CC,enr_CC),i))
    obs=hya %>% rowwise() %>% transmute(isitthere=unlist(enr_CC)[i] %in% unlist(GO_CC)) %>% dplyr::select(isitthere) %>% unlist(use.names = F) %>% sum()
    sflzl$z_CC[i]=(obs-mean(perms))/sd(perms)
  }
  #####Doing overall analysis on MF
  hya=loo %>% filter(lost_n>=10) %>%  rowwise() %>% filter(length(enr_MF_vals)>0 & length(GO_MF)>0) #removing those entries where both the proteins and interactors have atleast one GO terms.
  hya=hya[!sapply(hya$GO_MF, function(x) is.null(unlist(x))),] #Removing those proteins which didnot have any GO terms.
  for(i in c(1:nrow(sflzl)))
  { print(i)
    perms=sapply(c(1:10000),function(x) permit(hya %>% dplyr::select(GO_MF,enr_MF),i))
    obs=hya %>% rowwise() %>% transmute(isitthere=unlist(enr_MF)[i] %in% unlist(GO_MF)) %>% dplyr::select(isitthere) %>% unlist(use.names = F) %>% sum()
    sflzl$z_MF[i]=(obs-mean(perms))/sd(perms)
  }
  write.table(sflzl,"sflzl.tsv",sep = "\t")  
}##luka##
  
  #######Plotting Bar plots
  {
    
    sflzg=read.table("sflzg.tsv")
    sflzl=read.table("sflzl.tsv")  
  par(mfrow=c(2,1))
  
  z_matrix <- t(rbind(sflzg$z_BP, sflzg$z_CC, sflzg$z_MF))
  colnames(z_matrix) <- c("BP", "CC", "MF")
  rownames(z_matrix) <- sflzg$posn
  
  barplot(t(z_matrix), horiz=TRUE, beside=TRUE, las=1, col=c("skyblue", "lightgreen", "salmon"), xlab="z", main="Z Scores per position (Gained Interactors)", xlim=c(0,200))
  legend("topright", legend=c("BP", "CC", "MF"), fill=c("skyblue", "lightgreen", "salmon"), title="Category")

  z_matrix <- t(rbind(sflzl$z_BP, sflzl$z_CC, sflzl$z_MF))
  colnames(z_matrix) <- c("BP", "CC", "MF")
  rownames(z_matrix) <- sflzl$posn
  
  barplot(t(z_matrix), horiz=TRUE, beside=TRUE, las=1, col=c("skyblue", "lightgreen", "salmon"), xlab="z", main="Z Scores per position (Lost Interactors)", xlim=c(0,200))
  legend("topright", legend=c("BP", "CC", "MF"), fill=c("skyblue", "lightgreen", "salmon"), title="Category")
  
  }##luka##
  
  
  }##luka##


