

#First crucial step
detach("package:MASS", unload = TRUE) ##Should do  this because there was conflict with its select function while plotting cross-polymerase interaction plots below



##Loading some necesary files
matrix=readRDS("interaction_matrix.Rdata")
ECdz_min=readRDS("ECdz_min.Rdata")
########################################################################
########################################################################
################################Interaction between three rna polymerase
########################################################################
########################################################################
{
col_new="#FB8451"
col_old="#487A9C"

library(visNetwork)
library(dplyr)

unique_pol_subs=rownames(matrix)[(grepl("POLR",rownames(matrix)) & !rownames(matrix) %in% c("POLR1C","POLR1D","POLR2E","POLR2F","POLR2H","POLR2K","POLR2L","POLRMT")) | rownames(matrix) %in% c("TWISTNB","CD3EAP","ZNRD1","GRINL1A","BN51T","CRCP") ]
unique_pol_subs

temp1=matrix[unique_pol_subs,unique_pol_subs]  #actual matrix
temp2=((temp1==0)* ECdz_min[unique_pol_subs,unique_pol_subs])>=2 #Gained interactions
temp3= temp1>0 & ECdz_min[unique_pol_subs,unique_pol_subs]>=2 #Old interactions

##### Choose one of the two below  lines for old or gained interaction
if(F)
{
temp2=temp2 & F # Run this line to remove gained interactions in the plot
temp3=temp3 & F # Run this line to remove popularity corrected interactions.
}

# temp3=ECdz_min<2 # No interactions

#Remving isoforms of POLR2J
iso <- c("POLR2J", "POLR2J2", "POLR2J3", "POLR2J4")
temp2[, "POLR2J"] <- apply(temp2[, iso], 1, any);temp2["POLR2J", ] <- apply(temp2[iso, ], 2, any)
temp2 <- temp2[!rownames(temp2) %in% iso[-1], !colnames(temp2) %in% iso[-1]]
temp3[, "POLR2J"] <- apply(temp3[, iso], 1, any);temp3["POLR2J", ] <- apply(temp3[iso, ], 2, any)
temp3 <- temp3[!rownames(temp3) %in% iso[-1], !colnames(temp3) %in% iso[-1]]
unique_pol_subs=rownames(temp3)
print(ifelse(sum(temp3[,"POLR2J"] & temp2[,"POLR2J"])==0,"Great !!Everything all right","WAIT !!!!!!! Some interactions are termed as both OLD and GAINED Interactoin, !!Please fix it before moving ahead"))

pol1_subs= c("POLR1A","POLR1B","POLR1E","CD3EAP","TWISTNB","ZNRD1")
pol1_subs_abb=c("1A","1B","1E","1G","1F","1H")
pol2_subs= c("POLR2A","POLR2B","POLR2C","POLR2D","POLR2G","POLR2I","POLR2J","POLR2M") # "POLR2J2","POLR2J3","POLR2J4" are clubbed with POLR2J
pol2_subs_abb= c("2A","2B","2C","2D","2G","2I","2J","2M")
pol3_subs= c("POLR3A","POLR3B","POLR3C","POLR3D","POLR3E","POLR3F","POLR3G","POLR3GL","POLR3H","POLR3K","CRCP")
pol3_subs_abb= c("3A","3B","3C","3D","3E","3F","3G","3GL","3H","3K","3I")


# quartz(height = 4,width = 4)
#plot it better
####for high quality::Expand the Viewer tab then click Export>Save as Image, and zoom to fit the tab screen and then save
{
  # --- 1. Build node list ---
  nodes = bind_rows(
    data.frame(id = pol1_subs, group = "POL1", label = pol1_subs_abb),
    data.frame(id = pol2_subs, group = "POL2", label = pol2_subs_abb),
    data.frame(id = pol3_subs, group = "POL3", label = pol3_subs_abb)
  )
  pos = list(
    POL1 = c(x =   0, y =  -300),
    POL2 = c(x = -300, y = 150),
    POL3 = c(x =  300, y = 150)
  )
  
  set.seed(42)
  nodes = nodes %>%
    group_by(group) %>%
    mutate(
      n     = n(),
      idx   = row_number(),
      angle = 2 * pi * idx / n,
      x     = pos[[group[1]]]["x"] + 100 * cos(angle),
      y     = pos[[group[1]]]["y"] + 100 * sin(angle)
    ) %>%
    ungroup() %>%
    mutate(
      color   = case_when(
        group == "POL1" ~ "#B5503A",
        group == "POL2" ~ "#C9A641",
        group == "POL3" ~ "#2e6787"
      ),
      physics = FALSE,
      font.size = 25,
      font.vadjust = -50,
      font.color="white"
    )
  
  # --- 2. Build edge list ---
  # alpha and width scale with ECdz_min value (capped at 10 for visual sanity)
  matrix_to_edges = function(mask_mat, ecz_mat, group_map, rgb_base) {
    as.data.frame(as.table(mask_mat)) %>%
      rename(from = Var1, to = Var2, weight = Freq) %>%
      filter(weight > 0) %>%
      filter(as.character(from) != as.character(to)) %>%
      mutate(
        from = as.character(from),
        to   = as.character(to)
      ) %>%
      left_join(
        as.data.frame(as.table(ecz_mat)) %>%
          rename(from = Var1, to = Var2, ecz = Freq) %>%
          mutate(from = as.character(from), to = as.character(to)),
        by = c("from", "to")
      ) %>%
      mutate(
        grp_from = group_map[from],
        grp_to   = group_map[to]
      ) %>%
      filter(!is.na(grp_from) & !is.na(grp_to)) %>%
      filter(grp_from != grp_to) %>%
      filter(from < to) %>%
      mutate(
        ecz_scaled = pmin(ecz, 10),                          # cap at 10
        width      = 1,#1 + (ecz_scaled / 10) * 9,             # width: 1–10
        alpha      = 1,#0.3 + (ecz_scaled / 10) * 0.65,        # alpha: 0.3–0.95
        color      = sprintf("rgba(%s,%.2f)", rgb_base, alpha)
      ) %>%
      select(from, to, width, color)
  }
  
  group_map  = setNames(nodes$group, nodes$id)
  ecz_sub    = ECdz_min[unique_pol_subs, unique_pol_subs]
  old_mat    = temp3[unique_pol_subs, unique_pol_subs]
  new_mat    = temp2[unique_pol_subs, unique_pol_subs]
  
  edges_old = matrix_to_edges(old_mat, ecz_sub, group_map, rgb_base = paste0(col2rgb(col_old)[,1],collapse=","))   # blue
  edges_new = matrix_to_edges(new_mat, ecz_sub, group_map, rgb_base = paste0(col2rgb(col_new)[,1],collapse=","))   # orange
  
  edges = bind_rows(edges_old, edges_new)
  
  # --- 3. Plot ---
  visNetwork(nodes, edges, width = "100%", height = "1200px") %>%
    visOptions(highlightNearest = TRUE) %>%
    visInteraction(hover = TRUE) %>%
    visPhysics(enabled = FALSE)
  
}##luka##



}##Luka##




############################################################
############################################################
##############################Analyzing Mass spec result of Sihang
############################################################
############################################################
{
  
  ew=ECdz_min[,c("ATRN","FBRSL1","TTC13")] %>% as.data.frame()
  require(data.table)
  atr=fread("Sihang_maspec_data/ATRN.hit.txt",sep="\t",drop = 1) %>% mutate(pulled=name %in% rownames(ew)[ew$ATRN>=2])
  fbr=fread("Sihang_maspec_data/FBRSL1.hit.txt",sep="\t",drop = 1) %>% mutate(pulled=name %in% rownames(ew)[ew$FBRSL1>=2])
  ttc=fread("Sihang_maspec_data/TTC31.hit.txt",sep="\t",drop = 1) %>% mutate(pulled=name %in% rownames(ew)[ew$TTC31>=2]) # FIle name is mispelled as TTC31
  
  
  masp=fread("Sihang_maspec_data/25-170 Sihang_Mascot_AverQuant.txt",sep="\t")
  masp=masp %>% rowwise() %>% mutate(name=strsplit(strsplit(`Fasta headers`,"GN=")[[1]][2]," ")[[1]][1]) %>%
    ungroup() %>% filter(!is.na(name)) %>%
    # mutate(ATRN=`Intensity [ATRN]`/`Intensity [GFP]`,
    #                      FBRSL1=`Intensity [FBRSL1]`/`Intensity [GFP]`,
    #                      TTC31=`Intensity [TTC31]`/`Intensity [GFP]`)
    mutate(ATRN=`Peptide XICs [ATRN]`/`Peptide XICs [GFP]`,
           FBRSL1=`Peptide XICs [FBRSL1]`/`Peptide XICs [GFP]`,
           TTC13=`Peptide XICs [TTC31]`/`Peptide XICs [GFP]`)
  
  library(VennDiagram)
  gain_matrix2=((ECdz_min>=2) & matrix==0)
  lost_matrix2=((ECdz_min<2) & matrix>0)
  same_matrix2=((ECdz_min>=2) & matrix>0)
  
  par(mfcol=c(3,2))
  prot="ATRN"
  prot="FBRSL1"
  prot="TTC13"
  
  vennn <- venn.diagram(
    x = list(
      Gained    = names(which(gain_matrix2[,prot])),
      Old    = names(which(matrix[,prot]>0)),
      # Lost = names(which(lost_matrix2[,prot])),
      # Pulled = fread(paste0("~/Downloads/",prot,".hit.txt"),sep="\t",drop = 1) %>% pull(name)),
      Pulled = masp$name[which(masp[prot]>=Inf)]),
    filename = NULL,        # keep in R, instead of writing to a file
    fill = c("darkorange", "skyblue", "grey"),alpha=1)
  
  grid::grid.newpage()
  grid::grid.draw(vennn)
  
  
  #permutations 
  { 
    obs=masp$name[which(masp[prot]>=Inf)] #selecting the pulled down proteins
    obs=obs[obs %in% colnames(ECdz_min)] # only selecting those exist in our ppi network
    obs=obs[obs %in% names(which(matrix[,prot]==0))]
    length(obs)
    perms=sapply(c(1:10000),function(x) median(sample(ECdz_min[prot,],size = length(obs)),na.rm = T))
    plot(density(perms),xlim=c(0,0.2),main=prot)
    abline(v=median(ECdz_min[prot,obs],na.rm = T),lty=2)
    (median(ECdz_min[prot,obs],na.rm = T)-mean(perms))/sd(perms)
  }
  
  
  temp_fn1= function(prot)
  {  
    obs=masp$name[which(masp[prot]>=Inf)] #selecting the pulled down proteins
    obs=obs[obs %in% names(which(matrix[,prot]==0))] # only selecting the ones that are not known interactor of our poi in biogrid
    # obs=obs[obs %in% names(which(matrix[,prot]==0))]
    perms=sapply(c(1:100000),function(x) median(sample(ECdz_min[prot,],size = length(obs)),na.rm = T))
    obs=median(ECdz_min[prot,obs],na.rm = T)
    return(list(permutations=perms,observed=obs))
  }
  
  
  masspec_sihang=list(ATRN= temp_fn1("ATRN"),
                      TTC13= temp_fn1("TTC13"),
                      FBRSL1= temp_fn1("FBRSL1"))
  saveRDS(masspec_sihang,"masspec_sihang.Rdata")
  
  
}



########################################################################
########################################################################
################################Analyzing Mass spec result of Jenny
########################################################################
########################################################################
{
  
  ew=fread("./MS_Result_enrich_15oct_noMaf1.csv")
  prot="POLR3G"
  prot="POLR3GL"
  prot="POLR3C"
  #permutations 
  { 
    obs=ew %>% dplyr::select(contains(paste0(prot,"_")),gene_name) %>% filter(.[[1]]>0) %>% pull(gene_name) %>% unique()
    obs=obs[obs %in% colnames(ECdz_min)] # only selecting those exist in our ppi network
    obs=obs[obs %in% names(which(matrix[,prot]==0))] #Removing the already known interactors
    length(obs)
    perms=sapply(c(1:10000),function(x) median(sample(ECdz_min[prot,],size = length(obs)),na.rm = T))
    plot(density(perms),xlim=c(0,0.1),main=prot)
    abline(v=median(ECdz_min[prot,obs],na.rm = T),lty=2)
    (median(ECdz_min[prot,obs],na.rm = T)-mean(perms))/sd(perms)
  }
  
  temp_fn2= function(prot)
  {
    obs=ew %>% dplyr::select(contains(paste0(prot,"_")),gene_name) %>% filter(.[[1]]>0) %>% pull(gene_name) %>% unique() %>% sort()
    obs=obs[obs %in% names(which(matrix[,prot]==0))] # only selecting those that are not known interactors of poi in biogrid.
    # obs=obs[obs %in% names(which(matrix[,prot]==0))] #Removing the already known interactors
    perms=sapply(c(1:100000),function(x) median(sample(ECdz_min[prot,],size = length(obs)),na.rm = T))
    obs= median(ECdz_min[prot,obs],na.rm = T)
    return(list(permutations=perms,observed=obs))
    
  }
  
  masspec_jenny=list(POLR3C= temp_fn2("POLR3C"),
                     POLR3G= temp_fn2("POLR3G"),
                     POLR3GL= temp_fn2("POLR3GL"))
  saveRDS(masspec_jenny,"masspec_jenny.Rdata")
  
  
  
  temp_fn3= function(prot)
  {
    obs=readRDS(paste0("~/Downloads/",prot,"_ppi.RDS")) %>% pull(gene_name) %>% unique() %>% sort()
    obs=obs[obs %in% names(which(matrix[,prot]==0))] # only selecting those that are not known interactors of poi in biogrid.
    # obs=obs[obs %in% names(which(matrix[,prot]==0))] #Removing the already known interactors
    perms=sapply(c(1:100000),function(x) median(sample(ECdz_min[prot,],size = length(obs)),na.rm = T))
    obs= median(ECdz_min[prot,obs],na.rm = T)
    return(list(permutations=perms,observed=obs))
  }
  
  masspec_jenny_2=list(POLR3A= temp_fn3("POLR3A"),
                     POLR1B= temp_fn3("POLR1B"),
                     POLR2B= temp_fn3("POLR2B"))
  saveRDS(masspec_jenny_2,"masspec_jenny_2.Rdata")
  
  
}




