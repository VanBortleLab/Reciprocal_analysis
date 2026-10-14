
#First crucial step
detach("package:MASS", unload = TRUE) ##Should do  this because there was conflict with its select function while plotting cross-polymerase interaction plots below



 
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
   
   ct=readRDS("complexinfo_CORUM_human.Rdata")
   chosen_complex=  ct %>% filter(`Recommended name` %in% c("Polycomb repressive complex 2/3",
                                                "NuRD.1 complex",
                                                "DNMT3B complex"))
   shared_subs=names(which(sort(table(unlist(chosen_complex$unique_prots)))>1)) #shared subunits
   chosen_complex=chosen_complex %>% rowwise() %>%
     mutate(unique_prots=list(setdiff(unique_prots,shared_subs))) #removing shared subunits
   
   

   
   unique_pol_subs=unlist(chosen_complex$unique_prots)

   
   temp1=matrix[unique_pol_subs,unique_pol_subs]  #actual matrix
   temp2=((temp1==0)* ECdz_min[unique_pol_subs,unique_pol_subs])>=1.5 #Gained interactions
   temp3= temp1>0 & ECdz_min[unique_pol_subs,unique_pol_subs]>=2 #Old/popularity-corrected interactions
   
   ##### Choose one of the below for old or gained interaction
   if(F)
   {
     temp2=temp2 & F # This line removes gained interactions in the plot
     temp3=temp3 & F # This line removes popularity corrected interactions.
   }
   
   # temp3=ECdz_min<2 # No interactions
   

   pol1_subs= chosen_complex$unique_prots[1] %>% unlist()
   pol1_subs_abb=pol1_subs
   pol2_subs=chosen_complex$unique_prots[2] %>% unlist()
   pol2_subs_abb= pol2_subs
   pol3_subs=chosen_complex$unique_prots[3] %>% unlist()
   pol3_subs_abb= pol3_subs
   


   #plot it 
   {
     # --- 1. Build node list ---
     nodes = bind_rows(
       data.frame(id = pol1_subs, group = "POL1", label = pol1_subs_abb),
       data.frame(id = pol2_subs, group = "POL2", label = pol2_subs_abb),
       data.frame(id = pol3_subs, group = "POL3", label = pol3_subs_abb)
     )
     
     #making equilateral triagnle coordinate
     side <- 250
     h <- sqrt(3) / 2 * side
     
     pos = list(
       POL1 = c(x = 0,         y = -2*h/3),
       POL2 = c(x = -side/2,   y =  h/3),
       POL3 = c(x =  side/2,   y =  h/3)
       )
     
    cluster_radius=60 #radius of individual complexes in the plot.
     set.seed(42)
     nodes = nodes %>%
       group_by(group) %>%
       mutate(
         n     = n(),
         idx   = row_number(),
         angle = 2 * pi * idx / n,
         x     = pos[[group[1]]]["x"] + cluster_radius * cos(angle),
         y     = pos[[group[1]]]["y"] + cluster_radius * sin(angle)
       ) %>%
       ungroup() %>%
       mutate(
         color   = case_when(
           group == "POL1" ~ "#EEDDAA",
           group == "POL2" ~ "#FFAAAA",
           group == "POL3" ~ "#AADDFF"
         ),
         # shape="circle",
         # margin=13,
         physics = FALSE,
         font.size = 11,
         font.vadjust = -40,
         font.color="black"
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
 
 