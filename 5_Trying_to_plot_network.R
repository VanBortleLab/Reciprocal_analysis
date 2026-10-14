library(tidyverse)
library(Matrix)
library(pheatmap)           # or ComplexHeatmap
library(igraph)
library(ggraph)             # for pretty network plots
# install.packages("uwot")  # if not installed
library(uwot)   # UMAP
library(scales) 
library(ggrepel)



#########################################################
#########################################################
#### Making PPI network and starburst plots.
#########################################################
#########################################################
{

#Making initial things ready
library(igraph)
library(ggraph)
library(scales)
  
ECdz_min=readRDS("ECdz_min.Rdata")
matrix=readRDS("interaction_matrix.Rdata")
temp3=ECdz_min
table(colnames(matrix)==colnames(temp3))
temp3[temp3<=2]=0
temp3[temp3!=0]=1
temp3[matrix==0 & temp3==1]=-1 #temp3 is our interaction matrix where 0 is no edge, 1 is old edge and -1 is new edge



## 1) Build an undirected graph from temp3 (keep -1 if either side says -1)
A <- as.matrix(temp3); A[is.na(A)] <- 0; diag(A) <- 0
new_mask <- (A == -1) | (t(A) == -1)
old_mask <- (A ==  1) | (t(A) ==  1)
remove(temp3)

M <- matrix(0, nrow=nrow(A), ncol=ncol(A), dimnames=dimnames(A))
M[new_mask] <- -1
M[!new_mask & old_mask] <- 1

ij   <- which(M != 0 & upper.tri(M), arr.ind = TRUE)  # undirected = upper triangle
edges_df <- data.frame(
  from  = rownames(M)[ij[,1]],
  to    = colnames(M)[ij[,2]],
  is_new= (M[ij] == -1),
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(edges_df, directed = FALSE,
                           vertices = data.frame(name = rownames(M)))

## 2) Node degree (full graph) 
V(g)$deg_full <- degree(g)
g <- delete_vertices(g, V(g)[degree(g) == 0]) #removing nodes with 0 degree.


## 3) Finding layout/coordinates of the nodes (this is the super slow step, skip if you halready had done this before and load it from past)
if(F){set.seed(1)
xy <- layout_with_lgl(g) # same as create_layout(..., "lgl") but returns a matrix
saveRDS(xy,"layout.Rdata")} ##DO this One time job only.

xy=readRDS("layout.Rdata")
V(g)$x <- xy[,1]
V(g)$y <- xy[,2]

## 4) Sampling only 5k from old edges, to speed up preliminary netowrk. For final picture, we plot all the old edges.
# bg_idx <- which(!E(g)$is_new)
# bg_cap <- 5000L
# E(g)$bg_flag <- FALSE
# if (length(bg_idx)) {
#   E(g)$bg_flag[sample(bg_idx, min(bg_cap, length(bg_idx)))] <- TRUE
# }



##5) Selecting which node to label
V(g)$show_node_text <- FALSE
indx <- order(V(g)$deg_full, decreasing = TRUE)[1:10]
V(g)$show_node_text[indx] <- TRUE


## 6) Plot (all proteins; NEW edges blue; OLD edges faint & sampled)
library(scales)
library(ggrepel)
col_new="#FB8451"
col_poi="#E8E8E8"
col_old="#487A9C"

netwo=ggraph(g, layout = "manual", x = x, y = y) +
  geom_edge_link0(aes(filter = !is_new),
                  colour = col_old, alpha = 0.25, width = 0.1) +
  geom_edge_link0(aes(filter = is_new),
                  colour = col_new, alpha = 0.25, width = 0.1) +
  # top-50 hubs only, size rescaled for big contrast
  geom_node_point(aes(size = 1.015^deg_full,filter=(deg_full>=100)),   #only showing nodes with degree greater than 100
                  colour = col_poi, alpha = 0.55, show.legend = FALSE) +
  scale_size_continuous(range = c(0.01, 5),
                        name = "Degree",
                        breaks = 1.015 ^ (c(1:9)*100),   # choose degree values
                        labels = c(1:9)*100)        +
  theme(legend.position = c(0.8, 0.5))+
  geom_node_point(aes(size = 1.015^deg_full,filter=show_node_text),shape=21,colour="black",stroke=0.3,
                  fill = col_poi, alpha = 0.8, show.legend = T) +
  # labels for top-50 hubs
  # geom_node_text(aes(label = name, filter = show_node_text),
  #                size = 1, repel = F, alpha = 0.9,force=5) +
  geom_text_repel(
    data = igraph::as_data_frame(g, what = "vertices") |> dplyr::filter(show_node_text),
    aes(x = x, y = y, label = name),
    size = 2, force = 50, max.overlaps = Inf,
    segment.colour = "grey40", segment.size = 0.3,fontface="bold",
  ) +
  # geom_node_label(aes(label = name, filter = show_node_text),
  #                 size = 2, repel = TRUE, alpha = 0.9,
  #                 label.size = 0.15, fill = "transparent",force=5)+
  theme_void() +
  ggtitle("PPI: Top 10 hubs highlighted (others hidden)")
print(netwo)
ggsave("ppi_network.png", width = 8, height = 6, units = "in", dpi = 600,bg = "white")


####################################################################################
########################################################################
##Making STARBURST plot for different proteins########################
########################################################################
####################################################################################
# ── 2) Attach ECdz_min scores onto existing edges ─────────────────────────────
edge_ends <- ends(g, E(g), names = TRUE)  # returns protein names, not indices
E(g)$score <- apply(ends(g, E(g)), 1, function(e) ECdz_min[e[1], e[2]])
E(g)$score_norm <- pmin(pmax((E(g)$score - 2) / (20 - 2), 0), 1)

xy  <- readRDS("layout.Rdata")
V(g)$x <- xy[, 1]
V(g)$y <- xy[, 2]
# --- 1) Flags for RHOG-related edges & partners ------------------------------
##############blue orange for old and gained interactions
poi <- "POLR3GL"
highlight_node=NULL
label_size_new= 6
label_size_old=6

############### Plotting Starburst by Coloring as per old and new
{
# edges touching target and is it also highlighted
  E(g)$poi_related <- FALSE
  E(g)$is_highlighted     <- FALSE
  
  if (poi %in% V(g)$name) {
    inc <- incident(g, poi, mode = "all")
    E(g)$poi_related[inc] <- TRUE
    
    # among those edges, check if the OTHER endpoint is in highlight_node
    edge_ends <- ends(g, inc, names = TRUE)
    partners  <- ifelse(edge_ends[, 1] == poi, edge_ends[, 2], edge_ends[, 1])
    E(g)$is_highlighted[inc] <- partners %in% highlight_node
  }
  if(is.null(highlight_node)) {E(g)$is_highlighted=T} #HIhglighting everything if nothing given to highlihght
  
# split those edges by new/old
new_inc <- which(E(g)$poi_related & E(g)$is_new)
old_inc <- which(E(g)$poi_related & !E(g)$is_new)

# partners on each side (exclude poi itself)
new_partners <- setdiff(unique(as.vector(ends(g, new_inc))), poi)
old_partners <- setdiff(unique(as.vector(ends(g, old_inc))), poi)

# vertex flags
V(g)$is_poi         <- V(g)$name %in% poi
V(g)$is_new_partner  <- V(g)$name %in% new_partners
V(g)$is_old_partner  <- V(g)$name %in% old_partners
# V(g)$is_rhog_ego     <- V(g)$is_poi | V(g)$is_new_partner | V(g)$is_old_partner
V(g)$is_highlighted  <- if(!is.null(highlight_node)){V(g)$name %in% highlight_node} else T


# --- 2) Plot: only poi edges; highlight/label ONLY new-edge partners --------
newto2=ggraph(g, layout = "manual", x = x, y = y) +  theme_void() +
  # OLD edges to/from poi (whitesmoke)
  geom_edge_link(
    aes(filter     = poi_related & !is_new & is_highlighted,
        alpha      = score_norm,
        edge_width = score_norm),
    colour = col_old
  ) +
  scale_edge_alpha(range = c(0.5, 1), guide = "none") +
  
  # NEW edges (orange), alpha and width driven by score  
  geom_edge_link(
    aes(filter     = poi_related & is_new & is_highlighted,
        alpha      = score_norm,
        edge_width = score_norm),
    colour = col_new
  ) +
  scale_edge_width(range = c(0.5, 1.5), guide = "none")+
  
  # Old nodes
  geom_node_point(aes(filter =  is_old_partner & is_highlighted),
                  size = 1.2, colour =col_old, alpha = 0.9, show.legend = FALSE) +
  # New nodes
  geom_node_point(aes(filter = is_new_partner & is_highlighted),
                  colour = col_new, alpha = 0.95, show.legend = F,size=2) +
  
  
  #label highlihghted nodes but old
  # geom_node_text(aes(label = name, filter = is_highlighted & is_old_partner),
  #                size = label_size_old, colour = "black", alpha = 0.9, vjust = -0.5,fontface="bold") +
  geom_node_label(
    aes(label = name, filter = is_old_partner & is_highlighted),
    size = label_size_old*0.5, fontface = "bold",fill = col_old, colour = "black",    
    label.r = unit(0.3, "lines"), label.padding = unit(0.15, "lines"),alpha = 0.5,vjust=-0.8 )+
  
  
  # Label only new nodes
  # geom_node_text(aes(label = name, filter = is_new_partner & (deg_full>0) & is_highlighted),
  #                size = label_size_new, colour = "black", alpha = 0.9, vjust = -0.5,fontface="bold") +
  geom_node_label(
    aes(label = name, filter = is_new_partner & (deg_full>0) & is_highlighted),
    size = label_size_new*0.5, fontface = "bold",fill = col_new, colour = "black",    
    label.r = unit(0.3, "lines"), label.padding = unit(0.15, "lines"),alpha = 0.7,vjust=-0.8 )+
  

  
  # If you also want to label poi, uncomment:
  # geom_text_repel(
  #   data = igraph::as_data_frame(g, what = "vertices")[V(g)$is_poi, ],
  #   aes(x = x, y = y, label = name),
  #   size = 3, force = 12, max.overlaps = Inf,
  #   segment.colour = "grey40", segment.size = 0.4, fontface = "bold"
  # ) +
  # RHOG itself: emphasize
  geom_node_point(aes(filter = is_poi), size = 4,
                  colour = col_poi, alpha = 0.9, show.legend = FALSE,) +
  ggtitle(paste0(poi," - centered PPI"))+
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))

print(newto2)
# ggsave(paste0("ppi_network_",poi,".png"), plot=newto2,width = 8, height = 6, units = "in", dpi = 600)

}##luka

############### Plotting Starburst by Adding label colors as per polymerase
highlight_node = if(poi %in% pol1_subs) {c(pol2_subs, pol3_subs)} else if(poi %in% pol2_subs) {c(pol3_subs, pol1_subs)} else {c(pol2_subs, pol1_subs)}
col_pol1= "#B5503A"
col_pol2= "#C9A641"
col_pol3= "#557999"
label_size_new= 10
label_size_old=10
label_size_poi=10
{
  # edges touching target and is it also highlighted
  E(g)$poi_related <- FALSE
  E(g)$is_highlighted     <- FALSE
  
  if (poi %in% V(g)$name) {
    inc <- incident(g, poi, mode = "all")
    E(g)$poi_related[inc] <- TRUE
    
    # among those edges, check if the OTHER endpoint is in highlight_node
    edge_ends <- ends(g, inc, names = TRUE)
    partners  <- ifelse(edge_ends[, 1] == poi, edge_ends[, 2], edge_ends[, 1])
    E(g)$is_highlighted[inc] <- partners %in% highlight_node
  }
  if(is.null(highlight_node)) {E(g)$is_highlighted=T} #Highlighting everything if nothing given to highlight
  
  # split those edges by new/old
  new_inc <- which(E(g)$poi_related & E(g)$is_new)
  old_inc <- which(E(g)$poi_related & !E(g)$is_new)
  
  # partners on each side (exclude poi itself)
  new_partners <- setdiff(unique(as.vector(ends(g, new_inc))), poi)
  old_partners <- setdiff(unique(as.vector(ends(g, old_inc))), poi)
  
  # vertex flags
  V(g)$is_poi         <- V(g)$name %in% poi
  V(g)$is_new_partner  <- V(g)$name %in% new_partners
  V(g)$is_old_partner  <- V(g)$name %in% old_partners
  V(g)$is_highlighted  <- if(!is.null(highlight_node)){V(g)$name %in% highlight_node} else T
  
  # NEW: Assign polymerase group to each vertex
  V(g)$pol_group <- NA
  V(g)$pol_group[V(g)$name %in% pol1_subs] <- "pol1"
  V(g)$pol_group[V(g)$name %in% pol2_subs] <- "pol2"
  V(g)$pol_group[V(g)$name %in% pol3_subs] <- "pol3"
  
  # NEW: Assign label background color based on pol group
  V(g)$label_color <- NA
  V(g)$label_color[V(g)$pol_group == "pol1"] <- col_pol1
  V(g)$label_color[V(g)$pol_group == "pol2"] <- col_pol2
  V(g)$label_color[V(g)$pol_group == "pol3"] <- col_pol3
  
  # For nodes not in any pol group, use old/new colors
  V(g)$label_color[is.na(V(g)$label_color) & V(g)$is_old_partner] <- col_old
  V(g)$label_color[is.na(V(g)$label_color) & V(g)$is_new_partner] <- col_new
  
  # --- 2) Plot: only poi edges; highlight/label ONLY new-edge partners --------
  newto2 <- ggraph(g, layout = "manual", x = x, y = y) +  
    theme_void() +
    # OLD edges to/from poi (whitesmoke)
    geom_edge_link(
      aes(filter     = poi_related & !is_new & is_highlighted,
          alpha      = score_norm,
          edge_width = score_norm),
      colour = col_old
    ) +
    scale_edge_alpha(range = c(0.5, 1), guide = "none") +
    
    # NEW edges (orange), alpha and width driven by score  
    geom_edge_link(
      aes(filter     = poi_related & is_new & is_highlighted,
          alpha      = score_norm,
          edge_width = score_norm),
      colour = col_new
    ) +
    scale_edge_width(range = c(0.5, 1.5), guide = "none")+
    
    # Old nodes
    geom_node_point(aes(filter =  is_old_partner & is_highlighted),
                    size = 1.2, colour =col_old, alpha = 0.9, show.legend = FALSE) +
    # New nodes
    geom_node_point(aes(filter = is_new_partner & is_highlighted),
                    colour = col_new, alpha = 0.95, show.legend = F, size=2) +
    
    # Label old partners with pol-specific colors
    geom_node_label(
      aes(label = name, filter = is_old_partner & is_highlighted, fill = I(label_color)),
      size = label_size_old*0.5, fontface = "bold", colour = "black",    
      label.r = unit(0.3, "lines"), label.padding = unit(0.15, "lines"),
      alpha = 0.5, vjust=-0.8) +
    
    # Label new partners with pol-specific colors
    geom_node_label(
      aes(label = name, filter = is_new_partner & (deg_full>0) & is_highlighted, fill = I(label_color)),
      size = label_size_new*0.5, fontface = "bold", colour = "black",    
      label.r = unit(0.3, "lines"), label.padding = unit(0.15, "lines"),
      alpha = 0.7, vjust=-0.8) +
    
    # poi itself: emphasize
    geom_node_point(aes(filter = is_poi), size = 4,
                    colour = col_poi, alpha = 0.9, show.legend = FALSE) +
    ggtitle(paste0(poi," - centered PPI")) +
    theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))+
    geom_node_label(
      aes(label = name, filter = is_poi, fill = I(label_color)),
      size = label_size_poi*0.5, fontface = "bold", colour = "black",    
      label.r = unit(0.3, "lines"), label.padding = unit(0.15, "lines"),
      alpha = 0.5, vjust=-0.8)
  
  print(newto2)
  # ggsave(paste0("ppi_network_",poi,".png"), plot=newto2,width = 8, height = 6, units = "in", dpi = 600)
  
}##luka




}##luka##
