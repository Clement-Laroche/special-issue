library(data.table)
library(rlist)
library(stringr)
library(ClustGeo)
library(udunits2)
library(igraph)
library(sf)
load(file = "Preprocessed_data.Rdata")

# setting the hydrographic network spatial 
# projection to WGS84
st_crs(map_dep$geometry) <- 2154
map_dep$geometry <- st_transform(map_dep$geometry,4326)

# reformating the end and start points of hydrographic segments
map_dep$SEG_DEB <- st_zm(x = map_dep$SEG_DEB) 
map_dep$SEG_FIN <- st_zm(map_dep$SEG_FIN)

# stations that were active between 2017-02-07 and 2017-09-14
sta <- unique(tab_res$sta[tab_res$d >= "2017-02-07" & tab_res$d <= "2017-09-14"])
# creating a table with correct projection of
# the stations coordinates
sta_crs_o <- tab_sta[tab_sta$ï..CdStationMesureEauxSurface %in% sta,]
sta_crs_o <- st_as_sf(x = sta_crs_o, 
                      coords = c("CoordXStationMesureEauxSurface", "CoordYStationMesureEauxSurface"),
                      crs = 2154)
st_crs(sta_crs_o) <- 2154

# creating the graph with :
#      - nodes = segments
#      - edges = link between segments 
#     (if two segments are contiguous, we put an edge between their nodes)
# nodes data frame :
seg_nodes <- as.data.frame(1:nrow(map_dep))
# starting point of edges data frame
fr_om <- matrix(unlist(map_dep$SEG_FIN),ncol = 2,byrow = TRUE)
fr_om <- paste(fr_om[,1],fr_om[,2],sep = ";")
fr_om <- as.data.frame(fr_om)
# ending point of edges data frame
r_f <- matrix(unlist(map_dep$SEG_DEB),ncol = 2,byrow = TRUE)
r_f <- paste(r_f[,1],r_f[,2],sep = ";")
r_f <- as.data.frame(r_f)
# matching the corresponding end points and start points coordinates
t_o <- sapply(fr_om$fr_om, function(value) which(r_f$r_f == value))
FR_OM <- rep(x = 1:nrow(fr_om),unlist(lapply(X = t_o,FUN = "length")))
T_O <- unlist(t_o)
# creating the edges data frame with their end and start points
seg_edges <- as.data.frame(cbind(FR_OM,T_O))
colnames(seg_edges) <- c("from","to") 
# cleaning the edges that don't make sense 
if(sum(is.na(seg_edges$to)) > 0)
{
  seg_edges <- seg_edges[!is.na(seg_edges$to),]
}
if(length(which(seg_edges$from==seg_edges$to)) > 0)
{
  seg_edges <- seg_edges[-which(seg_edges$from==seg_edges$to),]
}
# creating the graph of the hydrographic segments
water_network <- graph_from_data_frame(d = seg_edges, vertices = seg_nodes, directed = FALSE)
# creating the dual graph defined as follows : 
#              - nodes : each joining point between segments is a node
#              - edges : each hydrographic segments between two joining points is an edge
#              - weight : the weight associated to an edge is the length of the segment in meters
# creating the nodes data frame :
nodes_water <- unique(c(unique(r_f$r_f),unique(fr_om$fr_om)))
nodes_water <- as.data.frame(nodes_water)
# creating the edges data frame : 
edge_water <- cbind(fr_om$fr_om,r_f$r_f)
edge_water <- as.data.frame(edge_water)
# adding the weight information to the edges data frame : 
map_dep$geometry <- st_transform(x = map_dep$geometry,crs = 2154)
edge_water$weight <- st_length(map_dep$geometry)
# creating the graph
new_water_network <- graph_from_data_frame(d = edge_water,directed = FALSE,vertices = nodes_water)
# save(new_water_network,file = "graph_seg_tot.Rdata")

# adding the coordinates of each nodes to the node data frame
nodes_water$lat <- str_extract(string = nodes_water$nodes_water,pattern = "(\\d+.\\d+)(?=;)")
nodes_water$lon <- str_extract(string = nodes_water$nodes_water,pattern = "(?<=;)(\\d+.\\d+)")
nodes_water$lon <- as.numeric(nodes_water$lon)
nodes_water$lat <- as.numeric(nodes_water$lat)


# assigning each station to the closest node of new_water_network graph
# working in the CRS Lambert 93
map_dep$geometry <- st_transform(map_dep$geometry,2154)
sta_crs_o <- st_transform(x = sta_crs_o,crs = 2154)
l_closest <- rep(NA,nrow(sta_crs_o))
for(j in 1:nrow(sta_crs_o))
{
  di <- st_distance(x = sta_crs_o$geometry[j],y = map_dep$geometry)
  d <- order(di)[1:10]
  d <- d[which(map_dep$LARGEUR[d] != "Inconnu" & map_dep$BRAS[d] != "Inconnu")]
  l_closest[j] <- r_f$r_f[d[1]]
}
sta_crs_o <- st_transform(x = sta_crs_o,crs = 4326)

# creating the matrix of distances in the hydrographic graph
# between the stations
M_dist <- matrix(NA,length(l_closest),length(l_closest))
for(i in 1:(length(l_closest)-1))
{
  for(j in (i+1):length(l_closest))
  {
    if(!is.na(l_closest[i])&!is.na(l_closest[j]))
    {
      pths <- get.shortest.paths(graph = new_water_network,
                                 from = l_closest[i],
                                 to = l_closest[j],mode = "all",
                                 weights = edge_water$weight,
                                 output = "epath",
                                 predecessors = FALSE,
                                 inbound.edges = FALSE)
      M_dist[i,j] <- sum(pths$epath[[1]]$weight)
      M_dist[j,i] <- sum(pths$epath[[1]]$weight)
    }
  }
}
M_dist[is.na(M_dist)] <- 0
M_dist[M_dist == 0] <- NA
diag(M_dist) <- 0
M_dist <- as.data.frame(M_dist)
colnames(M_dist) <- sta_crs_o$ï..CdStationMesureEauxSurface
row.names(M_dist) <- sta_crs_o$ï..CdStationMesureEauxSurface

# creating the station graph defined as follows :
#              - nodes = stations
#              - edges = existence of a path in the hydrographic network between two stations
#              - weight = length of the path in meters
nodes_sta_graph <- as.data.frame(colnames(M_dist))
fr_om <- c()
t_o <- c()
wei <- c()
for(i in 1:nrow(M_dist))
{
  pos <- which(M_dist[i,]!=0 & !is.na(M_dist[i,]))
  t_o <- c(t_o,colnames(M_dist)[pos])
  fr_om <- c(fr_om,rep(row.names(M_dist)[i],length(pos)))
  wei <- c(wei,as.numeric(M_dist[i,pos]))
}
edges_sta_graph <- as.data.frame(matrix(NA,ncol = 3,nrow = length(wei)))
edges_sta_graph[,1] <- fr_om
edges_sta_graph[,2] <- t_o
edges_sta_graph[,3] <- wei
colnames(edges_sta_graph) <- c("from","to","weight")
# creation of the graph
sta_graph <- graph_from_data_frame(d = edges_sta_graph,directed = FALSE,vertices = nodes_sta_graph)
sta_crs_o$comp <- as.factor(as.numeric(components(sta_graph)$membership))
# save(sta_graph,file = "graph_sta.Rdata")

# clustering from the station graph distance matrix  
M <- as_adjacency_matrix(graph = sta_graph,attr = "weight")
sta_crs_o <- st_transform(sta_crs_o,crs = 2154)
INFO_COMP <- components(sta_graph)
# creation of list with matrix distances of each component of the station graph
M_comp <- as.list(rep(NA,INFO_COMP$no))
for(i in 1:length(M_comp)) 
{
  M_comp[[i]] <- M[row.names(M)%in%names(INFO_COMP$membership)[INFO_COMP$membership == i],colnames(M)%in%names(INFO_COMP$membership)[INFO_COMP$membership == i]]
}
# maximum number of clusters allowed
KMAX <- 30
KLUST_MAT <- M_comp
KLUST_MAT <- lapply(KLUST_MAT,"as.dist")
LAB_MAT <- lapply(KLUST_MAT,"labels")
pos <- which(unlist(lapply(LAB_MAT,"length"))>1)
SCORE_FUN <- rep(0,KMAX)
KK <- as.list(rep(NA,KMAX))
# first iteration of the clustering process
# perform clustering on each component
KK[[1]] <-  lapply(X = KLUST_MAT[pos],FUN = "hclustgeo")
# cut the dendogram (at the first iteration no cut is performed on any component)
KK[[1]] <- lapply(X = KK[[1]],FUN = "cutree",k = 1)
# calculation of the inertia of the clustering
SCORE_FUN[1] <- sum(mapply(FUN = "withindiss",KLUST_MAT[pos],KK[[1]]))
nb_clust <- as.list(rep(NA,KMAX))
# number of clusters in each component
nb_clust[[1]] <- rep(1,length(pos))
for(i in 2:KMAX)
{
  KK[[i]] <- lapply(X = KLUST_MAT[pos],FUN = "hclustgeo")
  wd_score <- rep(0,length(pos))
  for(j in 1:length(pos))
  {
    nb_clust[[i]] <- nb_clust[[i-1]]
    nb_clust[[i]][j] <- nb_clust[[i]][j]+1
    KK1 <- mapply(FUN = "cutree",KK[[i]],k = nb_clust[[i]])
    wd_score[j] <- sum(mapply("withindiss",KLUST_MAT[pos],KK1))
  }
  p <- which.min(wd_score)
  nb_clust[[i]] <- nb_clust[[i-1]]
  # we keep the new cluster in the component 
  # where we observed the biggest drop in total inertia
  nb_clust[[i]][p] <- nb_clust[[i]][p] + 1
  KK[[i]] <- mapply(FUN = "cutree",KK[[i]],nb_clust[[i]])
  SCORE_FUN[i] <- sum(mapply(FUN = "withindiss",KLUST_MAT[pos],KK[[i]]))
}

# selection of the optimal number of clusters with
# the elbow method
elb <- cbind(SCORE_FUN,1:length(SCORE_FUN))
elb[,2] <- elb[,2]+length(M_comp)-1
elb <- as.data.frame(elb)
elbi <- rep(NA,length(SCORE_FUN)-2)
for(l in 2:(length(SCORE_FUN)-1))
{
  ml1 <- lm(formula = elb$SCORE_FUN[1:l]~elb$V2[1:l])
  ml2 <- lm(elb$SCORE_FUN[(l):length(SCORE_FUN)]~elb$V2[(l):length(SCORE_FUN)])
  elbi[l-1] <- sum(ml1$residuals^2)+sum(ml2$residuals^2)
}
K_ELB <- (2:(length(SCORE_FUN)-1))[which.min(elbi)]
# reformating
part_arg <- KK[[K_ELB]]
part_arg <- mapply(FUN = "cbind",KK[[K_ELB]],1:length(KK[[K_ELB]]))
part_arg <- list.rbind(part_arg)
new_clust <- as.numeric(as.factor(paste(part_arg[,1],part_arg[,2],sep = "")))
names(new_clust) <- row.names(part_arg)
# assigning stations to their clusters
sta_crs_o$clust <- new_clust[match(sta_crs_o$ï..CdStationMesureEauxSurface,names(new_clust))]
sta_crs_o$clust[which(is.na(sta_crs_o$clust))] <- max(new_clust)+1
clust_results <- sta_crs_o[,c("ï..CdStationMesureEauxSurface","clust")]
save(clust_results,file = "clust_results.Rdata")
