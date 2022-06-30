
# packages needed
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggmap)   
library(tidyr)
library(udunits2)
library(sf)
library(fitdistrplus)
library(Rcpp)
library(rlist)
library(stringr)
library(ClustGeo)
library(igraph)
library(transport)
library(rPref)

#########################################################
#################### Preprocessing ######################
#########################################################


# table containing informations on stations 
tab_sta <- read.csv2(file = "Data/Station.csv")

# loading tab_res and checking the format
tab_res <- read.csv2(file = "Data/Analyse.csv")
tab_res <- as.data.table(tab_res)
tab_res$DatePrel <- as.Date(tab_res$DatePrel)
tab_res <- tab_res[order(tab_res$DatePrel)]
tab_res[["RsAna"]] <- as.numeric(tab_res[["RsAna"]])

# selecting all data from the Centre-Val de Loire region
tab_sta <- data.table(tab_sta)
tab_sta$LbRegion <- as.character(tab_sta$LbRegion)
unique(tab_sta$LbRegion)
region <- "CENTRE-VAL DE LOIRE" 
station <- which(tab_sta$LbRegion == region)
station <- tab_sta$ï..CdStationMesureEauxSurface[station]
pos <- which(tab_res$ï..CdStationMesureEauxSurface %in% station)
tab_res <- tab_res[pos,]
tab_sta <- tab_sta[which(tab_sta$LbRegion == region),]
rm(pos,region,station)

# removing non usable samples
if(any(is.na(tab_res$RsAna)))
{
  tab_res <- tab_res%>%subset(!is.na(RsAna))
}
if(any(tab_res$RsAna==0))
{
  tab_res <- tab_res%>%subset(RsAna != 0)
}

# keeping informative column
tab_res <- tab_res %>% select(c("DatePrel","RsAna","SymUniteMesure","MnemoRqAna","ï..CdStationMesureEauxSurface"))

# creating the column of quantification
tab_res$col <- FALSE
# names(table(tab_res$MnemoRqAna))
tab_res$col[tab_res$MnemoRqAna == "RÃ©sultat < au seuil de quantification"] <- TRUE
tab_res$col[tab_res$MnemoRqAna == "RÃ©sultat < seuil de dÃ©tection"] <- TRUE
# ggplot()+geom_point(tab_res = tab_res,
#                     aes(x = DatePrel,y = RsAna,col = (col==FALSE)))+
#          scale_color_discrete("Quantification")+xlab("Date")+ylab("Concentrations (µg/L)")

# keeping column of interest
colnames(tab_res) <- c("d","C","unite","comment","sta","col")
tab_res <- tab_res[,c("d","C","col","sta")]


# hydrographic information
map_dep1 <- st_read(dsn = "Data/Cher/TRONCON_HYDROGRAPHIQUE.shp")
map_dep2 <- st_read(dsn = "Data/Eure/TRONCON_HYDROGRAPHIQUE.shp")
map_dep3 <- st_read(dsn = "Data/Indre/TRONCON_HYDROGRAPHIQUE.shp")
map_dep4 <- st_read(dsn = "Data/IndreetLoir/TRONCON_HYDROGRAPHIQUE.shp")
map_dep5 <- st_read(dsn = "Data/Loiret/TRONCON_HYDROGRAPHIQUE.shp")
map_dep6 <- st_read(dsn = "Data/Loiretcher/TRONCON_HYDROGRAPHIQUE.shp")

# extracting for each of the table the 
# geographical points that are starting 
# and ending the hydrographical segments


map_dep1$SEG_DEB <- NA
map_dep1$SEG_FIN <- NA
for(i in 1:nrow(map_dep1))
{
  remp <- st_cast(x = map_dep1$geometry[i],to = "POINT")
  map_dep1$SEG_DEB[i] <- remp[1]
  map_dep1$SEG_FIN[i] <- remp[length(remp)]
}
map_dep1 <- data.table(map_dep1)

map_dep2$SEG_DEB <- NA
map_dep2$SEG_FIN <- NA
for(i in 1:nrow(map_dep2))
{
  remp <- st_cast(x = map_dep2$geometry[i],to = "POINT")
  map_dep2$SEG_DEB[i] <- remp[1]
  map_dep2$SEG_FIN[i] <- remp[length(remp)]
}
map_dep2 <- data.table(map_dep2)

map_dep3$SEG_DEB <- NA
map_dep3$SEG_FIN <- NA
for(i in 1:nrow(map_dep3))
{
  remp <- st_cast(x = map_dep3$geometry[i],to = "POINT")
  map_dep3$SEG_DEB[i] <- remp[1]
  map_dep3$SEG_FIN[i] <- remp[length(remp)]
}
map_dep3 <- data.table(map_dep3)

map_dep4$SEG_DEB <- NA
map_dep4$SEG_FIN <- NA
for(i in 1:nrow(map_dep4))
{
  remp <- st_cast(x = map_dep4$geometry[i],to = "POINT")
  map_dep4$SEG_DEB[i] <- remp[1]
  map_dep4$SEG_FIN[i] <- remp[length(remp)]
}
map_dep4 <- data.table(map_dep4)


map_dep5$SEG_DEB <- NA
map_dep5$SEG_FIN <- NA
for(i in 1:nrow(map_dep5))
{
  remp <- st_cast(x = map_dep5$geometry[i],to = "POINT")
  map_dep5$SEG_DEB[i] <- remp[1]
  map_dep5$SEG_FIN[i] <- remp[length(remp)]
}
map_dep5 <- data.table(map_dep5)

map_dep6$SEG_DEB <- NA
map_dep6$SEG_FIN <- NA
for(i in 1:nrow(map_dep6))
{
  remp <- st_cast(x = map_dep6$geometry[i],to = "POINT")
  map_dep6$SEG_DEB[i] <- remp[1]
  map_dep6$SEG_FIN[i] <- remp[length(remp)]
}
map_dep6 <- data.table(map_dep6)

# Joining all informations into a single table
map_dep <- rbind(map_dep1,map_dep2)
map_dep <- rbind(map_dep,map_dep3)
map_dep <- rbind(map_dep,map_dep4)
map_dep <- rbind(map_dep,map_dep5)
map_dep <- rbind(map_dep,map_dep6)
rm(map_dep1,map_dep2,map_dep3,map_dep4,map_dep5,map_dep6,remp,i)

# output 
# save(tab_res,tab_sta,map_dep,file = "Preprocessed_data.Rdata")

#########################################################
################## Temporal segmentation ################
#########################################################

# necessary data for this step
# load(file = "Code/Preprocessed_data.Rdata")

# necessary functions for this step
sourceCpp("Code/PELT.cpp")
source(file = "Code/CROPS_FUNCTION.R")

# creating the daily maximum time series
temp_break <- tab_res%>%group_by(d)%>% summarise(C = max(C),col = min(col[C==max(C)]))
temp_break$col <- as.numeric(as.character(temp_break$col))

# shape parameter estimation 
tab <- cbind(temp_break$C,temp_break$C)
tab <- as.data.frame(tab)
colnames(tab) <- c("left","right")
tab$left[temp_break$col==1] <- 0 
res <- fitdistcens(censdata = tab,distr =  "weibull")
sigma <- as.numeric(res$estimate[1])
rm(tab,res)

# CROPS penalty interval
b_0 <- log(nrow(temp_break))/3
b_1 <- log(nrow(temp_break))*3

# CROPS results
crops_weibull <- CROPS_PELT(y = temp_break[,2:3],beta_0 = b_0,beta_1 = b_1,init = "mle",q_init = 0.9, k = sigma,min_seg = 50)
K_l50 <- unlist(lapply(crops_weibull$segments,FUN = "length"))
crops_weibull$cost <- crops_weibull$cost[-which(duplicated(K_l50))]
crops_weibull$beta <- crops_weibull$beta[-which(duplicated(K_l50))]
crops_weibull$segments <- crops_weibull$segments[-which(duplicated(K_l50))]
K_l50 <- K_l50[-which(duplicated(K_l50))]

# Choosing the optimal segmentations with the elbow method
slope <- rep(NA,length(K_l50)-2)
for(i in 2:(length(K_l50)-1))
{
  ml1 <- lm(crops_weibull$cost[1:i]~K_l50[1:i])
  ml2 <- lm(crops_weibull$cost[(i):length(K_l50)]~K_l50[(i):length(K_l50)])
  slope[i-1] <- sum((ml1$residuals)^2)+sum((ml2$residuals)^2)
}
K_star <- (2:(length(K_l50)-1))[which.min(slope)]

# illustration code 
ml1 <- lm(crops_weibull$cost[1:K_star]~K_l50[1:K_star])
ml2 <- lm(crops_weibull$cost[K_star:length(K_l50)]~K_l50[K_star:length(K_l50)])
g1 <- ggplot()+geom_point(aes(x = K_l50,y = crops_weibull$cost))+geom_abline(slope = ml1$coefficients[2],intercept = ml1$coefficients[1],col = "red")+geom_abline(slope = ml2$coefficients[2],intercept = ml2$coefficients[1],col = "red")+xlab("Number of breaks")+ylab("Cost of segmentation")
d_break <- temp_break$d[crops_weibull$segments[[K_star]]]
data_seg <- tab_res[tab_res$d > d_break[4] & tab_res$d < d_break[5],]
g2 <- ggplot()+geom_point(aes(x = temp_break$d,y = temp_break$C,col = as.factor(temp_break$col==0)))+
  geom_rect(aes(xmin = temp_break$d[c(1,crops_weibull[[3]][[K_star]]+1)],xmax = temp_break$d[c(crops_weibull[[3]][[K_star]],nrow(temp_break))],ymin= -0.1,ymax = 3),col = "white",fill = "grey",alpha = 0.4)+
  xlab("Date")+ylab("Concentration (ug/L)")+labs(col='Quantification')+
  geom_rect(aes(xmin = temp_break$d[crops_weibull[[3]][[K_star]][4]]+1,xmax = temp_break$d[crops_weibull[[3]][[K_star]][5]]),ymin= -0.1,ymax = 3,col = "black",fill = "grey",alpha = 0.0)+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0, units = 'in'))+
  guides(fill = guide_legend(label.position = "bottom"))

# output
# save(crops_weibull,sigma,K_star,file = "res_temp_break.Rdata")

#########################################################
################## Spatial Clustering ###################
#########################################################

# necessary data for this step
# load(file = "Preprocessed_data.Rdata")
# load(file = "res_temp_break.Rdata")

# setting the hydrographic network spatial 
# projection to WGS84
st_crs(map_dep$geometry) <- 2154
map_dep$geometry <- st_transform(map_dep$geometry,4326)

# reformating the end and start points of hydrographic segments
map_dep$SEG_DEB <- st_zm(x = map_dep$SEG_DEB) 
map_dep$SEG_FIN <- st_zm(map_dep$SEG_FIN)

# stations that were active between 2017-02-07 and 2017-09-14
sta <- unique(tab_res$sta[tab_res$d >= temp_break$d[crops_weibull$segments[[K_star]]][4] & tab_res$d <= temp_break$d[crops_weibull$segments[[K_star]]][5]])
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

# output
# save(clust_results,file = "clust_results.Rdata")

#########################################################
################## Anomaly Detection ####################
#########################################################

# necessary data for this step
# load(file = "Preprocessed_data.Rdata")
# load(file = "res_temp_break.Rdata")
# load(file = "clust_results.Rdata")

# necessary functions for this step
# sourceCpp("Code/PELT.cpp")

# temporal segment data selection
data_seg <- tab_res[tab_res$d >= temp_break$d[crops_weibull$segments[[K_star]]][4] & tab_res$d <= temp_break$d[crops_weibull$segments[[K_star]]][5],]
data_seg$col <- as.numeric(data_seg$col)

# Wasserstein distance between station concentrations
C_MAT <- as.list(rep(NA,max(clust_results$clust)))
for(k in 1:max(clust_results$clust))
{
  pos <- which(clust_results$clust == k)
  sta <- clust_results$ï..CdStationMesureEauxSurface[pos]
  mat_c <- matrix(NA,nc = length(sta),nr = length(sta))
  for(i in 1:nrow(mat_c))
  {
    for(j in 1:nrow(mat_c))
    {
      mat_c[i,j] <- wasserstein1d(a = data_seg$C[data_seg$sta == sta[i]],
                                  b = data_seg$C[data_seg$sta == sta[j]],
                                  p = 1)
    }
  }
  row.names(mat_c) <- sta
  colnames(mat_c) <- sta
  C_MAT[[k]] <- mat_c
}

# construction of the two criteria on each cluster
cluster_info <- matrix(NA,nrow = max(clust_results$clust),ncol = 3)
cluster_info <- as.data.frame(cluster_info)
colnames(cluster_info) <- c("cluster","abs","ord")
cluster_info$cluster <- 1:nrow(cluster_info)
for(i in 1:nrow(cluster_info))
{
  pos <- which(data_seg$sta %in% row.names(C_MAT[[i]]))
  abs <- mean(apply(X = C_MAT[[i]],MARGIN = 1,FUN = "mean"))
  ord <- costfcpp(data = data_seg[pos,2:3],
                  tstart = 0,
                  tstop = length(pos)-1,
                  init = "quantile",q_init = 0.5,k = sigma,prec = 10^-8,
                  Nmax = 100,
                  lmax = 100000)[1]
  cluster_info$abs[i] <- abs
  cluster_info$ord[i] <- ord
}
cluster_info$ord <- 1/cluster_info$ord

# anomaly detection through multicriterion optimisation
p <- high(cluster_info$abs) * high(cluster_info$ord) 
res_anomaly <- psel(cluster_info, p, top = nrow(cluster_info))
res_anomaly$.level <- as.factor(res_anomaly$.level)

# plot of the Pareto front
g <- ggplot(data = res_anomaly)+geom_point(aes(x = abs,y = ord,col = .level))

# output
# save(res_anomaly,"res_anom.Rdata")
