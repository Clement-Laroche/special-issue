library(transport)
library(Rcpp)
library(rPref)
library(ggplot2)
sourceCpp("PELT.cpp")
load("clust_results.Rdata")
load("Preprocessed_data.Rdata")

# sigma found in temporal_detect.R
sigma <- 0.3030359

# temporal segment data selection
data_seg <- tab_res[tab_res$d >= "2017-02-07" & tab_res$d <= "2017-09-14",]
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
