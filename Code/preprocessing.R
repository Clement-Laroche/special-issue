#########################################################
#################### Preprocessing ######################
#########################################################

# packages needed
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(sf)
library(udunits2)

# table containing informations on stations 
tab_sta <- read.csv2(file = "Data/Stations.csv")

# loading tab_res and checking the format
tab_res <- read.csv2(file = "Data/Analyses.csv")
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
map_dep1 <- st_read(dsn = "Data/TRONCON_HYDROGRAPHIQUE.shp")
map_dep2 <- st_read(dsn = "Data/TRONCON_HYDROGRAPHIQUE.shp")
map_dep3 <- st_read(dsn = "Data/TRONCON_HYDROGRAPHIQUE.shp")
map_dep4 <- st_read(dsn = "Data/TRONCON_HYDROGRAPHIQUE.shp")
map_dep5 <- st_read(dsn = "Data/TRONCON_HYDROGRAPHIQUE.shp")
map_dep6 <- st_read(dsn = "Data/TRONCON_HYDROGRAPHIQUE.shp")

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
save(tab_res,tab_sta,map_dep,file = "Preprocessed_data.Rdata")
