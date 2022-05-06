library(data.table)
library(dplyr)
library(fitdistrplus)
library(ggplot2)
library(Rcpp)

# loading Rcpp functions
sourceCpp("Code/PELT.cpp")

# loading data
load(file = "Code/Preprocessed_data.Rdata")

# CROPS code
CROPS_PELT <- function(y,beta_0,beta_1,init,q_init,k,min_seg)
{
  RESBETA <- c(beta_0,beta_1)
  res_0 <- peltcpp(sumstat = y,minseg = min_seg,pen = beta_0,init = init,q_init = q_init,k = k,prec = 10^-8,Nmax = 100,lmax = 100000)
  res_1 <- peltcpp(sumstat = y,minseg = min_seg,pen = beta_1,init = init,q_init = q_init,k = k,prec = 10^-8,Nmax = 100,lmax = 100000)
  cp0 <- sort(res_0)
  cp0 <- cp0[-1]
  cp1 <- sort(res_1)
  cp1 <- cp1[-1]
  RES <- c(list(cp0),list(cp1))
  q_0 <- sum(mapply(FUN = "costfcpp",c(0,(cp0+1)),c(cp0,nrow(y)-1),MoreArgs = list("data" = y,"init" = init,"q_init" = q_init,"k" = k,"prec" = 10^-8,"Nmax" = 100,"lmax" = 100000))[2,])
  q_1 <- sum(mapply(FUN = "costfcpp",c(0,(cp1+1)),c(cp1,nrow(y)-1),MoreArgs = list("data" = y,"init" = init,"q_init" = q_init,"k" = k,"prec" = 10^-8,"Nmax" = 100,"lmax" = 100000))[2,])
  COST <- c(q_0,q_1)
  setbeta <- list(c(beta_0,beta_1))
  while(length(setbeta)!=0)
  {
    betaz <- setbeta[[1]]
    pos0 <- which(RESBETA==betaz[1])
    pos1 <- which(RESBETA==betaz[2])
    if(length(RES[[pos0]]) > length(RES[[pos1]])+1)
    {
      cp0 <- RES[[pos0]]
      cp1 <- RES[[pos1]]
      q_0 <- COST[[pos0]]
      q_1 <- COST[[pos1]]
      beta_int <- (q_1-q_0)/(length(cp0)-length(cp1))
      if(beta_int < RESBETA[[pos1]] & beta_int > RESBETA[[pos0]])
      {
        res_int <- peltcpp(sumstat = y,minseg = min_seg,pen = beta_int,init = init,q_init = q_init,k = k,prec = 10^-8,Nmax = 100,lmax = 100000)
        res_int <- sort(res_int)
        res_int <- res_int[-1]
        cost_int <- sum(mapply(FUN = "costfcpp",c(0,(res_int+1)),c(res_int,nrow(y)-1),MoreArgs = list("data" = y,"init" = init,"q_init" = q_init,"k" = k,"prec" = 10^-8,"Nmax" = 100,"lmax" = 100000))[2,])
        COST <- c(COST,cost_int)
        RESBETA <- c(RESBETA,beta_int)
        RES <- c(RES,list(res_int))
        if(length(res_int)!=length(RES[[pos1]]))
        {
          setbeta <- c(setbeta,list(c(RESBETA[pos0],beta_int)),list(c(beta_int,RESBETA[pos1])))
        }
      }
    }
    setbeta <- setbeta[-1]
  }
  return(list("beta" = RESBETA[order(RESBETA)],"cost" = COST[order(RESBETA)],"segments" = RES[order(RESBETA)]))
}

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
data_seg <- data[data$d > d_break[4] & data$d < d_break[5],]
g2 <- ggplot()+geom_point(aes(x = temp_break$d,y = temp_break$C,col = as.factor(temp_break$col==0)))+
  geom_rect(aes(xmin = temp_break$d[c(1,crops_weibull[[3]][[K_star]]+1)],xmax = temp_break$d[c(crops_weibull[[3]][[K_star]],nrow(temp_break))],ymin= -0.1,ymax = 3),col = "white",fill = "grey",alpha = 0.4)+
  xlab("Date")+ylab("Concentration (ug/L)")+labs(col='Quantification')+
  geom_rect(aes(xmin = temp_break$d[crops_weibull[[3]][[K_star]][4]]+1,xmax = temp_break$d[crops_weibull[[3]][[K_star]][5]]),ymin= -0.1,ymax = 3,col = "black",fill = "grey",alpha = 0.0)+
  theme(legend.position="bottom",
        legend.spacing.x = unit(0, units = 'in'))+
  guides(fill = guide_legend(label.position = "bottom"))
