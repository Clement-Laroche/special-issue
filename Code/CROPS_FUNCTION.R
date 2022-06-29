# functions needed
# sourceCpp("Code/PELT.cpp")

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
