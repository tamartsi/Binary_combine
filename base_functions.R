require(poibin)

expit <- function(x){
  exp(x)/(1 + exp(x))
}

logit <- function(x){
  log(x/(1-x))
}


### 
calc.br.pval <- function(ncar, sum.d, phat){
  d.poibin <- dpoibin(0:ncar, phat)
  prob.cur <- d.poibin[sum.d + 1]
  pval <- prob.cur + sum(d.poibin[d.poibin < prob.cur])
  mid.pval <- pval - prob.cur/2
  
  return(c(pval = pval, mid.pval = mid.pval))


}

