

generate_case_control_data <- function(n, p, ncar){
  
  # sample number of cases. The first nD people in the data are cases. 
  nD <- rbinom(1, n, p)
  
  #sample indices of carries:
  ind_g <- sample(1:n, ncar)
  
  ## compute number of diseased carriers:
  ncar_D1 <- sum(ind_g <= nD)
  
  # controls:
  ind_cont <- (nD + 1):n
  
  ncar_D0 <- length(intersect(ind_g, ind_cont))
  
  prop <- nD/n
  
  return(list(n=n, ncar = ncar, 
              ncar_D0 = ncar_D0, ncar_D1 = ncar_D1, prop = prop))
  
}



generate_sampled_control_data <- function(n, p, ncar, control_case_ratio = 1){
  
  # sample number of cases. The first nD people in the data are cases. 
  nD <- rbinom(1, n, p)
  
  #sample indices of carries:
  ind_g <- sample(1:n, ncar)
  
  ## compute number of diseased carriers:
  ncar_D1 <- sum(ind_g <= nD)
  
  # sample controls:
  ind_samp_cont <- sample((nD + 1):n, nD*control_case_ratio)
  
  ncar_D0 <- length(intersect(ind_g, ind_samp_cont))
  
  ncar <- ncar_D1 + ncar_D0
  
  prop <- 1/(1 + control_case_ratio)
  
  return(list(n=nD + nD*control_case_ratio, ncar = ncar, 
              ncar_D0 = ncar_D0, ncar_D1 = ncar_D1, prop = prop))
  
}


compute_tests_two_pops <- function(n1, ncar1, ncar1_D0, ncar1_D1, prop1,
                                           n2, ncar2, ncar2_D0, ncar2_D1, prop2){
  if (ncar1 + ncar2 == 0)  return(list(pval_BR = NA, pval_score = NA, pval_SPA = NA))
  tests <- compute_tests_two_pops_efficient(n1, ncar1, ncar1_D0, ncar1_D1, prop1,
                                            n2, ncar2, ncar2_D0, ncar2_D1, prop2)
  
  if (tests$pval_score > 0.05){
    pval_SPA <- tests$pval_score
  } else{
    nD1 <- as.integer(n1*prop1) # number of diseased in study 1
    nD2 <- as.integer(n2*prop2) # number of disased in study 2
    nD1_non_car <- nD1 - ncar1_D1 # number of diseaesd non-carriers in study 1
    nD2_non_car <- nD2 - ncar2_D1 # number of diseaesd non-carriers in study 2
    
    # organize data for the ScoreTest_SPA function: 
    g <- c(rep(1, ncar1), rep(0, n1-ncar1), # study 1: carriers, then non-carriers
           rep(1, ncar2), rep(0, n2-ncar2)) # study 2: carriers, then non-carriers
    d <- c(rep(1, ncar1_D1), rep(0, ncar1_D0),  # study 1: diseased carriers, non-diseased carriers, 
           rep(1, nD1_non_car), rep(0, (n1-ncar1) - nD1_non_car), # study 1: diseased non-carriers, non-diseaesd non-carriers.
           rep(1, ncar2_D1), rep(0, ncar2_D0),  # study 2: diseased carriers, non-diseased carriers, 
           rep(1,nD2_non_car), rep(0, n2-ncar2 - nD2_non_car) # study 2: diseased non-carriers, non-diseaesd non-carriers.
    )
    
    pval_SPA <- ScoreTest_SPA(genos = g, 
                              pheno = d, 
                              cov = matrix(c(rep(1, n1), rep(0, n2))))$p.value
  }
  
  list(pval_BR = tests[["pval_BR"]], pval_BR_midp = tests[["pval_BR_midp"]], pval_score = tests[["pval_score"]], pval_SPA = pval_SPA)
  
}



compute_tests_two_pops_efficient <- function(n1, ncar1, ncar1_D0, ncar1_D1, prop1,
                                                     n2, ncar2, ncar2_D0, ncar2_D1, prop2){
  tests1 <- compute_tests_single_pop_efficient(n1, ncar1, ncar1_D0, ncar1_D1, prop1)
  tests2 <- compute_tests_single_pop_efficient(n2, ncar2, ncar2_D0, ncar2_D1, prop2)
  
  pval_BR <- calc.br.pval(ncar = tests1$ncar + tests2$ncar, 
                          sum.d = tests1$ncarD + tests2$ncarD, 
                          phat = c(rep(tests1$prop, tests1$ncar), rep(tests2$prop, tests2$ncar)))
  
  score <- tests1$score + tests2$score
  var_score <- tests1$var_score + tests2$var_score
  pval_score <- pchisq(score^2/var_score, df = 1, lower.tail = FALSE)
  
  return(list(pval_BR = pval_BR[["pval"]], pval_BR_midp = pval_BR[["mid.pval"]], pval_score = pval_score))
}





compute_tests_single_pop_efficient <- function(n, ncar, ncar_D0, ncar_D1, prop){
  
  if (ncar == 0){
    return(list(ncar = 0, ncarD = 0, prop = prop, score = 0, var_score =0, pval_score = NA, pval_BR = NA, pval_BR_midp = NA))
  }
  
  score <- (1-prop)*ncar_D1 -prop*ncar_D0
  
  int <- logit(prop)
  s <- exp(int)/((1 + exp(int))^2)
  var_score <- (s*ncar*(n-ncar))/n
  pval_score <- pchisq(score^2/var_score, df = 1, lower.tail = FALSE)
  
  pval_BR <- calc.br.pval(ncar, ncar_D1, rep(prop, ncar))
  
  return(list(ncar = ncar, ncarD = ncar_D1, prop = prop, score = score, var_score = var_score, 
  					pval_score = pval_score, pval_BR = pval_BR[["pval"]], pval_BR_midp = pval_BR[["mid.pval"]]))
  
}




