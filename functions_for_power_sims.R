

generate_case_control_data <- function(n, p, ncar, effect_size, model = "logistic"){
  
  # sample disease data. For computational efficiency, assume first 
  # ncar people are the carriers of the rare allele of interest. 
  
  if (model == "logistic"){
  	d <- c(rbinom(ncar, 1, expit(logit(p) + log(effect_size))), rbinom(n-ncar, 1, p))
  } else{ # model == "linear"
  	d <- c(rbinom(ncar, 1, p + effect_size), rbinom(n-ncar, 1, p))
  	}
  
  return(list(d=d, carrier_inds = 1:ncar))
  
}





### using only summaries -- not generating g vector. 
compute_tests_single_pop_efficient <- function(d, g_inds, pD = NULL){
  
  if (length(g_inds) == 0){
    return(list(ncar = 0, ncarD = 0, prop = mean(d), score = 0, var_score =0, pval_score = NA, pval_BR = NA, pval_BR_midp = NA))
  }
  
  n <- length(d)
  
  if (is.null(pD)){
  	prop <- mean(d)
  } else{ # disease probability in the population is provided. 
  	prop <- pD
  	}
  
  
  ncar <- length(g_inds)
  ncarD <- sum(d[g_inds])
  score <- sum((d - prop)[g_inds])
  
  int <- logit(prop)
  s <- exp(int)/((1 + exp(int))^2)
  var_score <- (s*ncar*(n-ncar))/n
  pval_score <- pchisq(score^2/var_score, df = 1, lower.tail = FALSE)
  
  pvals_BR <- calc.br.pval(ncar, ncarD, rep(prop, ncar))
  
  return(list(ncar = ncar, ncarD = ncarD, prop = prop, score = score, var_score = var_score, 
              pval_score = pval_score, pval_BR = pvals_BR[["pval"]], pval_BR_midp = pvals_BR[["mid.pval"]]))
  
}


compute_tests_two_pops_efficient <- function(d1, g_inds1, d2, g_inds2, pD1 = NULL, pD2 = NULL){
  tests1 <- compute_tests_single_pop_efficient(d1, g_inds1, pD = pD1)
  tests2 <- compute_tests_single_pop_efficient(d2, g_inds2, pD = pD2)
  
  if (tests1$ncar + tests2$ncar == 0) return(list(pval_BR = NA, pval_BR_midp = NA, pval_score = NA))
  
  pvals_BR <- calc.br.pval(ncar = tests1$ncar + tests2$ncar, 
                          sum.d = tests1$ncarD + tests2$ncarD, 
                          phat = c(rep(tests1$prop, tests1$ncar), rep(tests2$prop, tests2$ncar)))
  
  score <- tests1$score + tests2$score
  var_score <- tests1$var_score + tests2$var_score
  pval_score <- pchisq(score^2/var_score, df = 1, lower.tail = FALSE)
  
  return(list(pval_BR = pvals_BR[["pval"]], pval_BR_midp = pvals_BR[["mid.pval"]], pval_score = pval_score))
}


# includes option to take known disease probabilities (rather than compute them from the data).
# this will only apply for score and BinomiRare tests, not for SPA!
compute_tests_single_pop <- function(d, g_inds, pD = NULL){
  
  tests <- compute_tests_single_pop_efficient(d, g_inds, pD = pD)
  
  if (is.na(tests$pval_score)){
    return(list(pval_BR = tests$pval_BR, pval_BR_midp = tests$pval_BR_midp, 
                pval_score = tests$pval_score, pval_SPA = pval_SPA))
  }
  
  if (tests$pval_score > 0.05){
    pval_SPA <- tests$pval_score
  } else{
    g <- rep(0, length(d))
    g[g_inds] <- 1
    pval_SPA <- ScoreTest_SPA(genos = g, pheno = d)$p.value
  }
  list(pval_BR = tests$pval_BR, pval_BR_midp = tests$pval_BR_midp, pval_score = tests$pval_score, pval_SPA = pval_SPA)
  
}

# includes option to take known disease probabilities (rather than compute them from the data).
# this will only apply for score and BinomiRare tests, not for SPA!
compute_tests_two_pops <- function(d1, g_inds1, d2, g_inds2, pD1 = NULL, pD2 = NULL){
  if (length(c(g_inds1, g_inds2)) == 0)  return(list(pval_BR = NA, pval_BR_midp = NA, pval_score = NA, pval_SPA = NA ))
  tests <- compute_tests_two_pops_efficient(d1, g_inds1, d2, g_inds2, pD1 = pD1, pD2 = pD2)
  
  if (tests$pval_score > 0.05){
    pval_SPA <- tests$pval_score
  } else{
    d <- c(d1, d2)
    g1 <- rep(0, length(d1))
    g1[g_inds1] <- 1
    g2 <- rep(0, length(d2))
    g2[g_inds2] <- 1 
    g <- c(g1, g2)
    pval_SPA <- ScoreTest_SPA(genos = g, 
                              pheno = d, 
                              cov = matrix(c(rep(1, length(d1)), rep(0, length(d2)))))$p.value
  }
  
  list(pval_BR = tests[["pval_BR"]], pval_BR_midp = tests[["pval_BR_midp"]], pval_score = tests$pval_score, pval_SPA = pval_SPA)
  
}



sample_controls <- function(d, g_inds, cont_case_ratio = 3){
  num_case <- sum(d)
  g <- rep(0, length(d))
  g[g_inds] <- 1
  inds_case <- which(d == 1)
  inds_cont <- sample(setdiff(1:length(d), inds_case), size = num_case*cont_case_ratio)
  
  return(list(d = d[c(inds_case, inds_cont)], carrier_inds = which(g[c(inds_case, inds_cont)] == 1)))
  
}

