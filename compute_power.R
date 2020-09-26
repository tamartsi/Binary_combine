
### given a vector of disease probabilities (assumed under the null) for each of n individuals, 
### an odds ratio parameter
### a vector of p-value thresholds, and number of carriers, computes power for BinomiRare test
### (both p-value and mid-p-value). 
### Power is computed by sampling of individuals. Default number of samples is 1000.

compute_power_BinomiRare <- function(disease_prob, 
									 OR,
									 pval_thresh_vec, 
									 n_carrier,
									 n_samp = 1000){
									 
	
	n <- length(disease_prob)
	
	samp_pvals <- data.frame(pval = rep(NA, n_samp), 
							 mid.pval = rep(NA, n_samp))
	
	for (i in 1:n_samp){
		set.seed(i)
		samp_inds <- sample(1:n, n_carrier)
		phat_i_null <- disease_prob[samp_inds]
		phat_i_alt <- expit(logit(phat_i_null) + log(OR)) 
		pvals <- calc.br.pval(n_carrier, 
							  sum.d = sum(rbinom(n_carrier, 1, phat_i_alt)), 
							  phat = phat_i_null,
							  return.both.p = TRUE)
		samp_pvals[i,c("pval", "mid.pval")] <- pvals[c("pval", "mid.pval")]
	}
	
	power_by_thresh <- data.frame(pval_threshold = pval_thresh_vec,
								  power_pval = NA,
								  power_mid_pval = NA)
	
	for (i in 1:nrow(power_by_thresh)){
		power_by_thresh[i,"power_pval"] <- 
				mean(samp_pvals[,"pval"] <= power_by_thresh[i,"pval_threshold"])
		power_by_thresh[i,"power_mid_pval"] <- 
				mean(samp_pvals[,"mid.pval"] <= power_by_thresh[i,"pval_threshold"])

	}
	
	if (n_carrier == 1){
		power_by_thresh$power_mid_pval <- NULL
	}
	
	return(power_by_thresh)
	
}

# test: 
# n <- 1000
# disease_prob_null <- runif(n, 0.1, 0.4)
# disease_status <- rbinom(n, 1, disease_prob)
# pval_thresh_vec <- c(0.05, 0.01, 1e-3, 1e-4)
# compute_power_BinomiRare(disease_prob, 
						 # OR = 4.7, 
						 # pval_thresh_vec,
						 # n_carrier = 10)
						 
						 
						 
						 
