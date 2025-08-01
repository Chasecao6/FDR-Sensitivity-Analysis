rm(list = ls())

load("JL14_pvalue.rda") # Load the p-values from medical journal in Jager and Leek (2014)
source("fin_source.R")

pval_mat = matrix(as.numeric(pvalueData[,1:4]), ncol = 4)
colnames(pval_mat) <- colnames(pvalueData)[1:4]
rownames(pval_mat) <- rownames(pvalueData)

# Get the journals and years
journals <- c("Lancet", "JAMA", 
              "New England Journal of Medicine", 
              "BMJ", "American Journal of Epidemiology")
years <- 2000:2010

# Get the percent of P-values less than 0.05 by journal/year

num_mat <- percent05 <- matrix(NA, nrow = length(journals), ncol = length(years))

rownames(num_mat) <- rownames(percent05) <- journals
colnames(num_mat) <- colnames(percent05) <- years

for(i in 1:length(journals)){
  for(j in 1:length(years)){
    ind <- pval_mat[,4] == years[j] & rownames(pval_mat) == journals[i]
    num_mat[i,j] <- sum(ind)
    percent05[i,j] = mean(as.numeric(pval_mat[ind, 1]) < 0.05, na.rm = T)
  }
}


l <- 1 # for Lancet
ind <- rownames(pval_mat) == journals[l] # from a particular journal over all years
all_pval <- pval_mat[ind, 1] 
n <- length(all_pval)

sig_lv <- 0.05 # significance level
lam_vec <- c(0.01,0.015, 0.025, 0.035, 0.045) # sequence of lambda
fdr_vec <- rep(NA, length(lam_vec))
rho_vec <- seq(1.0, 0.1, by = -0.1)
bvec <- seq(1, 4, by = 1)

pi0 <- 0.8

pow <- 0.8
gam <- log(pow, base = sig_lv)

# compute the FDR hat and bias
bias_list <- bias_list2 <- list()

for(i in 1:length(lam_vec)){
  fdr_vec[i] <- min(comp_mme(all_pval, lam_vec[i], sig_lv), 1) 
  bias_vec <- sapply(c(1:length(rho_vec)), 
                     FUN = function(k){bias_two(lam_vec[i], pi0, rho = rho_vec[k], gam, sig_lv)})
  bias_list[[i]] <- data.frame("rho" = rho_vec, "bias" = bias_vec, "lambda" = lam_vec[i])
}

# FDR hat in Table 2
round(data.frame("lambda" = lam_vec, "FDR_hat" = fdr_vec), 4)

# the last four rows of Table 2
round(bias_list[[1]][c(1,3,5,7),], 4)
round(bias_list[[2]][c(1,3,5,7),], 4)
round(bias_list[[3]][c(1,3,5,7),], 4)
round(bias_list[[4]][c(1,3,5,7),], 4)
round(bias_list[[5]][c(1,3,5,7),], 4)



