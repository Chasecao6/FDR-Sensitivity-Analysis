require(zipfR)

########################################
#### 1. Function to compute the MME ####
########################################

# Inputs
# all_pval: the vector of all published p-values less than alpha
# lam: a cutoff to be used to compute the MME
# alp: significance level

# Output
# the resulting MME for pi0

comp_mme <- function(all_pval, lam, alp = 0.05){
  signif_pval <- all_pval[all_pval<=alp]
  adj_pval <- signif_pval/alp # as in Hung and Fithian (2020)
  mme <- sum(adj_pval > lam/alp)/(length(adj_pval)*(1 - lam/alp))
  return(mme)
}

############################################################################################################
##### 2. Functions to compute the bias of MME under each SPM with F_{1} = t^{r} under true alternative #####
############################################################################################################

###############################################################
#### 2.1 function to compute the bias under the simple SPM ####
###############################################################

# Inputs
# lam: lambda for MME
# pi0: null proportion
# gam: gamma appearing in the dist'n under true alternative
# alp: significance level

# output
# bias

bias_simple <- function(lam, pi0, gam = log(0.8, base = 0.05), alp = 0.05){
  pi1 <- 1 - pi0
  first <- (pi0*(alp - lam) + pi1*(alp^gam - lam^gam))/(pi0*alp + pi1*alp^gam)
  second <- (pi0*alp)/(pi0*alp + pi1*alp^gam)
  out <- (alp/(alp - lam))*first - second
  return(out)
  #return(alp*pi1*(alp^gam - lam^gam)/((alp - lam)*(pi0*alp + pi1*alp^gam)))
  #return((alp/(alp - lam))*(pi1*(alp^gam - lam^gam)/(pi0*alp + pi1*alp^gam)))
}


############################################################################
#### 2.2 function to compute the bias under the SPM with two components ####
############################################################################

# Inputs
# lam: lambda for MME
# pi0: null proportion
# rho: non-negative constant associated with the second indicator function I_{-}
# gam: gamma appearing in the dist'n under true alternative
# alp: significance level

# output
# bias

bias_two <- function(lam, pi0, rho, gam = log(0.8, base = 0.05), alp = 0.05){
  pi1 <- 1 - pi0
  if(lam >= alp/2){
    # numerator
    numer <- (alp/(alp - lam))*(pi0*rho*(alp - lam) + pi1*rho*(alp^gam - lam^gam)) - pi0*(1 + rho)*alp/2 
    # denominator
    denom <- pi0*(1 + rho)*alp/2 + pi1*((alp/2)^gam + rho*(alp^gam - (alp/2)^gam))
  } else{
    # numerator
    numer <- (alp/(alp - lam))*(pi0*(rho*alp/2 + alp/2 - lam) + pi1*(rho*(alp^gam - (alp/2)^gam) + (alp/2)^gam - lam^gam)) - pi0*(1 + rho)*alp/2
    # denominator
    denom <- pi0*(1 + rho)*alp/2 + pi1*((alp/2)^gam + rho*(alp^gam - (alp/2)^gam))
  }
  return(numer/denom)
}


####################################################################
#### 2.3 function to compute the bias under the Beta kernel SPM ####
####################################################################

# Inputs
# lam: lambda for MME
# pi0: null proportion
# b: the second shape parameter for the Betal kernel SPM which is expected to be greater or equal to 1
# gam: gamma appearing in the dist'n under true alternative
# alp: significance level

# output
# bias

# ib <- zipfR::Ibeta(0.05, a = 0.07, 2, lower = TRUE)

bias_beta <- function(lam, pi0, b, gam = log(0.8, base = 0.05), alp = 0.05){
  # first compute the incomplete beta function parts
  ib1 <- zipfR::Ibeta(alp, a = gam, b, lower = TRUE)
  ib2 <- zipfR::Ibeta(lam, a = gam, b, lower = TRUE)
  pi1 <- 1 - pi0
  # denominator
  denom <- (pi0/b)*(1 - (1 - alp)^b) + pi1*gam*ib1
  # first term of numerator
  numer1 <- (alp/(alp - lam))*(pi0*((1 - lam)^b - (1 - alp)^b)/b + pi1*gam*(ib1 - ib2))
  # second term of numerator
  numer2 <- (pi0/b)*(1 - (1 - alp)^b)
  # numerator
  numer <- numer1 - numer2
  return(numer/denom)
}


################################################################################################################
##### 3. Functions to compute the variance of MME under each SPM with F_{1} = t^{r} under true alternative #####
################################################################################################################


##########################################################################
#### 3.1 function to compute the bias under the SPM of step functions ####
##########################################################################


var_step <- function(lam, m, pi, rho, gam = log(0.8, base = sig_lv), sig_lv = 0.05){
  
  E_mu <- pi*(rho + 1)*sig_lv/2 + (1 - pi)*(rho*sig_lv^gam + (1 - rho)*(sig_lv/2)^gam)
  
  if(lam >= sig_lv/2){
    q1 <- (pi*(rho*lam + (1 - rho)*sig_lv/2) + (1 - pi)*(rho*sig_lv^gam + (1 - rho)*(sig_lv/2)^gam))/E_mu
    q2 <- (pi*rho*(sig_lv - lam) + (1 - pi)*rho*(sig_lv^gam - lam^gam))/E_mu  
  } else{
    q1 <- (pi*lam + (1 - pi)*lam^gam)/E_mu
    q2 <- (pi*((1 + rho)*sig_lv/2 - lam) + (1 - pi)*(rho*sig_lv^gam + (1 - rho)*(sig_lv/2)^gam - lam^gam))/E_mu  
  }
 
  out <- (sig_lv/(sig_lv - lam))^2*((q1*q2)/(m*E_mu*(q1 + q2)^3))
  return(out) 
}


#######################################################################
#### 3.2 function to compute the bias under the SPM of beta kernel ####
#######################################################################

var_beta <- function(lam, m, pi, b, gam = log(0.8, base = sig_lv), sig_lv = 0.05){
  #ib1 <- zipfR::Ibeta(alp, a = gam, b, lower = TRUE)
  E_mu <- pi/b + (1 - pi)*gam*zipfR::Ibeta(1, a = gam, b, lower = TRUE)
  beta_alp <- zipfR::Ibeta(sig_lv, a = gam, b, lower = TRUE)
  beta_lam <- zipfR::Ibeta(lam, a = gam, b, lower = TRUE)
  q1 <- ((pi/b)*(1 - (1 - lam)^b) + (1 - pi)*gam*beta_lam)/E_mu
  q2 <- ((pi/b)*((1 - lam)^b - (1 - sig_lv)^b) + (1 - pi)*gam*(beta_alp - beta_lam))/E_mu
  out <- (sig_lv/(sig_lv - lam))^2*((q1*q2)/(m*E_mu*(q1 + q2)^3))
  return(out) 
}


