
################################################################################
##### chebysev inequality

cheby_CI = function(x, alpha, K1)
{
  n = length(x)
  
  cntr = (mean(x)*alpha)/(alpha + (1/n))
  dev_cntr = (alpha*(K1 - mean(x)^2))/n + K1/n^2
  dev_cntr = sqrt( dev_cntr )/(alpha + (1/n))
  
  return(cntr + c(-1,1)*dev_cntr )
}


cheby_centered_CI = function(x, alpha, K1)
{
  n = length(x)
  
  cntr = mean(x)
  dev_cntr = sqrt(K1/(n*alpha))
  
  return(cntr + c(-1,1)*dev_cntr )
}


################################################################################
##### wald CI

wald_CI = function(x, alpha)
{
  return(t.test(x, conf.level = 1 - alpha)$conf.int[1:2])
}


################################################################################
##### 

phi_wide <- function(x) {
  
  temp = numeric(length(x))
  for( i in 1:length(x))
  {
    if (x[i] > 0) {
      temp[i] = (log(1 + x[i] + abs(x[i])^2 / 2))
    } else {
      temp[i] = (-log(1 - x[i] + abs(x[i])^2 / 2))
    }
  }
  return(temp)
}

score_phi = function(lambda, x, theta)
{
  return( mean(phi_wide(lambda*(x-theta))) )  
}


catoni_CI = function(x,alpha,v)
{
  n = length(x)
  if(n < log(2/alpha)){stop("Too small a sample size")}
  
  eta = sqrt((2*v*log(2/alpha))/(n - (2*log(2/alpha))))
  lambda = sqrt( (2*log(2/alpha))/(n*(v + eta^2)) )
  
  ctr = (uniroot( score_phi, lower = min(x) , upper = max(x),
                  lambda = lambda, x = x, tol = 10^-10 )$root)
  
  return( c(-eta, eta) + ctr )
}


#################################################################
#################################################################
# 
# n = 100
# al1 = 0.05
# 
# 
# n_sims = 100
# 
# cov_sims = data.frame("cheby" = rep(NA,n),
#                       "catoni" = rep(NA,n),
#                       "wald" = rep(NA,n))
# 
# width_sims = data.frame("cheby" = rep(NA,n),
#                         "catoni" = rep(NA,n),
#                         "wald" = rep(NA,n))
#
#
# mu1 = 0
# 
# for(sim_idx in 1:n_sims)
# {
#   print(sim_idx)
#   x1 = rnorm(n, mu1, sd = 1)
# 
#   CI_cheby = cheby_CI(x1,al1, mu1^2 + 1)
#   cov_sims$cheby[sim_idx] = ifelse(mu1 <= CI_cheby[2] & mu1 >= CI_cheby[1],1,0)
#   width_sims$cheby[sim_idx] = CI_cheby[2] - CI_cheby[1]
# 
# 
#   CI_catoni = catoni_CI(x1,al1,mu1^2 + 1)
#   cov_sims$catoni [sim_idx] = ifelse(mu1 <= CI_catoni[2] & mu1 >= CI_catoni[1],1,0)
#   width_sims$catoni[sim_idx] = CI_catoni[2] - CI_catoni[1]
# 
#   CI_wald = wald_CI(x1,al1)
#   cov_sims$wald [sim_idx] = ifelse(mu1 <= CI_wald[2] & mu1 >= CI_wald[1],1,0)
#   width_sims$wald[sim_idx] = CI_wald[2] - CI_wald[1]
# }
# 
# 
# apply(cov_sims,2,mean)
# apply(width_sims,2,mean)
