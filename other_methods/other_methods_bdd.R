### implementation of other methods
## Hoeffding, Chebyshev, Anderson's method, 
## wald CI

###### Dumbgen Wellner

library(doParallel)
library(gurobi)
library(Matrix)
library(pracma)
library(MASS)
library(dplyr)

########################################################################
########################################################################
# 
# KL_dw = function(p,q){
#   if(length(p)!= length(q))
#   {
#     print("p and q length not same")
#     stop()
#   }
#   
#   val = numeric(length(p))
#   
#   for( i in 1:length(p))
#   {
#     val[i] = p[i]*log(p[i]/q[i]) + (1-p[i])*log((1-p[i])/(1-q[i]))
#     if( p[i] == 0 ){val[i] = log(1/(1-q[i]))}
#     if( p[i] == 1 ){val[i] = log(1/q[i]) }
#     if( q[i] == 0 ){val[i] = Inf }
#     if( q[i] == 1 ){val[i] = Inf}
#     
#     if( p[i] == 0 & q[i] == 0 ){val[i] = 0}
#     if( p[i] == 1 & q[i] == 0 ){val[i] = Inf }
#     if( p[i] == 0 & q[i] == 1 ){val[i] = Inf }
#     if( p[i] == 1 & q[i] == 1 ){val[i] = 0}
#   }
#   
#   return(val)
# }
# 
# #########################
# C_dw = function(t)
# {
#   return( log( 1 - log(2*t) - log(2*(1-t)) ) )
# }
# 
# #########################
# D_dw = function(t)
# {
#   return( log(1 + (C_dw(t))^2) )
# }
# 
# C_v_dw = function(u,t, v)
# {
#   len_t = length(t)
#   if(length(u)!= len_t)
#   {
#     print("not same length")
#     stop
#   }
#   
#   temp = numeric(len_t)
#   for( i in 1:len_t)
#   {
#     if(max(u[i],t[i]) < 1/2){
#       temp[i] = C_dw(max(u[i],t[i])) + (v*D_dw(max(u[i],t[i])))
#     } else if(min(u[i],t[i]) > 1/2){
#       temp[i] = C_dw(min(u[i],t[i])) + (v*D_dw(min(u[i],t[i])))
#     } else {
#       temp[i] = 0
#     }
#   }
#   return(temp)
# }
# 
# #########################
# 
# dw_quant = function( v = 3/2, alpha = 0.05, n_runs = 20000, n )
# {
#   sam = numeric(n_runs)
#   t_n = (1:n)/(n) 
#   t_n_shift = (0:(n-1))/n
#   for( i in 1:n_runs)
#   {
#     U_sort = sort( runif(n))
#     
#     sam[i] = max( n*KL_dw(t_n_shift, U_sort) - C_v_dw(t_n_shift, U_sort,v = v),  
#                   n*KL_dw(t_n, U_sort) - C_v_dw(t_n, U_sort,  v = v) ) 
#   }
#   
#   n_cuttoff = ceiling(n_runs*(1-alpha))
#   
#   return(sort(sam)[n_cuttoff])
# }
# 
# ########################################################################
# ########################################################################
# 
# 
# fin_func_dw = function(p,q,quant,v = 3/2,n)
# {
#   return(KL_dw(p,q) - ((C_v_dw(p,q,v = v)+ quant)/n)  )
# }
# 
# get_bds_dw = function(alpha = 0.05, n, v = 3/2, quant_val = NA)
# {
#   if(is.na(quant_val))
#   { 
#     dw_quant_val = dw_quant(alpha = alpha, n = n, v = v)
#   }else{
#     dw_quant_val = quant_val
#   }
#   
#   bds_dw = data.frame(low = rep(NA, n),
#                       up = rep(NA, n))
#   
#   t_n = (1:n)/n
#   #t_n_shift = (0:(n-1))/n
#   
#   for(i in 1:(n-1))
#   {
#     ### upper bound
#     if( fin_func_dw(p = t_n[i], q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
#       temp_up = 1
#     }else{
#       temp_up = uniroot( fin_func_dw, lower = t_n[i] , upper = 1-10^(-10),
#                          p = t_n[i], quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
#     }
#     
#     bds_dw$up[i]  = temp_up
#     bds_dw$low[n-i] = 1 - temp_up
#   }
#   
#   #### upper for i = n
#   bds_dw$up[n] = 1
#   
#   ### upper bound for i = 0
#   if( fin_func_dw(p = 0, q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
#     temp_up_0 = 1
#   }else{
#     temp_up_0 = uniroot( fin_func_dw, lower = 0 , upper = 1-10^(-10),
#                          p = 0, quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
#   }
#   
#   ##### upper for i = 0
#   bds_dw$low[n] = 1 - temp_up_0
#   
#   return(bds_dw)
# }


########### DKW Simple


get_bds_dkw_simp = function(alpha = 0.05, n)
{
  bds_dkw_simp = data.frame(low = rep(NA, n),
                            up = rep(NA, n))
  
  bds_dkw_simp$low = (1:n)/n - sqrt(( log(2) - log(alpha))/(2*n))
  bds_dkw_simp$up = (1:n)/n + sqrt(( log(2) - log(alpha))/(2*n))
  
  bds_dkw_simp$low[ bds_dkw_simp$low < 0] = 0
  bds_dkw_simp$up[ bds_dkw_simp$up > 1] = 1
  
  return(bds_dkw_simp)
}


################################################################################
##### Hoeffding CI

hf_CI = function(x, alpha, low_bd, up_bd)
{
  n = length(x)
  dev_mean = sqrt( log(2/alpha)/(2*n))*(up_bd - low_bd)
  
  return(mean(x) + c(-1,1)*dev_mean)
}

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

################################################################################
#####  Anderson's method

anderson_CI = function(x, alpha, low_bd, up_bd)
{
  n = length(x)
  x_sc = (x - low_bd)/(up_bd-low_bd) 
  
  bds_dw = get_bds_dkw_simp( n = n, alpha = alpha)  
  #bds_dw = get_bds_dw( n = n, alpha = alpha)  
  
  x_sc_ord = c(0, sort(x_sc), 1 )
  x_sc_ord_diff = (x_sc_ord - lag(x_sc_ord,1))[2:(n+2)]
  
  up_bd_wts = c(1, 1 - bds_dw$low)
  CI_up_bd = sum(up_bd_wts*x_sc_ord_diff)
  
  low_bd_wts = c(1-bds_dw$up,0)
  CI_low_bd = sum(low_bd_wts*x_sc_ord_diff)
  
  return(  c(CI_low_bd,CI_up_bd)*(up_bd-low_bd) + low_bd  )
}
################################################################################
##### wald CI

wald_CI = function(x, alpha)
{
  return(t.test(x, conf.level = 1 - alpha)$conf.int[1:2])
}

################################################################################
##### Empirical Bernstein CI of Maurer and Pontil

emp_bern_CI = function(x, alpha, low_bd, up_bd)
{
  n = length(x)
  
  x_sc = (x - low_bd)/(up_bd-low_bd) 
  
  var_hat = var(x_sc)
  
  dev_mean = sqrt( (2*var_hat*log(4/alpha))/n )
  add_term = (7*log(4/alpha))/(3*(n-1))
  
  CI_low_bd = mean(x_sc)  - add_term - dev_mean
  CI_up_bd = mean(x_sc) + add_term + dev_mean
  
  return(  c(CI_low_bd,CI_up_bd)*(up_bd-low_bd) + low_bd  )
  
}

################################################################################
##### Plug-in empirical bernstein CI

PI_emp_bern_CI = function(x, alpha, low_bd, up_bd, c_bd = 1/2)
{
  n = length(x)
  
  x_sc = (x - low_bd)/(up_bd-low_bd) 
  
  var_seq = numeric(n)
  mean_seq = numeric(n)
  
  lambda_seq_eB = numeric(n)
  v_seq = numeric(n)
  
  ### first term in lambda sequence, since sigma_0 not definined.
  
  lambda_seq_eB[1] = 3/4
  v_seq[1] = 0
  
  for( i in 1:n)
  {
    mean_seq[i] = (1/2 + sum(x_sc[1:i]))/(i+1)
    var_seq[i] = (1/4 + sum( (x_sc[1:i] - mean_seq[1:i])^2 ))/(i+1)
    
    
    if(i > 1)
    {
      lambda_seq_eB[i] = min(sqrt((2*log(2/alpha))/(n*var_seq[i-1])), c_bd)
      v_seq[i] = 4*(x_sc[i]- mean_seq[i-1])^2
    }
    
  }
  
  
  mean_wtd = sum(lambda_seq_eB*x_sc)/sum( lambda_seq_eB)
  dev_mean = (log(2/alpha) + sum(v_seq*(-log(1 - lambda_seq_eB) - lambda_seq_eB)/4))/sum(lambda_seq_eB)
  
  CI_up_bd = (mean_wtd + dev_mean)
  CI_low_bd =(mean_wtd - dev_mean)
  
  return( c(CI_low_bd,CI_up_bd)*(up_bd - low_bd) + low_bd )    
}


################################################################################
##### Plug-in Hedged CI

K_plus = function(x, m, lambda_seq, c_bd = 1/2){
  
  ans = numeric(length(m))
  
  for( i in 1:length(m) ){
    ans[i] = prod( 1 + pmin(abs(lambda_seq), c_bd/m[i])*(x - m[i]))
  }
  
  return(ans)
}

K_minus = function(x, m, lambda_seq, c_bd  = 1/2){
  
  ans = numeric(length(m))
  
  for( i in 1:length(m)){
    ans[i] = prod( 1 - pmin(abs(lambda_seq), c_bd/(1-m[i]))*(x - m[i])) 
  }
  return(ans)
}

PI_hedged_CI = function(x, alpha, low_bd, up_bd, c_bd = 1/2, theta = 1/2, grid_precision = 0.01)
{
  
  n = length(x)
  
  x_sc = (x - low_bd)/(up_bd-low_bd) 
  
  var_seq = numeric(n)
  mean_seq = numeric(n)
  
  lambda_seq_H = numeric(n)
  v_seq = numeric(n)
  
  ### first term in lambda sequence, since sigma_0 not definined.
  
  lambda_seq_H[1] = sqrt((2*log(2/alpha))/(n*0.25))
  v_seq[1] = 0
  
  for( i in 1:n)
  {
    mean_seq[i] = (1/2 + sum(x_sc[1:i]))/(i+1)
    var_seq[i] = (1/4 + sum( (x_sc[1:i] - mean_seq[1:i])^2 ))/(i+1)
    
    if(i > 1)
    {
      lambda_seq_H[i] =sqrt((2*log(2/alpha))/(n*var_seq[i-1]))
      v_seq[i] = 4*(x_sc[i]- mean_seq[i-1])^2
    }
    
  }
  
  m_grid = seq(from = 0, to =  1, by = grid_precision)
  
  ans_plus = K_plus(x_sc, m_grid, lambda_seq_H)*(theta)
  ans_minus = K_minus(x_sc, m_grid, lambda_seq_H)*(1-theta)
  
  # plot(x = m_grid, y = ans_plus, type = "l", ylim = c(0,100))
  # points(x = m_grid, y = ans_minus, type = "l", col = "red")
  # abline(h = 1/alpha, col = "green")
  
  CI_up_bd = m_grid[max(which(ans_minus < 1/alpha))]
  CI_low_bd = m_grid[min(which(ans_plus < 1/alpha))]
  
  return(c(CI_low_bd, CI_up_bd)*(up_bd-low_bd) + low_bd)
}

##################################################################
##################################################################

# n = 100
# a1 = 2
# b1 = 40
# lb1 = 0
# ub1 = 1
# al1 = 0.05
# 
# 
# n_sims = 100
# cov_sims = data.frame("and" = rep(NA,n),
#                       "cheby" = rep(NA,n),
#                       "emp_bern" = rep(NA,n),
#                       "hf" = rep(NA,n),
#                       "PI_emp_bern" = rep(NA,n),
#                       "PI_hedged" = rep(NA,n),
#                       "wald" = rep(NA,n))
# 
# width_sims = data.frame("and" = rep(NA,n),
#                         "cheby" = rep(NA,n),
#                         "emp_bern" = rep(NA,n),
#                         "hf" = rep(NA,n),
#                         "PI_emp_bern" = rep(NA,n),
#                         "PI_hedged" = rep(NA,n),
#                         "wald" = rep(NA,n))
# 
# mu1 =   a1/(a1 + b1)
# 
# for(sim_idx in 1:n_sims)
# {
#   print(sim_idx)
#   x1 = rbeta(n, shape1 = a1, shape2 = b1 )
# 
#   CI_and = anderson_CI(x1,al1,lb1,ub1)
#   cov_sims$and[sim_idx] = ifelse(mu1 <= CI_and[2] & mu1 >= CI_and[1],1,0)
#   width_sims$and[sim_idx] = CI_and[2] - CI_and[1]
# 
#   CI_cheby = cheby_CI(x1,al1,1)
#   cov_sims$cheby[sim_idx] = ifelse(mu1 <= CI_cheby[2] & mu1 >= CI_cheby[1],1,0)
#   width_sims$cheby[sim_idx] = CI_cheby[2] - CI_cheby[1]
# 
# 
#   CI_emp_bern = emp_bern_CI(x1,al1,lb1,ub1)
#   cov_sims$emp_bern[sim_idx] = ifelse(mu1 <= CI_emp_bern[2] & mu1 >= CI_emp_bern[1],1,0)
#   width_sims$emp_bern[sim_idx] = CI_emp_bern[2] - CI_emp_bern[1]
# 
# 
#   CI_hf = hf_CI(x1,al1,lb1,ub1)
#   cov_sims$hf [sim_idx] = ifelse(mu1 <= CI_hf[2] & mu1 >= CI_hf[1],1,0)
#   width_sims$hf[sim_idx] = CI_hf[2] - CI_hf[1]
# 
#   CI_PI_emp_bern = PI_emp_bern_CI(x1,al1,lb1,ub1)
#   cov_sims$PI_emp_bern [sim_idx] = ifelse(mu1 <= CI_PI_emp_bern[2] & mu1 >= CI_PI_emp_bern[1],1,0)
#   width_sims$PI_emp_bern[sim_idx] = CI_PI_emp_bern[2] - CI_PI_emp_bern[1]
# 
#   CI_PI_hedged = PI_hedged_CI(x1,al1,lb1,ub1,0.001)
#   cov_sims$PI_hedged [sim_idx] = ifelse(mu1 <= CI_PI_hedged[2] & mu1 >= CI_PI_hedged[1],1,0)
#   width_sims$PI_hedged[sim_idx] = CI_PI_hedged[2] - CI_PI_hedged[1]
# 
#   CI_wald = wald_CI(x1,al1)
#   cov_sims$wald [sim_idx] = ifelse(mu1 <= CI_wald[2] & mu1 >= CI_wald[1],1,0)
#   width_sims$wald[sim_idx] = CI_wald[2] - CI_wald[1]
# }
# 
# 
# apply(cov_sims,2,mean)
# apply(width_sims,2,mean)
