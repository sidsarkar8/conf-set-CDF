##### centred 2-moments (aka variance)

###### Dumbgen Wellner

library(doParallel)
library(gurobi)
library(Matrix)
library(pracma)
library(MASS)
library(dplyr)

########################################################################
########################################################################

source("/home/siddhaarth/conf_tail/sims_new/other_methods_x2.R")

#########################

KL_dw = function(p,q){
  if(length(p)!= length(q))
  {
    print("p and q length not same")
    stop()
  }
  
  val = numeric(length(p))
  
  for( i in 1:length(p))
  {
    val[i] = p[i]*log(p[i]/q[i]) + (1-p[i])*log((1-p[i])/(1-q[i]))
    if( p[i] == 0 ){val[i] = log(1/(1-q[i]))}
    if( p[i] == 1 ){val[i] = log(1/q[i]) }
    if( q[i] == 0 ){val[i] = Inf }
    if( q[i] == 1 ){val[i] = Inf}
    
    if( p[i] == 0 & q[i] == 0 ){val[i] = 0}
    if( p[i] == 1 & q[i] == 0 ){val[i] = Inf }
    if( p[i] == 0 & q[i] == 1 ){val[i] = Inf }
    if( p[i] == 1 & q[i] == 1 ){val[i] = 0}
  }
  
  return(val)
}

#########################
C_dw = function(t)
{
  return( log( 1 - log(2*t) - log(2*(1-t)) ) )
}

#########################
D_dw = function(t)
{
  return( log(1 + (C_dw(t))^2) )
}

C_v_dw = function(u,t, v)
{
  len_t = length(t)
  if(length(u)!= len_t)
  {
    print("not same length")
    stop
  }
  
  temp = numeric(len_t)
  for( i in 1:len_t)
  {
    if(max(u[i],t[i]) < 1/2){
      temp[i] = C_dw(max(u[i],t[i])) + (v*D_dw(max(u[i],t[i])))
    } else if(min(u[i],t[i]) > 1/2){
      temp[i] = C_dw(min(u[i],t[i])) + (v*D_dw(min(u[i],t[i])))
    } else {
      temp[i] = 0
    }
  }
  return(temp)
}

#########################

dw_quant = function( v = 3/2, alpha = 0.05, n_runs = 20000, n )
{
  sam = numeric(n_runs)
  t_n = (1:n)/(n) 
  t_n_shift = (0:(n-1))/n
  for( i in 1:n_runs)
  {
    U_sort = sort( runif(n))
    
    sam[i] = max( n*KL_dw(t_n_shift, U_sort) - C_v_dw(t_n_shift, U_sort,v = v),  
                  n*KL_dw(t_n, U_sort) - C_v_dw(t_n, U_sort,  v = v) ) 
  }
  
  n_cuttoff = ceiling(n_runs*(1-alpha))
  
  return(sort(sam)[n_cuttoff])
}

########################################################################
########################################################################


fin_func_dw = function(p,q,quant,v = 3/2,n)
{
  return(KL_dw(p,q) - ((C_v_dw(p,q,v = v)+ quant)/n)  )
}

get_bds_dw = function(alpha = 0.05, n, v = 3/2, quant_val = NA)
{
  if(is.na(quant_val))
  { 
    dw_quant_val = dw_quant(alpha = alpha, n = n, v = v)
  }else{
    dw_quant_val = quant_val
  }
  
  bds_dw = data.frame(low = rep(NA, n),
                      up = rep(NA, n))
  
  t_n = (1:n)/n
  #t_n_shift = (0:(n-1))/n
  
  for(i in 1:(n-1))
  {
    ### upper bound
    if( fin_func_dw(p = t_n[i], q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
      temp_up = 1
    }else{
      temp_up = uniroot( fin_func_dw, lower = t_n[i] , upper = 1-10^(-10),
                         p = t_n[i], quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
    }
    
    bds_dw$up[i]  = temp_up
    bds_dw$low[n-i] = 1 - temp_up
  }
  
  #### upper for i = n
  bds_dw$up[n] = 1
  
  ### upper bound for i = 0
  if( fin_func_dw(p = 0, q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
    temp_up_0 = 1
  }else{
    temp_up_0 = uniroot( fin_func_dw, lower = 0 , upper = 1-10^(-10),
                         p = 0, quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
  }
  
  ##### upper for i = 0
  bds_dw$low[n] = 1 - temp_up_0
  
  return(bds_dw)
}


####### tight distr ###########

q1 = 2

F1 = function(x){
  temp = rep(NA,length(x))
  
  for( i in 1: length(x))
  {
    temp[i] = 1/(abs(x[i])^q1*(max(log(abs(x[i])),1))^2)
  }
  return(temp)
}

F1_inv = function(x){
  
  temp = rep(NA,length(x))
  
  for( i in 1: length(x))
  {
    if( x[i] < -exp(1)){
      temp[i] = ( -exp(lambertWp( q1/(2*sqrt(x[i]) )*(2/q1)) ) )
    }else{
      temp[i] = ( -(x[i])^(-1/q1) )     
    }
  }
  return(temp)
}

rdistr = function(n)
{
  toss = 2*rbinom(n,1, 0.5) - 1
  U = runif(n)
  return(toss*F1_inv(U))
}

F_distr = function(x)
{
  temp = numeric(length(x))
  
  for(i in 1:length(x))
  {
    if(x[i]<=-1){
      temp[i] = F1(x[i])/2
    }else if(x[i]>1){
      temp[i] =1 - F1(x[i])/2
    }else{
      temp[i] = 1/2
    }
  }
  return(temp)
}


####################################################################################
################# int lb and ub ###################

H = function(x, mu){return((x-mu)^2)}
## returns both inverses
H_inv = function(x, mu){return( c(mu - sqrt(x), mu + sqrt(x)) )}
H_der = function(x, mu){return(2*(x - mu))}
H_der_inv = function(x, mu){return(x/2 + mu)}


#############################################################################
################ A_dual across methods ######################################
#############################################################################

########################################################################
############### A_dual sparse ########################################


A_dual_dw = function(x, data)
{
  ans_pos = NULL
  ans_neg = NULL
  
  data_sort = sort(data)
  n_constr = length(data)
  
  for( i in 1:length(x))
  {
    #### which datapoints are 
    temp_idx = which(data_sort >= x[i])
    
    temp_spMat = spMatrix(nrow = 1, ncol = n_constr,
                          i = rep(1,length(temp_idx)),
                          j = temp_idx, 
                          x = rep(1,length(temp_idx)))
    
    ans_pos = rbind(ans_pos, temp_spMat)
    
    ans_neg = rbind(ans_neg, -temp_spMat)
  }
  
  ans = cbind(ans_pos, ans_neg, H(x, mean(data)), rep(1, length(x)))
  
  return(ans)
}

#############################################################################
################ dual algo across methods ######################################
#############################################################################

########################################################################
############### dual lp sparse ##############################


dual_lp_dw = function(x_start, K, data, genCDF_n, vbasis,
                      cbasis)
{
  model_dual = list()
  
  model_dual$A = A_dual_dw( x = x_start, data = data)
  
  model_dual$obj  = c(genCDF_n$conf_up, -genCDF_n$conf_low, K, 1)
  model_dual$modelsense = "min"
  model_dual$rhs = x_start
  model_dual$sense =  rep(">=", length(x_start))
  
  model_dual$vbasis = vbasis
  model_dual$cbasis = cbasis
  model_dual$LPWarmStart = 1
  model_dual$Method = 2
  
  model_dual$lb = c(rep(0, (2*nrow(genCDF_n)) + 1), -Inf)
  model_dual$ub = rep(Inf, (2*nrow(genCDF_n)) + 2)
  
  result_dual = gurobi(model_dual, list(OutputFlag = 0))
  
  return(result_dual)
}  

########################################################################
######## arg_min ends only  for of functions with [low_bd, up_bd] #######

arg_min_end = function(lambda, data, genCDF_n , eps = 10^(-7), analyze = F, 
                       low_bd, up_bd)
{
  x_temp = NULL
  
  ##### end points of the interval
  
  ##### end points of our serach query
  x_temp = c(x_temp, low_bd, up_bd)
  
  ##### if lambda_H!= 0, we have another candidate to consider
  
  lambda_H = lambda[2*nrow(genCDF_n) + 1]
  
  if(lambda_H != 0){
    x_temp = c(x_temp,H_der_inv(1/lambda_H, mean(data)))
  }
  
  A_dual_temp = A_dual_dw(x = x_temp, data = data)
  slack_x_temp = (A_dual_temp%*%lambda) %>% as.vector() - x_temp
  
  #### if not to analyze
  if(!analyze){
    if(min(slack_x_temp) >= -10^(-6)){
      print("optimized!")
      return(NA)
    }else{
      return( x_temp[which.min(slack_x_temp)])
    }
  }
  
  #### if to analyze
  if(analyze)
  {
    return(data.frame(x_min = x_temp[ which.min(slack_x_temp)], 
                      min = min(slack_x_temp), 
                      idx = which.min(slack_x_temp), 
                      lambda_H = lambda_H))
  }
}

#############################################################################
################  Algos ######################################
#############################################################################
# 
# x1 = x2
# K1 = K2 
# alpha = 0.1
# genCDF_n1 = genCDF_n2

up_bd_dw_centered = function(x1, K1, alpha = 0.05, analyze = F,
                             eps = 10^-7, genCDF_n1)
{  
  n = length(x1)
  x1_sort = sort(x1)
  
  delta = alpha/log(n)
  alpha_corr = alpha - delta
  # chebysev correction
  K1_corr = K1 + sqrt(K1/(n*delta))
  
  # alpha_corr = alpha
  # K1_corr = K1
  
  #### starting optimization
  
  x1_start = x1_sort[1]
  for(i in 1:(n-1) )
  {
    x1_start = c(x1_start, x1_sort[i] + eps, x1_sort[i+1]) 
  }
  
  
  iter = 0
  multi_factor = 0
  
  temp_vbasis = NULL
  temp_cbasis = NULL
  
  while( TRUE)
  {
    iter = iter + 1
    dual_temp = dual_lp_dw(x_start = x1_start, K = K1_corr, data = x1, 
                           genCDF_n = genCDF_n1, 
                           vbasis = temp_vbasis,
                           cbasis = temp_cbasis)
    
    
    if(dual_temp$status == "INF_OR_UNBD"){
      return(mean(x1) + sqrt(K1/(n*alpha)))
    }
    
    lambda_temp = dual_temp$x
    lambda_H_temp = lambda_temp[2*nrow(genCDF_n1) + 1]
    
    print(paste("objval:", round(dual_temp$objval,4)))
    print(paste("lambda_H:", lambda_H_temp))
    
    multi_factor = multi_factor + 1
    
    ## expanding the search region 
    # low = mu + multi_factor*(H_inv_low(mu)-mu)
    # up = mu + multi_factor*(H_inv_up(mu)-mu) 
    ## this keeps expanding as multi_factor grows
    
    (x1_new = arg_min_end(lambda = lambda_temp,
                          data = x1,
                          genCDF_n = genCDF_n1,
                          low_bd = multi_factor*(H_inv(K1, mean(x1))[1] - mean(x1)) + mean(x1),
                          up_bd = multi_factor*(H_inv(K1, mean(x1))[2] - mean(x1)) + mean(x1),
                          analyze = F))
    
    print(paste("new point", round(x1_new,4)))
    
    if(is.na(x1_new)){break}
    
    x1_start = c(x1_start, x1_new)
    
    temp_vbasis = dual_temp$vbasis
    temp_cbasis = c(dual_temp$cbasis, 0)
    
    print( paste("iter =", iter) )
    #points(x = iter, y = lambda_H_temp)
  }
  
  if(!analyze)
  { 
    return(dual_temp$objval)
  }else{
    ## TODO
    #primal_temp = primal_lp(x_start = x1_start, K = K1, data = x1, 
    #                        genCDF_n = genCDF_n1)
    #return(list(
    #  lambda = lambda_temp, x1_opt = x1_start, 
    #  p = primal_temp$x, objval = dual_temp$objval
    #))
  }
  
}


############################################################
############################################################
#######################################################

parallel_fin_func = function(simul_idx, n2, genCDF_n2, mu2)
{
  set.seed(simul_idx + 2025)
  
  x2 = rdistr(n2)
  x2 = H_inv(x2^q1,0)[1]*(2*rbinom(n2, size = 1, prob = 1/2) -1) + mu2
  K2 = 2*q1 + 1
  
  up_dw = up_bd_dw_centered(x1 = x2, K1 = K2, alpha = alpha, 
                            genCDF_n1 = genCDF_n2)
  low_dw = -up_bd_dw_centered(x1 = -x2, K1 = K2, alpha = alpha, 
                              genCDF_n1 = genCDF_n2)
  
  CI_t = wald_CI(x2, alpha)
  
  CI_cat = catoni_CI(x2,alpha, K2)
  
  CI_cheby = cheby_centered_CI(x2, alpha, K2) 
  
  qq1 = data.frame(
    "simul_idx"   = simul_idx,
    "up_dw"       = up_dw,
    "low_dw"      = low_dw,
    "up_t_test"   = CI_t[2],
    "low_t_test"  = CI_t[1],
    "up_catoni"   = CI_cat[2],
    "low_catoni"  = CI_cat[1],
    "up_cheby"    = CI_cheby[2],
    "low_cheby"   = CI_cheby[1]
  )
  
  write.csv(qq1,
            file = paste("/home/siddhaarth/conf_tail/sims_new/x2/indiv_files/x2_",
                         n2_seq[n2_idx],"_",simul_idx,".csv",sep = "" ))
  
  return(qq1)
  
}


n2_seq = c(20,50,100,250,500)

M = 500

numCores = detectCores() - 1
registerDoParallel(numCores) 

for( n2_idx in 1:length(n2_seq))
{
  print(paste(n2_idx, ": run sims"))
  
  alpha = 0.1
  mu2 = 3
  #delta = min(log( n2_seq[n2_idx])/ n2_seq[n2_idx], alpha/2)
  delta = alpha/log(n2_seq[n2_idx])
  alpha_corr = alpha - delta

  get_bds_CDF_n2 = get_bds_dw(n = n2_seq[n2_idx],
                              alpha = alpha_corr)
  genCDF_n2 = data.frame(conf_low = get_bds_CDF_n2$low,
                         conf_up = get_bds_CDF_n2$up)
  
  
  df_n2_x2 = foreach( simul_idx = 1:M,  .combine = rbind ) %dopar%{
    parallel_fin_func(simul_idx, 
                      n2 = n2_seq[n2_idx],
                      genCDF_n2 = genCDF_n2,
                      mu2 = mu2)
  }
  
  
  write.csv(df_n2_x2,
            file = paste("/home/siddhaarth/conf_tail/sims_new/x2/x2_",
                         n2_seq[n2_idx],".csv",sep = "" ))
}

stopImplicitCluster()
