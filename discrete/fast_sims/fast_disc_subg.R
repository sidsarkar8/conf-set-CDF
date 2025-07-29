##### method
# sparse matrices
# pass warmstart info
# have a good start, and check only H_inv(1/lambda_H)

###### Dumbgen Wellner

library(confseq)
library(gurobi)
library(Matrix)
library(MASS)

# source("/Users/sidsarkar/Documents/Projects/conf_tail/bellmanford/get_bds_essHist.R")
# source("/Users/sidsarkar/Documents/Projects/conf_tail/bellmanford/get_bds_eH_mc.R")
source("~/Documents/Projects/Conf/dumbgun_wellner_new/get_bds_dw1.R")
source("~/Documents/Projects/Conf/dkw/get_bds_dkw.R")

#######################
######### xlogx

H = function(x){return( exp(x^2) ) }
H_inv = function(x){return( sqrt(log(x)) )}
H_der = function(x){return( (2*x)*exp(x^2) )}
H_der_inv = function(x){return( sqrt( lambertWp(x^2/2)/2 )  )}


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
  
  ans = cbind(ans_pos, ans_neg, H(x), rep(1, length(x)))
  
  return(ans)
}

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
  
  #### adding upper bound on objective
  # model_dual$A = rbind(model_dual$obj, model_dual$A)
  # model_dual$rhs = c(H_inv(K), model_dual$rhs)
  # model_dual$sense = c("<=", model_dual$sense)
  
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
    x_temp = c(x_temp,H_der_inv(1/lambda_H))
  }
  
  A_dual_temp = A_dual_dw(x = x_temp, data = data)
  slack_x_temp = (A_dual_temp%*%lambda) %>% as.vector() - x_temp
  
  #### if not to analyze
  if(!analyze){
    if(min(slack_x_temp) >= -10^(-6)){
      print("optimized!")
      return(NA)
    }else{
      return( x_temp[ which.min(slack_x_temp)])
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

################# spare + basis + ends_only ###################
# x1 = x2
# K1 = K2
# alpha = 0.1
# eps = 10^-7
# genCDF_n1 = genCDF_n2

# genCDF_n1  = data.frame(conf_low = rep(0, n2),
#                        conf_up = rep(1,n2))
# genCDF_n1 = data.frame(conf_low = rep(0,n2),
#                        conf_up = genCDF_n2$conf_up)
# genCDF_n1 = data.frame(conf_low = genCDF_n2$conf_low,
#                        conf_up = rep(1,n2))

up_bd_dw = function(x1, K1, alpha = 0.05, analyze = F,
                    eps = 10^-7, genCDF_n1)
{  
  
  n = length(x1)
  x1_sort = sort(x1)
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
    dual_temp = dual_lp_dw(x_start = x1_start, K = K1, data = x1, 
                           genCDF_n = genCDF_n1, 
                           vbasis = temp_vbasis,
                           cbasis = temp_cbasis)
    
    lambda_temp = dual_temp$x
    lambda_H_temp = lambda_temp[2*nrow(genCDF_n1) + 1]
    
    print(paste("objval:", round(dual_temp$objval,4)))
    print(paste("lambda_H:", lambda_H_temp))
    
    multi_factor = multi_factor + 1
    
    (x1_new = arg_min_end(lambda = lambda_temp,
                          data = x1,
                          genCDF_n = genCDF_n1,
                          low_bd = -multi_factor*H_inv(K1),
                          up_bd = multi_factor*H_inv(K1),
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

####################################################################################
####################################################################################
####################################################################################
################# int X_1 and X_n ###################

########################################################################
############### A_dual bdd ########################################


A_dual_dw_bdd = function(x, data)
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
  
  ans = cbind(ans_pos, ans_neg, rep(1, length(x)))
  
  return(ans)
}
########################################################################
############### dual sparse lp bdd #####################################

dual_lp_dw_bdd = function(x_start, data, genCDF_n, 
                          output_flag = 0)
{
  model_dual = list()
  
  model_dual$A = A_dual_dw_bdd( x = x_start, data = data)
  
  model_dual$obj  = c(genCDF_n$conf_up, -genCDF_n$conf_low, 1)
  
  model_dual$modelsense = "min"
  model_dual$rhs = x_start
  model_dual$sense =  rep(">=", length(x_start))
  
  model_dual$lb = c( rep(0, (2*nrow(genCDF_n))), -Inf)
  model_dual$ub = rep(Inf, (2*nrow(genCDF_n)) + 1)
  
  result_dual = gurobi(model_dual, list(OutputFlag = output_flag))
  
  return(result_dual)
}  

########################################################################
############### up bd algo #####################################


up_bd_dw_x1_to_xn = function(x1, alpha = 0.05, analyze = F,
                             genCDF_n1, eps = 10^-6)
{
  n = length(x1)
  x1_sort = sort(x1)
  
  #### starting optimization
  ### sorted data point only sufficient to check
  x1_start = c(x1_sort, x1_sort[1:(n-1)] + eps)
  
  dual_temp = dual_lp_dw_bdd(x_start = x1_start, data = x1, genCDF_n = genCDF_n1)
  
  return( dual_temp$objval)
}

#########################################################################
#########################################################################

n2 = 100
alpha = 0.1

x2 = rdistr(n2)
x2 = H_inv(x2^q1)*(2*rbinom(n2, size = 1, prob = 1/2) -1)
K2 = 2*q1 + 1

# x2 = rnorm(n2)
# K2 = 1

# x2 = 2*runif(n2) - 1
# K2 = 1

get_bds_CDF_n2 = get_bds_dw(n = n2, alpha = alpha)
genCDF_n2 = data.frame(conf_low = get_bds_CDF_n2$low,
                       conf_up = get_bds_CDF_n2$up)


{
  t_dw = Sys.time()
  up_dw = up_bd_dw(x1 = x2, K1 = K2, alpha = alpha, 
                   genCDF_n1 = genCDF_n2)
  low_dw = -up_bd_dw(x1 = -x2, K1 = K2, alpha = alpha, 
                     genCDF_n1 = genCDF_n2)
  
  
  t_dw = Sys.time() - t_dw
}

{
  t_dw_x1_to_xn = Sys.time()
  up_dw_x1_to_xn = up_bd_dw_x1_to_xn(x1 = x2,alpha = alpha, 
                                     genCDF_n1 = genCDF_n2)
  low_dw_x1_to_xn = -up_bd_dw_x1_to_xn(x1 = -x2, alpha = alpha, 
                                       genCDF_n1 = genCDF_n2)
  
  
  t_dw_x1_to_xn = Sys.time() - t_dw_x1_to_xn
}

up_tail = (sort(x2)[n2]/H(sort(x2)[n2]))*K2
low_tail = -(abs(sort(x2)[1])/H(sort(x2)[1]))*K2

CI_t = t.test(x2, conf.level = 1 - alpha)
CI_t = CI_t$conf.int[1:2]

#c(t_dw, t_dw_x1_to_xn)
c(low_dw, up_dw)
c(low_dw_x1_to_xn, up_dw_x1_to_xn)
c(low_dw_x1_to_xn + low_tail, 
  up_dw_x1_to_xn + up_tail)
CI_t
