##### trying to just add points, discretization

source("/Users/sidsarkar/Documents/Projects/conf_tail/bellmanford/get_bds_essHist.R")
source("/Users/sidsarkar/Documents/Projects/conf_tail/bellmanford/get_bds_eH_mc.R")

library(devtools)
library(confseq)
library(CVXR)
library(pracma)
library(gurobi)
library(Matrix)
library(MASS)

######### x^2

# H = function(x){return(x^2)}
# H_inv = function(x){return(sqrt(x))}
# H_der = function(x){return(2*x)}
# H_der_inv = function(x){return(x/2)}

######## x^1.5

# H = function(x){return(abs(x)^(3/2))}
# H_inv = function(x){return(abs(x)^(2/3))}
# H_der = function(x){return( (3/2)*sign(x)*(abs(x)^(1/2)) ) }
# H_der_inv = function(x){return( (4/9)*sign(x)*(x^2)  )}

######### exp(x^2/2)

# H = function(x){return(exp(x^2/2))}
# H_inv = function(x){return(sqrt(2*log(x)))}
# H_der = function(x){return(exp(x^2/2)*(x))}
# H_der_inv = function(x){return( sqrt(lambertWp(x^2))  )}

### linear translation

# b1 = 1.4
# a1 = 1
# H = function(x){return(x^2 - 2*b1*x)}
# H_inv = function(x){return(sqrt(x + b1^2) + b1)}
# H_der = function(x){return(2*x - 2*b1)}
# H_der_inv = function(x){return((x/2) + b1)}

######### prior functions needed for LP construction and optimisation ####
###### function for A

A_dual = function(x, data, genIntv_n)
{
  ans = matrix(ncol = (2*nrow(genIntv_n)) + 2, nrow = length(x))
  
  data_sort = sort(data)
  
  for( i in 1:length(x))
  {
    for(j in 1:nrow(genIntv_n))
    {
      ans[i,j] = ifelse( data_sort[genIntv_n$left[j]]< x[i] && data_sort[genIntv_n$right[j]] >= x[i],
                         1, 0)
      ans[i, j+nrow(genIntv_n)] = (-1)*ans[i,j]
    }
    
    ans[i, 2*nrow(genIntv_n) + 1] = H(x[i])
    ans[i, 2*nrow(genIntv_n) + 2] = 1
  }
  return(ans)
}

############### dual lp

dual_lp = function(x_start, K, data, genIntv_n)
{
  model_dual = list()
  
  model_dual$A = A_dual( x = x_start, data = data, genIntv_n = genIntv_n)
  
  model_dual$obj  = c(genIntv_n$conf_up, -genIntv_n$conf_low, K, 1)
  model_dual$modelsense = "min"
  model_dual$rhs = x_start
  model_dual$sense =  rep(">=", length(x_start))
  
  model_dual$lb = c( rep(0, (2*nrow(genIntv_n)) + 1), -Inf)
  model_dual$ub = rep(Inf, (2*nrow(genIntv_n)) + 2)
  
  result_dual = gurobi(model_dual, list(OutputFlag = 0))
  
  return(result_dual)
}  

############### primal lp

primal_lp =  function(x_start, K, data, genIntv_n)
{
  model_primal = list()
  
  model_primal$A = t( A_dual( x = x_start, data = data, genIntv_n = genIntv_n) )
  
  model_primal$obj  = x_start
  model_primal$modelsense = "max"
  model_primal$rhs = c(genIntv_n$conf_up, -genIntv_n$conf_low, K, 1)
  model_primal$sense =  c(rep("<=",(2*nrow(genIntv_n)) + 1), "=")
  
  model_primal$lb = rep(0, length(x_start))
  model_primal$ub = rep(1, length(x_start))
  
  result_primal = gurobi(model_primal, list(OutputFlag = 0))
  
  return(result_primal)
}

############### arg_min for of functions with [low_bd, up_bd]

arg_min = function(lambda, data, genIntv_n , eps = 10^(-7), analyze = F, 
                       low_bd, up_bd)
{
  data_sort = sort(data)
  
  x_temp = NULL
  
  ##### end points of the interval
  
  for(i in 1:(length(data)-1) )
  {
    x_temp = c(x_temp,data_sort[i]+ eps, data_sort[i+1]) 
  }
  
  ##### end points of our serach query
  x_temp = c(x_temp, low_bd, up_bd)
  
  ##### if lambda_H!= 0, we have another candidate to consider
  
  lambda_H = lambda[2*nrow(genIntv_n) + 1]
  
  if(lambda_H != 0){
    x_temp = c(x_temp,H_der_inv(1/lambda_H))
  }
  
  A_dual_temp = A_dual( x = x_temp, 
                        data = data, genIntv_n = genIntv_n)
  slack_x_temp = as.numeric(A_dual_temp%*%lambda - x_temp)
  
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
    return(data.frame(x_min = x_temp[ which.min(slack_x_temp)], min = min(slack_x_temp), 
                      idx = which.min(slack_x_temp), lambda_H = lambda_H))
  }
}

# ############ cutting plane ##########
# 
# n1 = 100
# K1 = 1
# alpha = 0.1
# 
# x1 = rnorm(n1, mean = 0)
# 
# plot(type = "n", x = 0, xlim = c(0,30), ylim = c(0,0.1))
# 
# x1_sort = sort(x1)
# 
# ######### Values of F_n(I) on intervals given by Rivera Walther(2010)
# 
# genIntv_n1 = genIntv(n)
# 
# genIntv_n1$F_n_int = (genIntv_n1$right - genIntv_n1$left )/n
# 
# get_bds_n1 = get_bds_essHist(genIntv_n1$F_n_int , n = n, alpha = alpha)
# genIntv_n1$conf_low = get_bds_n1$low 
# genIntv_n1$conf_up = get_bds_n1$up 
# 
# 
# #### starting optimization
# 
# x1_start = x1_sort
# 
# iter = 0
# multi_factor = 0
# 
# while( TRUE)
# {
#   iter = iter + 1
#   dual_temp = dual_lp(x_start = x1_start, K = K1, data = x1, genIntv_n = genIntv_n1)
#   primal_temp = primal_lp(x_start = x1_start, K = K1, data = x1, genIntv_n = genIntv_n1)
#   
#   lambda_temp = dual_temp$x
#   lambda_H_temp = lambda_temp[2*nrow(genIntv_n1) + 1]
#   
#   print(paste("objval:", round(primal_temp$objval,4)))
#   print(paste("lambda_H:", lambda_H_temp))
#   
#   multi_factor = multi_factor + 1
#   
#   p_temp = primal_temp$x
#   (x1_new = arg_min(lambda = lambda_temp, 
#                         data = x1, 
#                         genIntv_n = genIntv_n1,
#                         low_bd = -multi_factor*H_inv(K1),
#                         up_bd = multi_factor*H_inv(K1),
#                         analyze = F))
#   
#   print(paste("new point", round(x1_new,4)))
#   
#   if(is.na(x1_new)){break}
#   
#   
#   x1_start = c(x1_start, x1_new)
#   
#   print( paste("iter =", iter) )
#   points(x = iter, y = lambda_H_temp)
# }
# 
# primal_temp$objval
# 
# ########################################################
# ########################################################
