##### discretization for general
### comparison for method
### xlogx

######### xlogx

H = function(x){return(abs(x)*log(abs(x))) }
H_inv = function(x){return( exp(lambertWp(x)) )}
H_der = function(x){return( (1 + log(abs(x)))*sign(x) )}
H_der_inv = function(x){return( exp(x-1) )}

######## CDF methods #########

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

#######################
# 
# x11 = rdistr(1000)
# 
# mean(x11^2)
# 
# plot(ecdf(x11))
# points(sort(x11), F_distr(sort(x11)), col = "red", type = "l")
# 
source("~/Documents/Projects/conf_tail/discrete/disc_algo_general.R")

############ cutting plane ##########

optimize_up_bd = function(x1, K1, alpha = 0.05, analyze = F)
{  
  
  n = length(x1)
  x1_sort = sort(x1)
  
  ######### Values of F_n(I) on intervals given by Rivera Walther(2010)
  
  genIntv_n1 = genIntv(n)
  
  genIntv_n1$F_n_int = (genIntv_n1$right - genIntv_n1$left )/n
  
  get_bds_n1 = get_bds_essHist(genIntv_n1$F_n_int , n = n, alpha = alpha)
  genIntv_n1$conf_low = get_bds_n1$low
  genIntv_n1$conf_up = get_bds_n1$up
  
  #### starting optimization
  
  x1_start = x1_sort
  
  iter = 0
  multi_factor = 0
  
  while( TRUE)
  {
    iter = iter + 1
    dual_temp = dual_lp(x_start = x1_start, K = K1, data = x1, 
                        genIntv_n = genIntv_n1)
    #primal_temp = primal_lp(x_start = x1_start, K = K1, data = x1, genIntv_n = genIntv_n1)
    
    lambda_temp = dual_temp$x
    lambda_H_temp = lambda_temp[2*nrow(genIntv_n1) + 1]
    
    print(paste("objval:", round(dual_temp$objval,4)))
    print(paste("lambda_H:", lambda_H_temp))
    
    multi_factor = multi_factor + 1
    
    (x1_new = arg_min(lambda = lambda_temp,
                      data = x1,
                      genIntv_n = genIntv_n1,
                      low_bd = -multi_factor*H_inv(K1),
                      up_bd = multi_factor*H_inv(K1),
                      analyze = F))
    
    print(paste("new point", round(x1_new,4)))
    
    if(is.na(x1_new)){break}
    
    
    x1_start = c(x1_start, x1_new)
    
    print( paste("iter =", iter) )
    #points(x = iter, y = lambda_H_temp)
  }
  
  if(!analyze)
  { 
    return(dual_temp$objval)
  }else{
    primal_temp = primal_lp(x_start = x1_start, K = K1, data = x1, genIntv_n = genIntv_n1)
    return(list(
      lambda = lambda_temp, x1_opt = x1_start, 
      p = primal_temp$x, objval = dual_temp$objval
    ))
  }
}

########################################################
########################################################


##############

N_sam = 100

width_all_xlogx = NULL

n2 = 50
alpha = 0.1
K2 = 2*q1 + 1 


for(sam_idx in 1:N_sam)
{
  x2 = rdistr(n2)
  x2 = H_inv(x2^q1)*(2*rbinom(n2, size = 1, prob = 1/2) -1)
  #hist(x2)
  
  up_mu = optimize_up_bd(x1 = x2, K1 = K2, alpha = alpha)
  low_mu = -optimize_up_bd(x1 = -x2, K1 = K2, alpha = alpha)
  
  width_all_xlogx = c(width_all_xlogx, up_mu - low_mu)
}

write.csv(width_all_xlogx, "~/Documents/Projects/conf_tail/discrete/sims/width_comb_bdd_xlog.csv")

