##### discretization for general rv 
### comparison for method
### x^2
### combining confidence sets

######### x^2

H = function(x){return(x^2)}
H_inv = function(x){return(sqrt(x))}
H_der = function(x){return(2*x)}
H_der_inv = function(x){return(x/2)}

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

x11 = rdistr(100000)

mean(x11^2)

plot(ecdf(x11))
points(sort(x11), F_distr(sort(x11)), col = "red", type = "l")

source("~/Documents/Projects/conf_tail/discrete/disc_algo_comb.R")


######## starting optimization

optimize_comb_up_bd = function(x1, K1, alpha = 0.05, analyze = F, genIntv_n1, genCDF_n1)
{
  
  n = length(x1)
  x1_sort = sort(x1)
  x1_start = x1_sort

  iter = 0
  multi_factor = 0
  
  while( TRUE)
  {
    iter = iter + 1
    dual_temp = dual_lp_comb(x_start = x1_start, K = K1, data = x1,
                             genIntv_n = genIntv_n1,
                             genCDF_n = genCDF_n1)
    
    lambda_temp = dual_temp$x
    lambda_H_temp = lambda_temp[(2*nrow(genIntv_n1)) + (2*nrow(genCDF_n1)) +  1]
    
    print(paste("objval:", round(dual_temp$objval,4)))
    print(paste("lambda_H:", lambda_H_temp))
    
    multi_factor = multi_factor + 1
    
    #p_temp = primal_temp$x
    (x1_new = arg_min_comb(lambda = lambda_temp,
                           data = x1,
                           genIntv_n = genIntv_n1,
                           genCDF_n = genCDF_n1,
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
    primal_temp = primal_lp_comb(x_start = x1_start, K = K1, data = x1,
                                genIntv_n = genIntv_n1,
                                genCDF_n = genCDF_n1)
    return(list(
      lambda = lambda_temp, x1_opt = x1_start, 
      p = primal_temp$x, objval = dual_temp$objval
    ))
  }
}

############################ 

n2 = 100
alpha = 0.1

x2 = rdistr(n2)
K2 = 1 + 2*q1 


############ combined genIntv_n 

genIntv_comb = genIntv(n2)
genIntv_comb$F_n_int = (genIntv_comb$right - genIntv_comb$left )/n2
get_bds_eH_comb = get_bds_essHist(genIntv_comb$F_n_int, n = n2, alpha = alpha/2)
genIntv_comb$conf_low = get_bds_eH_comb$low
genIntv_comb$conf_up = get_bds_eH_comb$up

get_bds_CDF_comb = get_bds_dkw_simp( n = n2, alpha = alpha/2)
genCDF_comb = data.frame(conf_low = get_bds_CDF_comb$low,
                       conf_up = get_bds_CDF_comb$up)

####### ess Hist

genIntv_eH = genIntv(n2)
genIntv_eH$F_n_int = (genIntv_eH$right - genIntv_eH$left )/n2
get_bds_eH_eH = get_bds_essHist(genIntv_eH$F_n_int , n = n2, alpha = alpha)
genIntv_eH$conf_low = get_bds_eH_eH$low
genIntv_eH$conf_up = get_bds_eH_eH$up

genCDF_eH = data.frame(conf_low = rep(0,n2),
                       conf_up = rep(1,n2))

#### Anderson

genIntv_and = genIntv(n2)
genIntv_and$F_n_int = (genIntv_and$right - genIntv_and$left )/n2
genIntv_and$conf_low = rep(0, nrow(genIntv_and))
genIntv_and$conf_up = rep(1, nrow(genIntv_and))

get_bds_CDF_and = get_bds_dkw_simp( n = n2, alpha = alpha)
genCDF_and = data.frame(conf_low = get_bds_CDF_and$low,
                      conf_up = get_bds_CDF_and$up)


#######

up_comb = optimize_comb_up_bd(x2, K2, alpha = 0.05, analyze = F, genIntv_comb, genCDF_comb)
low_comb = -optimize_comb_up_bd(-x2, K2, alpha = 0.05, analyze = F, genIntv_comb, genCDF_comb)

up_eH = optimize_comb_up_bd(x2, K2, alpha = 0.05, analyze = F, genIntv_eH, genCDF_eH)
low_eH = -optimize_comb_up_bd(-x2, K2, alpha = 0.05, analyze = F, genIntv_eH, genCDF_eH)

up_and = optimize_comb_up_bd(x2, K2, alpha = 0.05, analyze = F, genIntv_and, genCDF_and)
low_and = -optimize_comb_up_bd(-x2, K2, alpha = 0.05, analyze = F, genIntv_and, genCDF_and)

#######

CI_comb = c(low_comb, up_comb)
CI_eH = c(low_eH, up_eH)
CI_and = c(low_and, up_and)

CI_comb
CI_eH
CI_and