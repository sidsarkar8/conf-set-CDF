##### discretization for bounded RV 
### comparison for method

######## CDF methods #########

source("~/Documents/Projects/conf_tail/discrete/disc_algo_bdd.R")


bdd_optimize_up_bd = function(x1, low_bd, up_bd, alpha = 0.05)
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
  
  while( TRUE)
  {
    iter = iter + 1
    dual_temp = dual_lp_bdd(x_start = x1_start, data = x1, genIntv_n = genIntv_n1)
    primal_temp = primal_lp_bdd(x_start = x1_start, data = x1, genIntv_n = genIntv_n1)
    
    lambda_temp = dual_temp$x
    
    print(paste("objval:", round(primal_temp$objval,4)))
    
    p_temp = primal_temp$x
    (x1_new = arg_min_bdd(lambda = lambda_temp,
                          data = x1,
                          genIntv_n = genIntv_n1,
                          low_bd = low_bd,
                          up_bd = up_bd,
                          analyze = F))
    
    print(paste("new point", round(x1_new,4)))
    
    if(is.na(x1_new)){break}
    
    
    x1_start = c(x1_start, x1_new)
    
    print( paste("iter =", iter) )
    
  }
  
  return( primal_temp$objval)
}

######### Hulc #########

B_alpha_func = function(alpha = 0.05)
{
  B.alpha = ceiling(log(2/alpha, base = 2))
  
  ### rand_p is p such that alpha/2 = p*( 2^B ) + (1-p)*( 2^{B-1} )
  rand_p = 2 - (alpha/2)*2^(B.alpha)
  if( runif(1) > rand_p) {
    B.alpha = B.alpha - 1
  }
  return(B.alpha)
}

CI_hulc = function(x, alpha = 0.05)
{
  B_alpha = B_alpha_func(alpha = alpha) 
  
  f_hat_intervals = rep(0, B_alpha)
  split_points = split( sample(length(x)), (1:length(x))%%B_alpha )
  
  
  for(b_idx in 1:B_alpha){
    f_hat_intervals[b_idx] = mean( x[split_points[[b_idx]]] )
  }
  
  return(range(f_hat_intervals))
}

######### 

################# Anderson's method

source("~/Documents/Projects/Conf/dumbgun_wellner_new/get_bds_dw1.R")
source("~/Documents/Projects/Conf/dkw/get_bds_dkw.R")

CI_anderson = function(x, alpha = 0.05, low_bd, up_bd)
{
  n = length(x)
  x_sc = (x - low_bd)/(up_bd-low_bd) 
  
  bds_dw = get_bds_dkw_simp( n = n, alpha = alpha)  
  
  x_sc_ord = c(0, sort(x_sc), 1 )
  x_sc_ord_diff = (x_sc_ord - lag(x_sc_ord,1))[2:(n+2)]
  
  up_bd_wts = c(1, 1 - bds_dw$low)
  CI_up_bd = sum(up_bd_wts*x_sc_ord_diff)
  
  low_bd_wts = c(1-bds_dw$up,0)
  CI_low_bd = sum(low_bd_wts*x_sc_ord_diff)
  
  return(  c(CI_low_bd,CI_up_bd)*(up_bd-low_bd) + low_bd  )
}

########################################################
########################################################

CI_cheby = function(x, alpha = 0.05, K1)
{
  n = length(x)
  
  cntr = (mean(x)*alpha)/(alpha + (1/n))
  dev_cntr = (alpha*(K1 - mean(x)^2))/n + K1/n^2
  dev_cntr = sqrt( dev_cntr )/(alpha + (1/n))
  
  return(cntr + c(-1,1)*dev_cntr )
}

######################################################################################
######################################################################################

n2 = 100
alpha = 0.1

low_bd1 = 0
up_bd1 = 1

#x2 = runif(n2, min = -1, max = 1)
#x2 = rbeta(n2, shape1 = 1000, shape2 = 2)
x2 = numeric(n2)
mix1 = rbeta(n2, shape1 = 1000, shape2 = 2)
mix2 = rbeta(n2, shape2 = 1000, shape1 = 2)
toss_mix = runif(n2)
for( i in 1:n2)
{
  x2[i] = ifelse(toss_mix[i]< 0.5, mix1[i], mix2[i])
}

hist(x2)


##########################################################

up_mu = bdd_optimize_up_bd(x1 = x2, 
                           low_bd = low_bd1, up_bd = up_bd1, alpha = alpha)
low_mu = -bdd_optimize_up_bd(x1 = -x2, 
                             low_bd = -up_bd1, up_bd = -low_bd1, alpha = alpha)

CI_optim = c(low_mu, up_mu)

CI_t = t.test(x2, conf.level = 1 - alpha)
CI_t = CI_t$conf.int[1:2]

CI_cheby1 = CI_cheby(x = x2, alpha = alpha, K1 = max(abs(up_bd1), abs(low_bd1)))

#CI_hulc_1 = CI_hulc(x2, alpha)

CI_and = CI_anderson(x2 , alpha = alpha, low_bd = low_bd1, up_bd = 1)


df_compare = rbind(CI_optim,CI_cheby1,CI_and,CI_t) %>% 
  as.data.frame() %>% mutate(V1 = ifelse(V1 < low_bd1, low_bd1, V1), 
                             V2 = ifelse(V2 > up_bd1, up_bd1, V2)) %>% 
  mutate(width = V2 - V1) %>% rename(low = V1, up = V2)

df_compare
