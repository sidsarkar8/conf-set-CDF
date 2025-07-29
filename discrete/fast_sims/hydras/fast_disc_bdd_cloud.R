##### method
# sparse matrices
# pass warmstart info
# have a good start, and check only H_inv(1/lambda_H)

###### Dumbgen Wellner

library(doParallel)
library(gurobi)
library(Matrix)
library(MASS)
library(dplyr)

########################################################################
########################################################################

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

####################################################################################
################# int lb and ub ###################

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


up_bd_dw_bdd = function(x1, alpha = 0.05, analyze = F,
                        genCDF_n1, eps = 10^-6, low_bd, up_bd)
{
  n = length(x1)
  x1_sort = sort(x1)
  
  #### starting optimization
  ### sorted data point only sufficient to check
  x1_start = c(low_bd, x1_sort - eps, x1_sort, x1_sort + eps, 
               up_bd)
  
  dual_temp = dual_lp_dw_bdd(x_start = x1_start, data = x1, genCDF_n = genCDF_n1)
  
  return( dual_temp$objval)
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

parallel_fin_func = function(simul_idx, n2, genCDF_n2)
{
  set.seed(simul_idx)
  
  x2 = numeric(n2)
  mix1 = rbeta(n2, shape1 = 1000, shape2 = 2)
  mix2 = rbeta(n2, shape2 = 1000, shape1 = 2)
  toss_mix = runif(n2)
  
  for( i in 1:n2)
  {
    x2[i] = ifelse(toss_mix[i]< 0.5, mix1[i], mix2[i])
  }
  
  lb1 = 0
  ub1 = 1
  
  
  {
    t_dw_bdd = Sys.time()
    up_dw_bdd = up_bd_dw_bdd(x1 = x2, low_bd = lb1, up_bd = ub1,
                             alpha = alpha, genCDF_n1 = genCDF_n2)
    low_dw_bdd = -up_bd_dw_bdd(x1 = -x2, low_bd = -ub1, up_bd = -lb1,
                               alpha = alpha, genCDF_n1 = genCDF_n2)
    t_dw_bdd = Sys.time() - t_dw_bdd
  }
  
  {
    t_dw_x1_to_xn = Sys.time()
    up_dw_x1_to_xn = up_bd_dw_x1_to_xn(x1 = x2,alpha = alpha, 
                                       genCDF_n1 = genCDF_n2)
    low_dw_x1_to_xn = -up_bd_dw_x1_to_xn(x1 = -x2, alpha = alpha, 
                                         genCDF_n1 = genCDF_n2)
    
    
    t_dw_x1_to_xn = Sys.time() - t_dw_x1_to_xn
  }
  
  up_tail = (sort(x2)[1])*genCDF_n2$conf_up[1] + ub1*(1 - genCDF_n2$conf_low[n2])
    
  low_tail = (lb1)*genCDF_n2$conf_low[1] + (sort(x2)[n2])*(1 - genCDF_n2$conf_up[n2])

  
  { 
    t_t_test = Sys.time()
    CI_t = t.test(x2, conf.level = 1 - alpha)
    CI_t = CI_t$conf.int[1:2]
    t_t_test = Sys.time() - t_t_test
  }
  
  
  return(data.frame("simul_idx" = simul_idx, 
                    "up_dw" = up_dw_bdd, 
                    "low_dw" = low_dw_bdd,
                    "t_dw" = t_dw_bdd,
                    "up_dw_x1_to_xn" = up_dw_x1_to_xn,
                    "low_dw_x1_to_xn" = low_dw_x1_to_xn,
                    "t_dw_x1_to_xn" = t_dw_x1_to_xn,
                    "up_dw_x1_to_xn_tails" = up_dw_x1_to_xn + 
                      up_tail,
                    "low_dw_x1_to_xn_tails" = low_dw_x1_to_xn + 
                      low_tail,
                    "t_dw_x1_to_xn_tails" = t_dw_x1_to_xn,
                    "up_t_test" = CI_t[2], 
                    "low_t_test" = CI_t[1],
                    "t_t_test" = t_t_test))
}


n2_seq = c(20,50,100,250,500)

M = 500

numCores = detectCores() - 1
registerDoParallel(numCores) 

for( n2_idx in 1:length(n2_seq))
{
  print(paste(n2_idx, ": run sims"))
  
  alpha = 0.1
  
  get_bds_CDF_n2 = get_bds_dw(n = n2_seq[n2_idx],
                              alpha = alpha)
  genCDF_n2 = data.frame(conf_low = get_bds_CDF_n2$low,
                         conf_up = get_bds_CDF_n2$up)
  
  
  df_n2_bdd = foreach( simul_idx = 1:M,  .combine = rbind ) %dopar%{
    parallel_fin_func(simul_idx, 
                      n2 = n2_seq[n2_idx],
                      genCDF_n2 = genCDF_n2)
  }
  
  write.csv(df_n2_bdd,
            file = paste("/home/siddhaarth/conf_tail/fast_sims/bdd_",
                         n2_seq[n2_idx],".csv",sep = "" ))
}

stopImplicitCluster()

