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
  x1_start = c(low_bd, x1_sort - eps, x1_sort, x1_sort + eps, up_bd)
  
  dual_temp = dual_lp_dw_bdd(x_start = x1_start, data = x1, genCDF_n = genCDF_n1)
  
  return( dual_temp$objval)
}


#########################################################################
#########################################################################

n2 = 100
alpha = 0.1

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

get_bds_CDF_n2 = get_bds_dw(n = n2, alpha = alpha)
genCDF_n2 = data.frame(conf_low = get_bds_CDF_n2$low,
                       conf_up = get_bds_CDF_n2$up)


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

CI_t = t.test(x2, conf.level = 1 - alpha)
CI_t = CI_t$conf.int[1:2]

c( low_dw_bdd, up_dw_bdd)
c( low_dw_x1_to_xn, up_dw_x1_to_xn)
CI_t