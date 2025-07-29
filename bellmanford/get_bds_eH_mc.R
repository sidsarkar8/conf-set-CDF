library(essHist)
library(dplyr)




## l(x) function
l_eh = function(x)
{
  if( x > 1 || x<0){ print(" invalid x")
  }else{
    return( sqrt(2-2*log(x)-2*log(1-x)) ) 
  }
}


KL = function(p,q){
  #val = suppressMessages(distance(rbind( c(p,1-p), c(q,1-q) ), method = "kullback-leibler"))
  val = p*log(p/q) + (1-p)*log((1-p)/(1-q))
  if( p == 0 ){val = log(1/(1-q))}
  if( p == 1 ){val = log(1/q) }
  if( q == 0 ){val = Inf }
  if( q == 1 ){val = Inf}
  
  if( p == 0 & q == 0 ){val = 0}
  if( p == 1 & q == 0 ){val = Inf }
  if( p == 0 & q == 1 ){val = Inf }
  if( p == 1 & q == 1 ){val = 0}
  
  
  
  return( as.numeric(val) )
}

## final function whose roots need to be solved.
fin_func = function(p,q,quant,n)
{
  (2*n*KL(p,q))^(0.5) - l_eh(q) - quant
}

##### quantile of intervals

quant_gen = function( alpha = 0.05, n_runs = 10000, n, F_n_full )
{
  
  sam = rep(-Inf,n_runs)
  
  for( i in 1:n_runs)
  {
    if(!(i%%100)){ print(i) }
    
    U_sort = sort(runif(n))
    t_n = (1:n)/(n)
    
    for(j in 1:nrow(F_n_full))
    {
      p = U_sort[F_n_full$right[j]] - U_sort[F_n_full$left[j]]
      q = t_n[F_n_full$right[j]] - t_n[F_n_full$left[j]]
      sam[i] = max(sam[i], ( (2*n*KL(p,q))^(0.5) - l_eh(q) )  )
    }
  }
  
  n_cuttoff = ceiling(n_runs*(1-alpha))
  
  return(sort(sam)[n_cuttoff])
}


######

get_bds_eH_mc= function(F_n_full, alpha = 0.05, n)
{
  bds = data.frame(low = rep(NA, nrow(F_n_full)),
                   up = rep(NA, nrow(F_n_full)))
  
  quant = quant_gen(F_n_full = F_n_full, alpha = alpha, n = n)
  
  for( i in 1:nrow(F_n_full))
  {
    print(i)
    ### if no p>=0 achieve equality, lower bd = 0
    if(fin_func( p = 0, q = F_n_full$F_n_int[i], quant = quant, n = n  ) < 0 ){
      bds$low[i] = 0
    }else{
      bds$low[i] = uniroot(fin_func , lower = 0 , upper = F_n_full$F_n_int[i], 
                           q = F_n_full$F_n_int[i], quant = quant, n = n, tol = 10^-10)$root
    }
    
    ### if no p<=1 achieve equality, upper bd = 1
    if(fin_func(p = 1, q = F_n_full$F_n_int[i], quant = quant, n = n) < 0 ){
      bds$up[i] = 1
    }else{
      bds$up[i] = uniroot(fin_func , lower = F_n_full$F_n_int[i] , upper = 1, 
                          q = F_n_full$F_n_int[i], quant = quant, n = n, tol = 10^-10)$root
    }
  }
  return( bds )
}


### example
# get_bds(1:9/10,n = 10)
# Quantile simulation might take a while ... ... end!
#   low        up
# 1  0.00000000 0.5158972
# 2  0.00000000 0.6479451
# 3  0.00000000 0.7538756
# 4  0.02029976 0.8439009
# 5  0.07995866 0.9200413
# 6  0.15609906 0.9797002
# 7  0.24612444 1.0000000
# 8  0.35205489 1.0000000
# 9  0.48410284 1.0000000

