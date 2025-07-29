#### comparing bellman ford to exhaustive search

setwd("/Users/sidsarkar/Documents/Projects/conf_tail/bellmanford")
source("get_bds_essHist.R")
source("get_bds_eH_mc.R")
source("bellmanford_genfunc.R")
library(devtools)
library(confseq)
library(CVXR)
#set.seed(1)

n = 65
#x = runif(n, min = -2, max = 2)
x = rnorm(n)
x_sort = sort(x)

func = function(x){return(x)}

## Values of F_n(I) on intervals given by Rivera Walther(2010)
genIntv_n = genIntv(n)

genIntv_n$F_n_int = (genIntv_n$right - genIntv_n$left )/n

get_bds_n = get_bds_essHist(genIntv_n$F_n_int , n = n, alpha = 0.05)
#get_bds_n = get_bds_eH_mc(genIntv_n , n = n, alpha = 0.05)
genIntv_n$conf_low = get_bds_n$low 
genIntv_n$conf_up = get_bds_n$up 


genIntv_n$func_max = NA
genIntv_n$func_min = NA

### storing max of functions of interval
for( row_idx in 1:nrow(genIntv_n))
{
  i = genIntv_n$left[row_idx]
  j = genIntv_n$right[row_idx]
  
  x_ij_grid = seq(from = x_sort[i], to = x_sort[j], by = 1/(n^2))
  
  genIntv_n$func_max[row_idx] = max(func(x_ij_grid))
  genIntv_n$func_min[row_idx] = min(func(x_ij_grid))
}
###########

bf_up_n = bf_up(genIntv = genIntv_n, n = n)
bf_low_n = bf_low(genIntv = genIntv_n, n = n)

pred_up = bf_up_n$pred
pred_up[is.na(pred_up)] = 1

pred_low = bf_low_n$pred
pred_low[is.na(pred_low)] = 1

### tracking path
bf_up_path = NULL
row_idx_bf_up_path = NULL

bf_low_path = NULL
row_idx_bf_low_path = NULL

temp_pred_up = n
while(temp_pred_up != 1)
{
  bf_up_path = rbind(bf_up_path,
                     c(pred_up[temp_pred_up], temp_pred_up))
  row_idx_bf_up_path = c(row_idx_bf_up_path,
                         which(genIntv_n$left == pred_up[temp_pred_up] & genIntv_n$right == temp_pred_up ))
  ## move to next path
  temp_pred_up = pred_up[temp_pred_up]
}

###############


temp_pred_low = n
while(temp_pred_low != 1)
{
  bf_low_path = rbind(bf_low_path,
                      c(pred_low[temp_pred_low], temp_pred_low))
  row_idx_bf_low_path = c(row_idx_bf_low_path,
                          which(genIntv_n$left == pred_low[temp_pred_low] & genIntv_n$right == temp_pred_low ))
  ## move to next path
  temp_pred_low = pred_low[temp_pred_low]
}


### to plot the best path
plot(x_sort, func(x_sort), main = "norm")

for( row_idx in row_idx_bf_up_path)
{
  i = genIntv_n$left[row_idx]
  j = genIntv_n$right[row_idx]
  
  x_ij_grid = seq(from = x_sort[i], to = x_sort[j], by = 1/n)
  
  points(x = x_ij_grid, y = rep(genIntv_n$func_max[row_idx],length(x_ij_grid)), 
         col = "red", type = "l")
}


for( row_idx in row_idx_bf_low_path)
{
  i = genIntv_n$left[row_idx]
  j = genIntv_n$right[row_idx]
  
  x_ij_grid = seq(from = x_sort[i], to = x_sort[j], by = 1/n)
  
  points(x = x_ij_grid, y = rep(genIntv_n$func_min[row_idx],length(x_ij_grid)), 
         col = "blue", type = "l")
}

#bf_up_path

#########################
#########################

library(igraph)
library(dplyr)

g_intv = graph_from_edgelist(genIntv_n[,c("left","right")] %>% as.matrix, directed = T)
#plot(g_intv)
paths = all_simple_paths(g_intv, from = 1, to = n)

####### max
best_up_sum = Inf
best_up_path_idx = 0
for( i in 1:length(paths))
{ temp_path = paths[[i]] %>% as.numeric()
temp_path_sum = 0
for(j in 1:(length(temp_path)-1))
{
  row_idx= which(genIntv_n$left == temp_path[j] & genIntv_n$right == temp_path[j+1] )
  
  dist_row_idx = ifelse(genIntv_n$func_max[row_idx]> 0,
                       (genIntv_n$func_max[row_idx]*genIntv_n$conf_up[row_idx]),
                       (genIntv_n$func_max[row_idx]*genIntv_n$conf_low[row_idx]))
  temp_path_sum = temp_path_sum + dist_row_idx
}

if(temp_path_sum < best_up_sum)
{
  best_up_sum = temp_path_sum
  best_up_path_idx = i
}
}
############ min

best_low_sum = -Inf
best_low_path_idx = 0
for( i in 1:length(paths))
{ temp_path = paths[[i]] %>% as.numeric()
temp_path_sum = 0
for(j in 1:(length(temp_path)-1))
{
  row_idx= which(genIntv_n$left == temp_path[j] & genIntv_n$right == temp_path[j+1] )
  dist_row_idx = ifelse(genIntv_n$func_min[row_idx]> 0,
                       (genIntv_n$func_min[row_idx]*genIntv_n$conf_low[row_idx]),
                       (genIntv_n$func_min[row_idx]*genIntv_n$conf_up[row_idx]))
  
  
  temp_path_sum = temp_path_sum + dist_row_idx
}

if(temp_path_sum > best_low_sum)
{
  best_low_sum = temp_path_sum
  best_low_path_idx = i
}
}

#############

bf_up_path
paths[[best_up_path_idx]]

bf_low_path

paths[[best_low_path_idx]]

c(bf_low_n$dist[n], bf_up_n$dist[n])
mean(x)

( bernoulli_confidence_interval(sum((x)), n, alpha = 0.05, t_opt = n))
