setwd("/Users/sidsarkar/Documents/Projects/conf_tail/bellmanford")
source("get_bds_essHist.R")
source("get_bds_eH_mc.R")
source("bellmanford_genfunc.R")
library(devtools)
library(confseq)
set.seed(1)

n = 100
#x = runif(n, min = -2, max = 2)
x = rnorm(n)
#func = function(x){return( rep(1,length(x)) ) }
#func = function(x){ return((sin(x) + 1)/2)}
func = function(x){return(x)}
#func = function(x){return(exp(x)/exp(1))}
#func = exp
# func = function(x)
# {
#   intervals = c(0.1,0.3,0.4,0.8)
#   temp = numeric(length(x))
#   
#   for( i in 1:length(x))
#   {
#     #print(i)
#     if(x[i]<intervals[1]){
#       temp[i] = 1.2
#     }else if( x[i]< intervals[2]){
#       temp[i] =2.42
#     }else if(x[i]< intervals[3]){
#       temp[i] =5
#     }else if(x[i]< intervals[4]){
#       temp[i] = (0.2)
#     }else{ 
#       temp[i] = 2}
#   }
#   return(temp/7)
# }


#x = exp(x)
#func = function(x){return(x)}

#x = pnorm(x)
#func = function(x){return(exp(qnorm(x)))}

## Values of F_n(I) on intervals given by Rivera Walther(2010)
genIntv_n = genIntv(n)

## Values of F_n(I) on intervals given by all intervals
#genIntv_n = genIntv(n, type= "Full")

### values to get normal conf band
#genIntv_n = data.frame(left = rep(1, n-1), right = 2:n)

#genIntv_n = data.frame(left = bf_up_path[,1], right = bf_up_path[,2])

genIntv_n$F_n_int = (genIntv_n$right - genIntv_n$left )/n

get_bds_n = get_bds_essHist(genIntv_n$F_n_int , n = n, alpha = 0.05)
#get_bds_n = get_bds_eH_mc(genIntv_n , n = n, alpha = 0.05)
genIntv_n$conf_low = get_bds_n$low 
genIntv_n$conf_up = get_bds_n$up 

x_sort = sort(x)

genIntv_n$func_max = NA
genIntv_n$func_min = NA

### storing max of functions of interval
for( row_idx in 1:nrow(genIntv_n))
{
  i = genIntv_n$left[row_idx]
  j = genIntv_n$right[row_idx]
  
  x_ij_grid = seq(from = x_sort[i], to = x_sort[j], by = 1/n)
  
  genIntv_n$func_max[row_idx] = max(func(x_ij_grid))
  genIntv_n$func_min[row_idx] = min(func(x_ij_grid))
  
  genIntv_n$func_avg[row_idx] = sum(func(x_sort[(i+1):j]))/(j-i)
  
}

#### plotting to check stuff

plot(x_sort, func(x_sort))

for( row_idx in 1:nrow(genIntv_n))
{
  i = genIntv_n$left[row_idx]
  j = genIntv_n$right[row_idx]
  
  x_ij_grid = seq(from = x_sort[i], to = x_sort[j], by = 1/n)
  
  points(x = x_ij_grid, y = rep(genIntv_n$func_max[row_idx],length(x_ij_grid)), 
         col = "red", type = "l")
  
  points(x = x_ij_grid, y = rep(genIntv_n$func_min[row_idx],length(x_ij_grid)), 
         col = "blue", type = "l")
}

###########

bf_up_n = bf_up(genIntv = genIntv_n, n = n)
pred = bf_up_n$pred
pred[is.na(pred)] = 1

### tracking path
bf_up_path = NULL
row_idx_bf_up_path = NULL

temp_pred = n
while(temp_pred != 1)
{
  bf_up_path = rbind(bf_up_path,
                     c(pred[temp_pred], temp_pred))
  row_idx_bf_up_path = c(row_idx_bf_up_path,
                         which(genIntv_n$left == pred[temp_pred] & genIntv_n$right == temp_pred ))
  ## move to next path
  temp_pred = pred[temp_pred]
}


#######
### if we look at step 1->2->3 ... -> n, what value are we getting
row_idx_step = NULL

idx_path = c(1,25,50,75,100)
for( i in 1:length(idx_path))
{
  row_idx_step = c(row_idx_step,
                   which(genIntv_n$left == idx_path[i] & genIntv_n$right == idx_path[(i+1)] ))
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

# temp_sum = 0
# for( row_idx in row_idx_step)
# {
#   i = genIntv_n$left[row_idx]
#   j = genIntv_n$right[row_idx]
#   
#   x_ij_grid = seq(from = x_sort[i], to = x_sort[j], by = 1/n)
#   
#   points(x = x_ij_grid, y = rep(genIntv_n$func_max[row_idx],length(x_ij_grid)),
#          col = "blue", type = "l", lty = 2)
#   temp_sum = temp_sum + (genIntv_n$func_max[row_idx]*genIntv_n$conf_up[row_idx])
# }
# temp_sum

bf_up_path

confint(lm(y~1, data.frame(y = func(x_sort))))

bf_up_n$dist[n]

#( bernoulli_confidence_interval(sum(func(x_sort)), n, alpha = 0.05, t_opt = n))

