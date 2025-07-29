bf_low = function( genIntv_n,n )
{
  dist = rep(-Inf, n)
  pred = rep(NA,n)
  dist[1] = 0


  for( rep in 1:(n-1))
  {
    for( row_idx in 1:nrow(genIntv_n))
    {
      i = genIntv_n$left[row_idx]
      j = genIntv_n$right[row_idx]
      
      dist_ij = ifelse(genIntv_n$func_min[row_idx]> 0,
                       (genIntv_n$func_min[row_idx]*genIntv_n$conf_low[row_idx]),
                       (genIntv_n$func_min[row_idx]*genIntv_n$conf_up[row_idx]))
      

      if(dist[j] < dist[i] + dist_ij ){pred[j] = i}
      dist[j] = max(dist[j], dist[i] + dist_ij )

    }
  }
  ans = data.frame(dist = dist, pred = pred)
  return(ans)
}

##################

bf_up = function( genIntv_n, n)
{
  dist = rep(Inf, n)
  pred = rep(NA,n)
  
  dist[1] = 0
  
  for( rep in 1:(n-1))
  {
    for( row_idx in 1:nrow(genIntv_n))
    {
      i = genIntv_n$left[row_idx]
      j = genIntv_n$right[row_idx]
      
      dist_ij = ifelse(genIntv_n$func_max[row_idx]> 0,
                       (genIntv_n$func_max[row_idx]*genIntv_n$conf_up[row_idx]),
                       (genIntv_n$func_max[row_idx]*genIntv_n$conf_low[row_idx]))
      # if(dist[j] > dist[i] + (genIntv_n$func_max[row_idx]*genIntv_n$conf_up[row_idx]) ){
      #   pred[j] = i}
      # dist[j] = min(dist[j],dist[i] + (genIntv_n$func_max[row_idx]*genIntv_n$conf_up[row_idx]) )
      
      if(dist[j] > dist[i] + dist_ij ){
         pred[j] = i}
       dist[j] = min(dist[j],dist[i] + dist_ij )
      
   }
  }
  ans = data.frame(dist = dist, pred = pred)
  return(ans)
}


# ##################
# 
# bf_low = function( genIntv_n, n)
# {
#   dist = rep(Inf, n)
#   pred = rep(NA,n)
#   
#   dist[1] = 0
#   
#   for( rep in 1:(n-1))
#   {
#     for( row_idx in 1:nrow(genIntv_n))
#     {
#       i = genIntv_n$left[row_idx]
#       j = genIntv_n$right[row_idx]
#       
#       
#       
#       if(dist[j] > dist[i] - (genIntv_n$func_min[row_idx]*genIntv_n$conf_low[row_idx]) ){
#         pred[j] = i}
#       dist[j] = min(dist[j],dist[i] - (genIntv_n$func_min[row_idx]*genIntv_n$conf_low[row_idx]) )
#       
#       
#     }
#   }
#   ans = data.frame(dist = -dist, pred = pred)
#   return(ans)
# }
# 

