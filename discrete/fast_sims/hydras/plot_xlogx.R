library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/conf_tail/discrete/fast_sims/hydras/")

xlogx_20 = read.csv("xlogx/xlogx_20.csv")[,-1]
xlogx_50 = read.csv("xlogx/xlogx_50.csv")[,-1]
xlogx_100 = read.csv("xlogx/xlogx_100.csv")[,-1]
xlogx_250 = read.csv("xlogx/xlogx_250.csv")[,-1]
xlogx_500 = read.csv("xlogx/xlogx_500.csv")[,-1]

nrow(xlogx_20)
nrow(xlogx_50)
nrow(xlogx_100)
nrow(xlogx_250)
nrow(xlogx_500)

########## 

xlogx_20_fin = xlogx_20 %>% mutate(width_dw = up_dw - low_dw,
                             width_dw_x1_to_xn = up_dw_x1_to_xn - 
                               low_dw_x1_to_xn,
                             width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                               low_dw_x1_to_xn_tails,
                             width_t_test = up_t_test - low_t_test,
                             cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                             cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                             cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                             cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 

xlogx_50_fin = xlogx_50  %>% mutate(width_dw = up_dw - low_dw,
                              width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                low_dw_x1_to_xn,
                              width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                low_dw_x1_to_xn_tails,
                              width_t_test = up_t_test - low_t_test,
                              cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                              cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                              cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                              cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 

xlogx_100_fin = xlogx_100  %>% mutate(width_dw = up_dw - low_dw,
                                width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                  low_dw_x1_to_xn,
                                width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                  low_dw_x1_to_xn_tails,
                                width_t_test = up_t_test - low_t_test,
                                cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                                cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                                cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                                cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 

xlogx_250_fin = xlogx_250 %>% mutate(width_dw = up_dw - low_dw,
                               width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                 low_dw_x1_to_xn,
                               width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                 low_dw_x1_to_xn_tails,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                               cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                               cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                               cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")


########## 

xlogx_100_fin = xlogx_100  %>% mutate(width_dw = up_dw - low_dw,
                                width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                  low_dw_x1_to_xn,
                                width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                  low_dw_x1_to_xn_tails,
                                width_t_test = up_t_test - low_t_test,
                                cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                                cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                                cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                                cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 

xlogx_250_fin = xlogx_250  %>% mutate(width_dw = up_dw - low_dw,
                                width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                  low_dw_x1_to_xn,
                                width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                  low_dw_x1_to_xn_tails,
                                width_t_test = up_t_test - low_t_test,
                                cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                                cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                                cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                                cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")


xlogx_500_fin = xlogx_500  %>% mutate(width_dw = up_dw - low_dw,
                                width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                  low_dw_x1_to_xn,
                                width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                  low_dw_x1_to_xn_tails,
                                width_t_test = up_t_test - low_t_test,
                                cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                                cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0 & low_dw_x1_to_xn <0, 1, 0),
                                cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0 & low_dw_x1_to_xn_tails <0, 1, 0),
                                cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 


n2_seq = c(20,50,100,250,500)

fin_xlogx = cbind(
  apply( xlogx_20_fin, 2, mean),
  apply( xlogx_50_fin, 2, mean),
  apply( xlogx_100_fin, 2, mean),
  apply( xlogx_250_fin, 2, mean),
  apply( xlogx_500_fin, 2, mean)
)

colnames(fin_xlogx) = n2_seq



# plot(log(n2_seq), fin_xlogx[1,], type = "l", ylim = c(min(fin_xlogx[c(1,2,4),]), max(fin_xlogx[c(1,2,4),])),
#      xlab = "log(n)", ylab = "Width", main = "Width comparison for bounded E[|X|log|X|]")
# points(log(n2_seq), fin_xlogx[2,], type = "l", col = "red")
# points(log(n2_seq), fin_xlogx[4,], type = "l", col = "blue") 
# legend(x = 5, y =  4, fill = c("black", "red", "blue"),
#        legend = c("CDF method","CDF method (data)", "t-test"))

