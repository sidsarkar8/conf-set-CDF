library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/conf_tail/discrete/fast_sims/hydras/")

x2_20 = read.csv("x2/x2_20.csv")[,-1]
x2_50 = read.csv("x2/x2_50.csv")[,-1]
x2_100 = read.csv("x2/x2_100.csv")[,-1]
x2_250 = read.csv("x2/x2_250.csv")[,-1]
x2_500 = read.csv("x2/x2_500.csv")[,-1]

nrow(x2_20)
nrow(x2_50)
nrow(x2_100)
nrow(x2_250)
nrow(x2_500)

########## 

x2_20_fin = x2_20 %>% mutate(width_dw = up_dw - low_dw,
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

x2_50_fin = x2_50  %>% mutate(width_dw = up_dw - low_dw,
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

x2_100_fin = x2_100  %>% mutate(width_dw = up_dw - low_dw,
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

x2_250_fin = x2_250 %>% mutate(width_dw = up_dw - low_dw,
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

x2_100_fin = x2_100  %>% mutate(width_dw = up_dw - low_dw,
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

x2_250_fin = x2_250  %>% mutate(width_dw = up_dw - low_dw,
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


x2_500_fin = x2_500  %>% mutate(width_dw = up_dw - low_dw,
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

fin_x2 = cbind(
  apply( x2_20_fin, 2, mean),
  apply( x2_50_fin, 2, mean),
  apply( x2_100_fin, 2, mean),
  apply( x2_250_fin, 2, mean),
  apply( x2_500_fin, 2, mean)
)

colnames(fin_x2) = n2_seq


# plot(log(n2_seq), fin_x2[1,], type = "l", ylim = c(min(fin_x2[c(1,2,4),]), max(fin_x2[c(1,2,4),])),
#      xlab = "log(n)", ylab = "Width", main = "Width comparison for bounded second moment")
# points(log(n2_seq), fin_x2[2,], type = "l", col = "red")
# points(log(n2_seq), fin_x2[4,], type = "l", col = "blue") 
# legend(x = 5, y =  4, fill = c("black", "red", "blue"),
#        legend = c("CDF method","CDF method (data)", "t-test"))
  
