library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/conf_tail/discrete/fast_sims/hydras/")

var_20 = read.csv("var/var_20.csv")[,-1]
var_50 = read.csv("var/var_50.csv")[,-1]
var_100 = read.csv("var/var_100.csv")[,-1]
var_250 = read.csv("var/var_250.csv")[,-1]
var_500 = read.csv("var/var_500.csv")[,-1]

nrow(var_20)
nrow(var_50)
nrow(var_100)
nrow(var_250)
nrow(var_500)

########## 

var_20_fin = var_20 %>% mutate(width_dw = up_dw - low_dw,
                             width_t_test = up_t_test - low_t_test,
                             cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                             cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_t_test", 
                "cov_dw", "cov_t_test")

########## 

var_50_fin = var_50 %>% mutate(width_dw = up_dw - low_dw,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                               cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_t_test", 
                "cov_dw", "cov_t_test")
########## 

var_100_fin = var_100 %>% mutate(width_dw = up_dw - low_dw,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                               cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_t_test", 
                "cov_dw", "cov_t_test")
########## 

var_250_fin = var_250 %>% mutate(width_dw = up_dw - low_dw,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                               cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_t_test", 
                "cov_dw", "cov_t_test")

########## 
var_500_fin = var_500 %>% mutate(width_dw = up_dw - low_dw,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0 & low_dw <0, 1, 0),
                               cov_t_test = ifelse(up_t_test >0 & low_t_test <0, 1, 0)) %>%
  dplyr::select("width_dw", "width_t_test", 
                "cov_dw", "cov_t_test")
########## 

n2_seq = c(20,50,100,250,500)

fin_var = cbind(
  apply( var_20_fin, 2, mean),
  apply( var_50_fin, 2, mean),
  apply( var_100_fin, 2, mean),
  apply( var_250_fin, 2, mean),
  apply( var_500_fin, 2, mean)
)

colnames(fin_var) = n2_seq
