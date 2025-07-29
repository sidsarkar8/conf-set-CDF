library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/conf_tail/discrete/fast_sims/hydras/")

subg_20 = read.csv("subg/subg_20.csv")[,-1]
subg_50 = read.csv("subg/subg_50.csv")[,-1]
subg_100 = read.csv("subg/subg_100.csv")[,-1]
subg_250 = read.csv("subg/subg_250.csv")[,-1]
subg_500 = read.csv("subg/subg_500.csv")[,-1]

nrow(subg_20)
nrow(subg_50)
nrow(subg_100)
nrow(subg_250)
nrow(subg_500)

########## 

subg_20_fin = subg_20 %>% mutate(width_dw = up_dw - low_dw,
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

subg_50_fin = subg_50  %>% mutate(width_dw = up_dw - low_dw,
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

subg_100_fin = subg_100  %>% mutate(width_dw = up_dw - low_dw,
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

subg_250_fin = subg_250 %>% mutate(width_dw = up_dw - low_dw,
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

subg_100_fin = subg_100  %>% mutate(width_dw = up_dw - low_dw,
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

subg_250_fin = subg_250  %>% mutate(width_dw = up_dw - low_dw,
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


subg_500_fin = subg_500  %>% mutate(width_dw = up_dw - low_dw,
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

fin_subg = cbind(
  apply( subg_20_fin, 2, mean),
  apply( subg_50_fin, 2, mean),
  apply( subg_100_fin, 2, mean),
  apply( subg_250_fin, 2, mean),
  apply( subg_500_fin, 2, mean)
)

colnames(fin_subg) = n2_seq


library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/conf_tail/discrete/fast_sims/hydras/")

bdd_20 = read.csv("bdd/bdd_20.csv")[,-1]
bdd_50 = read.csv("bdd/bdd_50.csv")[,-1]
bdd_100 = read.csv("bdd/bdd_100.csv")[,-1]
bdd_250 = read.csv("bdd/bdd_250.csv")[,-1]
bdd_500 = read.csv("bdd/bdd_500.csv")[,-1]

nrow(bdd_20)
nrow(bdd_50)
nrow(bdd_100)
nrow(bdd_250)
nrow(bdd_500)

########## 

bdd_20_fin = bdd_20 %>% mutate(width_dw = up_dw - low_dw,
                               width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                 low_dw_x1_to_xn,
                               width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                 low_dw_x1_to_xn_tails,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0.5 & low_dw <0.5, 1, 0),
                               cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0.5 & low_dw_x1_to_xn <0.5, 1, 0),
                               cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0.5 & low_dw_x1_to_xn_tails <0.5, 1, 0),
                               cov_t_test = ifelse(up_t_test >0.5 & low_t_test <0.5, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 

bdd_50_fin = bdd_50 %>% mutate(width_dw = up_dw - low_dw,
                               width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                 low_dw_x1_to_xn,
                               width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                 low_dw_x1_to_xn_tails,
                               width_t_test = up_t_test - low_t_test,
                               cov_dw = ifelse(up_dw >0.5 & low_dw <0.5, 1, 0),
                               cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0.5 & low_dw_x1_to_xn <0.5, 1, 0),
                               cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0.5 & low_dw_x1_to_xn_tails <0.5, 1, 0),
                               cov_t_test = ifelse(up_t_test >0.5 & low_t_test <0.5, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")
########## 

bdd_100_fin = bdd_100 %>% mutate(width_dw = up_dw - low_dw,
                                 width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                   low_dw_x1_to_xn,
                                 width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                   low_dw_x1_to_xn_tails,
                                 width_t_test = up_t_test - low_t_test,
                                 cov_dw = ifelse(up_dw >0.5 & low_dw <0.5, 1, 0),
                                 cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0.5 & low_dw_x1_to_xn <0.5, 1, 0),
                                 cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0.5 & low_dw_x1_to_xn_tails <0.5, 1, 0),
                                 cov_t_test = ifelse(up_t_test >0.5 & low_t_test <0.5, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

apply( bdd_100_fin, 2, mean)

########## 

bdd_250_fin = bdd_250 %>% mutate(width_dw = up_dw - low_dw,
                                 width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                   low_dw_x1_to_xn,
                                 width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                   low_dw_x1_to_xn_tails,
                                 width_t_test = up_t_test - low_t_test,
                                 cov_dw = ifelse(up_dw >0.5 & low_dw <0.5, 1, 0),
                                 cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0.5 & low_dw_x1_to_xn <0.5, 1, 0),
                                 cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0.5 & low_dw_x1_to_xn_tails <0.5, 1, 0),
                                 cov_t_test = ifelse(up_t_test >0.5 & low_t_test <0.5, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")

########## 

bdd_500_fin = bdd_500 %>% mutate(width_dw = up_dw - low_dw,
                                 width_dw_x1_to_xn = up_dw_x1_to_xn - 
                                   low_dw_x1_to_xn,
                                 width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - 
                                   low_dw_x1_to_xn_tails,
                                 width_t_test = up_t_test - low_t_test,
                                 cov_dw = ifelse(up_dw >0.5 & low_dw <0.5, 1, 0),
                                 cov_dw_x1_to_xn =  ifelse(up_dw_x1_to_xn >0.5 & low_dw_x1_to_xn <0.5, 1, 0),
                                 cov_dw_x1_to_xn_tails =  ifelse(up_dw_x1_to_xn_tails >0.5 & low_dw_x1_to_xn_tails <0.5, 1, 0),
                                 cov_t_test = ifelse(up_t_test >0.5 & low_t_test <0.5, 1, 0)) %>%
  dplyr::select("width_dw", "width_dw_x1_to_xn", "width_dw_x1_to_xn_tails", "width_t_test", 
                "cov_dw", "cov_dw_x1_to_xn", "cov_dw_x1_to_xn_tails", "cov_t_test")
########## 

n2_seq = c(20,50,100,250,500)

fin_bdd = cbind(
  apply( bdd_20_fin, 2, mean),
  apply( bdd_50_fin, 2, mean),
  apply( bdd_100_fin, 2, mean),
  apply( bdd_250_fin, 2, mean),
  apply( bdd_500_fin, 2, mean)
)

colnames(fin_bdd) = n2_seq



# plot(log(n2_seq), fin_subg[1,], type = "l", ylim = c(min(fin_subg[c(1,2,4),]), max(fin_subg[c(1,2,4),])),
#      xlab = "log(n)", ylab = "Width", main = "Width comparison for sub-gaussian case")
# points(log(n2_seq), fin_subg[2,], type = "l", col = "red")
# points(log(n2_seq), fin_subg[4,], type = "l", col = "blue")
#legend(x = 5, y =  0.5, fill = c("black", "red", "blue"),
#       legend = c("CDF method","CDF method (data)", "t-test"))
#plot(n2_seq, fin_bdd[3,])


