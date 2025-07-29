setwd("~/Documents/Projects/conf_tail/discrete/fast_sims/hydras")

library(ggplot2)
source("plot_x2.R")
source("plot_xlogx.R")
source("plot_bdd.R")
source("plot_subg.R")
source("plot_var.R")


n2_seq = c(20,50,100,250,500)


fin1 = expand.grid(n = n2_seq, Method = c("CDF", "Wald"),
                   Property = c("Width", "Coverage"), 
                   Type = c("Bounded", "2-moment","Var", "Sub-G", "|X|log|X|")) 

# layout( rbind(matrix(c(1:6), ncol = 3, byrow = F),
#               rep(7,3)), heights = c(4,4,0.7) )

fin1 = mutate(fin1, Val = NA)
fin_bdd

#bdd
fin1[1:10,5] = c(fin_bdd[1,] , fin_bdd[4,])   
fin1[11:20 ,5] = c(fin_bdd[5,], fin_bdd[8,]) 

## second mom
fin1[1:10 + 20,5] = c(fin_x2[1,], fin_x2[4,])     
fin1[11:20 + 20,5] = c(fin_x2[5,], fin_x2[8,])

## var
fin1[1:10 + 40,5] = c(fin_var[1,], fin_var[2,])     
fin1[11:20 + 40,5] = c(1 - fin_var[3,], 1 - fin_x2[4,])

## SUB-G
fin1[1:10 + 60,5] = c(fin_subg[1,] , fin_subg[4,] )   
fin1[11:20 + 60,5] = c(fin_subg[5,], fin_subg[8,])

## xlogx
fin1[1:10 + 80,5] = c(fin_xlogx[1,], fin_xlogx[4,] )    
fin1[11:20 + 80,5] = c(fin_xlogx[5,], fin_xlogx[8,])


fin1


ggplot(fin1, 
       aes(x = log(n), y = Val, col = Method)) +
  # facet_wrap(~Type + Property) 
  geom_hline(data = data.frame(xint=.90, Property ='Coverage'),
             aes(yintercept = xint), linetype = 'dashed', color = 'gray41') +
  
  geom_line() +
  geom_point() + 
  
  facet_grid(cols = vars(Type), rows = vars(Property), scales = 'free_y') +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'bottom',
        strip.background = element_rect(fill = 'white'), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggplot2::theme_set(theme_bw() + 
                     theme(strip.background = element_rect(fill = 'white'), 
                           panel.grid.minor = element_blank()))

# ggplot( fin1 %>% filter(Property == "Width") ) +
#   geom_point( aes(x = log(n), y = Val, col = Method)) + 
#   geom_line(aes(x = log(n), y = Val, col = Method)) + 
#   facet_wrap(~Type) 



# 
# 
# #par(mfrow = c(4,2))
# par(mai= c(rep(0.6,2),rep(0.5,2)))
# 
# plot(log(n2_seq), fin_bdd[1,], type = "l", ylim = c(min(fin_bdd[c(1,2,4),]), 
#                                                     max(fin_bdd[c(1,2,4),])),
#      xlab ="log(n)", ylab = "Width", main = "Width(CI) for bounded r.v.", 
#      cex.lab = 1.2, lwd = 1)
# points(log(n2_seq), fin_bdd[1,], lwd = 1)
# points(log(n2_seq), fin_bdd[4,], type = "l", col = "red") 
# points(log(n2_seq), fin_bdd[4,], col = "red") 
# 
# par(mai= c(rep(0.6,2),rep(0.5,2)))
# 
# plot(log(n2_seq), fin_bdd[5,], ylim = c(0.75,1.1),
#      xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for bounded r.v.",
#      cex.lab = 1.2, lwd = 1)
# points(log(n2_seq), fin_bdd[8,], col = "red") 
# abline(h = 0.9, lty = 2, col = 1)
# 
# par(mai= c(rep(0.6,2),rep(0.5,2)))
# 
# plot(log(n2_seq), fin_x2[1,], type = "l", ylim = c(min(fin_x2[c(1,2,4),]), max(fin_x2[c(1,2,4),])),
#      xlab ="log(n)", ylab = "Width", main = "Width(CI) for bounded second moment", cex.lab = 1.2)
# points(log(n2_seq), fin_x2[1,])
# points(log(n2_seq), fin_x2[4,], type = "l", col = "red") 
# points(log(n2_seq), fin_x2[4,], col = "red") 
# 
# plot(log(n2_seq), fin_x2[5,], ylim = c(0.75,1.1),
#      xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for bounded second moment", 
#      cex.lab = 1.2, lwd = 1)
# points(log(n2_seq), fin_x2[8,], col = "red") 
# abline(h = 0.9, lty = 2, col = 1)
# 
# # plot(log(n2_seq), fin_subg[1,], type = "l", ylim = c(min(fin_subg[c(1,2,4),]), max(fin_subg[c(1,2,4),])),
# #      xlab ="log(n)", ylab = "Width", main = "Width(CI) for sub-gaussian r.v.", cex.lab = 1.2)
# # points(log(n2_seq), fin_subg[1,])
# # points(log(n2_seq), fin_subg[4,], type = "l", col = "red")
# # points(log(n2_seq), fin_subg[4,], col = "red")
# # 
# # 
# # 
# # plot(log(n2_seq), fin_subg[5,], ylim = c(0.75,1.1),
# #      xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for sub-gaussian r.v.",
# #      cex.lab = 1.2, lwd = 1)
# # points(log(n2_seq), fin_subg[8,], col = "red") 
# # abline(h = 0.9, lty = 2, col = 1)
# # 
# # par(mai= c(rep(0.6,2),rep(0.5,2)))
# 
# plot(log(n2_seq), fin_xlogx[1,], type = "l", ylim = c(min(fin_xlogx[c(1,2,4),]), max(fin_xlogx[c(1,2,4),])),
#      xlab ="log(n)", ylab = "Width", main = "Width(CI) for bounded E[|X|log|X|]", cex.lab = 1.2)
# points(log(n2_seq), fin_xlogx[1,])
# points(log(n2_seq), fin_xlogx[4,], type = "l", col = "red") 
# points(log(n2_seq), fin_xlogx[4,], col = "red")
# 
# plot(log(n2_seq), fin_xlogx[5,], ylim = c(0.75,1.1),
#      xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for bounded E[|X|log|X|]",
#      cex.lab = 1.2, lwd = 1)
# points(log(n2_seq), fin_xlogx[8,], col = "red") 
# abline(h = 0.9, lty = 2, col = 1)
# 
# 
# par(mai=c(0,0,0,0), cex = 0.8)
# plot.new()
# legend(x = "center",  horiz=T, fill = c("black", "red"),
#        legend = c("CDF method", "Wald CI"),  bty = "n")



