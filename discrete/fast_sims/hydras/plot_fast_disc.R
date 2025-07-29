setwd("~/Documents/Projects/conf_tail/discrete/fast_sims/hydras")

source("plot_x2.R")
source("plot_xlogx.R")
source("plot_bdd.R")
source("plot_subg.R")

n2_seq = c(20,50,100,250,500)
layout( rbind(matrix(c(1:8), ncol = 4, byrow = F), rep(9,4)), heights = c(4,4,0.7))

#par(mfrow = c(4,2))
par(mai= c(rep(0.6,2),rep(0.5,2)))

plot(log(n2_seq), fin_bdd[1,], type = "l", ylim = c(min(fin_bdd[c(1,2,4),]), 
                                                    max(fin_bdd[c(1,2,4),])),
     xlab ="log(n)", ylab = "Width", main = "Width(CI) for bounded r.v.", 
     cex.lab = 1.2, lwd = 2.5)
points(log(n2_seq), fin_bdd[1,], lwd = 2.5)
points(log(n2_seq), fin_bdd[2,], type = "l", col = "red")
points(log(n2_seq), fin_bdd[2,], col = "red")
points(log(n2_seq), fin_bdd[3,], type = "l", col = "green")
points(log(n2_seq), fin_bdd[3,], col = "green")
points(log(n2_seq), fin_bdd[4,], type = "l", col = "blue") 
points(log(n2_seq), fin_bdd[4,], col = "blue") 

par(mai= c(rep(0.6,2),rep(0.5,2)))

plot(log(n2_seq), fin_bdd[5,], ylim = c(0.75,1.1),
     xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for bounded r.v.",
     cex.lab = 1.2, lwd = 4)
points(log(n2_seq), fin_bdd[6,], col = "red", lwd = 2)
points(log(n2_seq), fin_bdd[7,], col = "green")
points(log(n2_seq), fin_bdd[8,], col = "blue") 
abline(h = 0.9, lty = 2, col = 1)

par(mai= c(rep(0.6,2),rep(0.5,2)))

plot(log(n2_seq), fin_x2[1,], type = "l", ylim = c(min(fin_x2[c(1,2,4),]), max(fin_x2[c(1,2,4),])),
     xlab ="log(n)", ylab = "Width", main = "Width(CI) for bounded second moment", cex.lab = 1.2)
points(log(n2_seq), fin_x2[1,])
points(log(n2_seq), fin_x2[2,], type = "l", col = "red")
points(log(n2_seq), fin_x2[2,], col = "red")
points(log(n2_seq), fin_x2[3,], type = "l", col = "green")
points(log(n2_seq), fin_x2[3,], col = "green")
points(log(n2_seq), fin_x2[4,], type = "l", col = "blue") 
points(log(n2_seq), fin_x2[4,], col = "blue") 

plot(log(n2_seq), fin_x2[5,], ylim = c(0.75,1.1),
     xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for bounded second moment", 
     cex.lab = 1.2, lwd = 4)
points(log(n2_seq), fin_x2[6,], col = "red", lwd = 2)
points(log(n2_seq), fin_x2[7,], col = "green")
points(log(n2_seq), fin_x2[8,], col = "blue") 
abline(h = 0.9, lty = 2, col = 1)

plot(log(n2_seq), fin_subg[1,], type = "l", ylim = c(min(fin_subg[c(1,2,4),]), max(fin_subg[c(1,2,4),])),
     xlab ="log(n)", ylab = "Width", main = "Width(CI) for sub-gaussian r.v.", cex.lab = 1.2)
points(log(n2_seq), fin_subg[1,])
points(log(n2_seq), fin_subg[2,], type = "l", col = "red")
points(log(n2_seq), fin_subg[2,], col = "red")
points(log(n2_seq), fin_subg[3,], type = "l", col = "green")
points(log(n2_seq), fin_subg[3,], col = "green")
points(log(n2_seq), fin_subg[4,], type = "l", col = "blue")
points(log(n2_seq), fin_subg[4,], col = "blue")



plot(log(n2_seq), fin_subg[5,], ylim = c(0.75,1.1),
     xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for sub-gaussian r.v.",
     cex.lab = 1.2, lwd = 4)
points(log(n2_seq), fin_subg[6,], col = "red", lwd = 2)
points(log(n2_seq), fin_subg[7,], col = "green")
points(log(n2_seq), fin_subg[8,], col = "blue") 
abline(h = 0.9, lty = 2, col = 1)

par(mai= c(rep(0.6,2),rep(0.5,2)))

plot(log(n2_seq), fin_xlogx[1,], type = "l", ylim = c(min(fin_xlogx[c(1,2,4),]), max(fin_xlogx[c(1,2,4),])),
     xlab ="log(n)", ylab = "Width", main = "Width(CI) for bounded E[|X|log|X|]", cex.lab = 1.2)
points(log(n2_seq), fin_xlogx[1,])
points(log(n2_seq), fin_xlogx[2,], type = "l", col = "red")
points(log(n2_seq), fin_xlogx[2,], col = "red")
points(log(n2_seq), fin_xlogx[3,], type = "l", col = "green")
points(log(n2_seq), fin_xlogx[3,], col = "green")
points(log(n2_seq), fin_xlogx[4,], type = "l", col = "blue") 
points(log(n2_seq), fin_xlogx[4,], col = "blue")

plot(log(n2_seq), fin_xlogx[5,], ylim = c(0.75,1.1),
     xlab ="log(n)", ylab = "Coverage", main = "Coverage(CI) for bounded E[|X|log|X|]",
     cex.lab = 1.2, lwd = 4)
points(log(n2_seq), fin_xlogx[6,], col = "red", lwd = 2)
points(log(n2_seq), fin_xlogx[7,], col = "green")
points(log(n2_seq), fin_xlogx[8,], col = "blue") 
abline(h = 0.9, lty = 2, col = 1)


par(mai=c(0,0,0,0), cex = 0.8)
plot.new()
legend(x = "center",  horiz=T, fill = c("black", "blue", "red", "green"),
      legend = c("CDF method", "Wald CI", "CDF method \n (range of data)", "CDF method \n (range + tail)"),  bty = "n")



