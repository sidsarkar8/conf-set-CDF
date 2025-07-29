library(dplyr)

# Set working directory
setwd("~/Documents/Projects/conf_tail/other_methods/hydra")

# Set true mean for coverage evaluation
mu2 <- 3

# Read in x2 datasets (excluding 500 case)
x2_20  = read.csv("x2/x2_20.csv")[,-1]
x2_50  = read.csv("x2/x2_50.csv")[,-1]
x2_100 = read.csv("x2/x2_100.csv")[,-1]
x2_250 = read.csv("x2/x2_250.csv")[,-1]

# Function to compute widths and coverage
process_df <- function(df) {
  df %>% mutate(
    width_dw = up_dw - low_dw,
    width_t_test = up_t_test - low_t_test,
    width_catoni = up_catoni - low_catoni,
    width_cheby = up_cheby - low_cheby,
    
    cov_dw = ifelse(up_dw > mu2 & low_dw < mu2, 1, 0),
    cov_t_test = ifelse(up_t_test > mu2 & low_t_test < mu2, 1, 0),
    cov_catoni = ifelse(up_catoni > mu2 & low_catoni < mu2, 1, 0),
    cov_cheby = ifelse(up_cheby > mu2 & low_cheby < mu2, 1, 0)
  ) %>% select(starts_with("width_"), starts_with("cov_"))
}

# Process datasets
x2_20_fin  = process_df(x2_20)
x2_50_fin  = process_df(x2_50)
x2_100_fin = process_df(x2_100)
x2_250_fin = process_df(x2_250)

# Combine results into a matrix
n2_seq = c(20, 50, 100, 250)
fin_x2 = cbind(
  apply(x2_20_fin, 2, mean),
  apply(x2_50_fin, 2, mean),
  apply(x2_100_fin, 2, mean),
  apply(x2_250_fin, 2, mean)
)
colnames(fin_x2) = n2_seq

# Plot: Width comparison for x2
plot(log(n2_seq), fin_x2["width_dw",], type = "l",
     ylim = range(fin_x2[c("width_dw", "width_t_test", "width_catoni"),]),
     xlab = "log(n)", ylab = "Width",
     main = "Width comparison for x2 centered case (mu = 3)")

lines(log(n2_seq), fin_x2["width_t_test",], col = "blue")
lines(log(n2_seq), fin_x2["width_catoni",], col = "red")

legend("topright", legend = c("CDF method", "t-test", "Catoni"),
       col = c("black", "blue", "red"), lty = 1)

# Plot: Coverage comparison for x2
plot(log(n2_seq), fin_x2["cov_dw",], type = "l",
     ylim = c(0.8, 1),
     xlab = "log(n)", ylab = "Coverage",
     main = "Coverage comparison for x2 centered case (mu = 3)")

lines(log(n2_seq), fin_x2["cov_t_test",], col = "blue")
lines(log(n2_seq), fin_x2["cov_catoni",], col = "red")
abline(h = 0.9, lty = 2, col = "gray")  # Nominal reference line

legend("bottomright", legend = c("CDF method", "t-test", "Catoni", "Target = 0.9"),
       col = c("black", "blue", "red", "gray"), lty = c(1, 1, 1, 2))
