library(dplyr)

# Set working directory
setwd("~/Documents/Projects/conf_tail/other_methods/hydra")

# Set true mean for coverage calculations
mu1 <- 0.5

# Read data
bdd_20  = read.csv("bdd/bdd_20.csv")[,-1]
bdd_50  = read.csv("bdd/bdd_50.csv")[,-1]
bdd_100 = read.csv("bdd/bdd_100.csv")[,-1]
bdd_250 = read.csv("bdd/bdd_250.csv")[,-1]
bdd_500 = read.csv("bdd/bdd_500.csv")[,-1]

# Unified processing function
process_df <- function(df, mu) {
  df %>%
    mutate(
      width_dw = up_dw - low_dw,
      width_dw_x1_to_xn = up_dw_x1_to_xn - low_dw_x1_to_xn,
      width_dw_x1_to_xn_tails = up_dw_x1_to_xn_tails - low_dw_x1_to_xn_tails,
      width_t_test = up_t_test - low_t_test,
      width_anderson = up_anderson - low_anderson,
      width_cheby = up_cheby - low_cheby,
      width_emp_bern = up_emp_bern - low_emp_bern,
      width_hoeffding = up_hoeffding - low_hoeffding,
      width_PI_emp_bern = up_PI_emp_bern - low_PI_emp_bern,
      width_PI_hedged = up_PI_hedged - low_PI_hedged,
      
      cov_dw = ifelse(up_dw > mu & low_dw < mu, 1, 0),
      cov_dw_x1_to_xn = ifelse(up_dw_x1_to_xn > mu & low_dw_x1_to_xn < mu, 1, 0),
      cov_dw_x1_to_xn_tails = ifelse(up_dw_x1_to_xn_tails > mu & low_dw_x1_to_xn_tails < mu, 1, 0),
      cov_t_test = ifelse(up_t_test > mu & low_t_test < mu, 1, 0),
      cov_anderson = ifelse(up_anderson > mu & low_anderson < mu, 1, 0),
      cov_cheby = ifelse(up_cheby > mu & low_cheby < mu, 1, 0),
      cov_emp_bern = ifelse(up_emp_bern > mu & low_emp_bern < mu, 1, 0),
      cov_hoeffding = ifelse(up_hoeffding > mu & low_hoeffding < mu, 1, 0),
      cov_PI_emp_bern = ifelse(up_PI_emp_bern > mu & low_PI_emp_bern < mu, 1, 0),
      cov_PI_hedged = ifelse(up_PI_hedged > mu & low_PI_hedged < mu, 1, 0)
    ) %>%
    select(starts_with("width_"), starts_with("cov_"))
}

# Process datasets
bdd_20_fin  = process_df(bdd_20, mu1)
bdd_50_fin  = process_df(bdd_50, mu1)
bdd_100_fin = process_df(bdd_100, mu1)
bdd_250_fin = process_df(bdd_250, mu1)
bdd_500_fin = process_df(bdd_500, mu1)

# Combine results
n2_seq = c(20, 50, 100, 250, 500)
fin_bdd = cbind(
  apply(bdd_20_fin, 2, mean),
  apply(bdd_50_fin, 2, mean),
  apply(bdd_100_fin, 2, mean),
  apply(bdd_250_fin, 2, mean),
  apply(bdd_500_fin, 2, mean)
)
colnames(fin_bdd) = n2_seq

# Plot width comparison
plot(log(n2_seq), fin_bdd["width_dw",], type = "l",
     ylim = range(fin_bdd[c("width_dw", "width_t_test", "width_PI_hedged"),]),
     xlab = "log(n)", ylab = "Width",
     main = "Width comparison for bounded case [0,1]")

lines(log(n2_seq), fin_bdd["width_PI_hedged",], col = "red")
lines(log(n2_seq), fin_bdd["width_t_test",], col = "blue")

legend("topright", legend = c("CDF method", "Hedged", "t-test"),
       col = c("black", "red", "blue"), lty = 1)

# Plot coverage comparison
plot(log(n2_seq), fin_bdd["cov_dw",], type = "l",
     ylim = c(0.8, 1),
     xlab = "log(n)", ylab = "Coverage",
     main = "Coverage comparison for bounded case [0,1]")

lines(log(n2_seq), fin_bdd["cov_PI_hedged",], col = "red")
lines(log(n2_seq), fin_bdd["cov_t_test",], col = "blue")

abline(h = 0.9, lty = 2, col = "gray")  # reference line

legend("bottomright", legend = c("CDF method", "Hedged", "t-test", "Target = 0.9"),
       col = c("black", "red", "blue", "gray"), lty = c(1, 1, 1, 2))
