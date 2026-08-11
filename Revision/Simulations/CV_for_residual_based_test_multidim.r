# ==========================================================
# FINAL SIMULATION: Quantiles for residual-based test applied to empirical study
# Setup: d1 = 1, d2 = 2 (Fixed)
# ==========================================================
set.seed(12345) 
probs <- seq(0, 1, by = 0.001); d2_fixed <- 2; set.seed(9999); n_iter <- 100000 ; n_step <- 200 ; d1_list <- 1

## Critical values for the test H_0: d_v=3,2,1 can be obtained by setting d1_list to 3,2,1. Critical value for 5% significance level is reported at the end of this code.  
quantile_list <- list()
quantile_list2 <- list()
cat("Starting detailed simulation for d1 =", paste(d1_list, collapse=", "), "...\n")

for (d1 in d1_list) {
  d_total <- d1 + d2_fixed
  trace_dist <- numeric(n_iter)
  max_dist <- numeric(n_iter)
  for (i in 1:n_iter) {
    # 1. Generate d_total-dimensional BM and demean
    dw <- matrix(rnorm(n_step * d_total, mean=0, sd = sqrt(1/n_step)), n_step, d_total)
    w <- apply(dw, 2, cumsum)
    b <- sweep(w, 2, colMeans(w))
    
    b1 <- cbind(b[, 1:d1])      
    b2 <- cbind(b[, (d1+1):d_total]) 
    
    # 2. Calculate Residual BM (G): B1.2
    int_b2b2 <- t(b2) %*% b2 / n_step
    int_b2b1 <- t(b2) %*% b1 / n_step
    g <- b1 - b2 %*% solve(int_b2b2, int_b2b1)
    
    # 3. Calculate Integrated Residual BM (F): F = int G
    f <- apply(g, 2, function(x) cumsum(x) / n_step)
    
    # 4. Calculate Matrix Integrals and Statistic
    M_GG <- t(g) %*% g / n_step
    M_FF <- t(f) %*% f / n_step
    
    # Trace(inv(M_FF) %*% M_GG)
    trace_dist[i] <- sum(diag(solve(M_FF, M_GG)))
	max_dist[i] <- eigen(solve(M_FF, M_GG))$values[1]

  }
  
# Calculate and store quantiles at 0.1% intervals
  quantile_list[[paste0("d1_", d1)]] <- quantile(trace_dist, probs = probs)
  quantile_list2[[paste0("d1_", d1)]] <- quantile(max_dist, probs = probs)
  
# Check progress (print 5% CV as a sample)
  cat("Completed d1 =", d1, "| 5% CV =", round(quantile(trace_dist, 0.5), 2), "\n")
    cat("Completed d1 =", d1, "| 5% CV =", round(quantile(max_dist, 0.5), 2), "\n")
}

# Convert final results into a single large data frame (for easier viewing)
final_quantile_table <- as.data.frame(quantile_list)
final_quantile_table2 <- as.data.frame(quantile_list2)
rownames(final_quantile_table) <- paste0("p_", probs * 100, "%")
rownames(final_quantile_table2) <- paste0("p_", probs * 100, "%")
#cat("\n--- Detailed Quantile Table (Top 10 rows) ---\n")
final_quantile_table[951,]## 5% quantile for the trace test
#final_quantile_table2[951,]## 5% quantile

