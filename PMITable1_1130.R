# Load necessary libraries
library(pracma)
library(MASS)

# Define functions

# Function to generate data
data_generate <- function(tm, p, alpha, V) {
  # V: d x p matrix, p = [log(d)^2]
  z <- rep(0, p)
  to_seq <- rep(0, tm)
  r <- matrix(rnorm(tm * p, sd = 1 / sqrt(p)), tm, p)
  for (i in 1:tm) {
    z <- sqrt(alpha) * z + sqrt(1 - alpha) * r[i, ]  # p x 1
    se <- exp(V %*% as.matrix(z))  # d x 1
    to_seq[i] <- sample(1:d, 1, replace = TRUE, prob = se / sum(se))
  }
  return(to_seq)
}

# Function to calculate co-occurrence matrix for one patient, returns first f rows
coocur_cal_row <- function(wd_seq, d, q = 1, f) {
  co_m <- matrix(0, d, d)
  l <- length(wd_seq)
  for (i in 1:(l - q)) {
    for (j in 1:min(q, l - i)) {
      a <- wd_seq[i]
      b <- wd_seq[i + j]
      if (a == b) {
        co_m[a, a] <- co_m[a, a] + 2
      } else {
        co_m[a, b] <- co_m[a, b] + 1
        co_m[b, a] <- co_m[b, a] + 1
      }
    }
  }
  co_m <- co_m[1:f, ]
  return(co_m)
}

# Function to calculate co-occurrence matrix for one patient
coocur_cal <- function(wd_seq, d, q = 1) {
  co_m <- matrix(0, d, d)
  l <- length(wd_seq)
  for (i in 1:(l - q)) {
    for (j in 1:min(q, l - i)) {
      a <- wd_seq[i]
      b <- wd_seq[i + j]
      if (a == b) {
        co_m[a, a] <- co_m[a, a] + 2
      } else {
        co_m[a, b] <- co_m[a, b] + 1
        co_m[b, a] <- co_m[b, a] + 1
      }
    }
  }
  return(co_m)
}

# Function to compute SPPMI entry
SPPMI_calc_entry <- function(co, q = 1, g, len, n) {
  # co : f x d co-occurrence matrix
  r <- rowSums(co)
  f <- length(r)
  if (g > f) stop("g should be less than or equal to f!")
  c <- (2 * q * len - 2 * q^2) * n  # Total co-occurrence counts
  SPPMI <- rep(0, f)
  for (j in 1:f) {
    SPPMI[j] <- log(max(1e-3, c * co[g, j] / (r[g] * r[j])))
  }
  return(SPPMI)
}

# Function to compute SPPMI matrix
SPPMI_calc <- function(co, q = 1, d, len, n) {
  r <- rowSums(co)
  c <- (2 * q * len - 2 * q^2) * n  # Total co-occurrence counts
  SPPMI <- matrix(0, d, d)
  for (i in 1:d) {
    for (j in 1:i) {
      SPPMI[i, j] <- log(max(1e-3, c * co[i, j] / (r[i] * r[j])))
      SPPMI[j, i] <- SPPMI[i, j]
    }
  }
  return(SPPMI)
}

# Function to compute low-rank SPPMI approximation
SPPMI_calc_lr <- function(co, q = 1, d, len, n, p) {
  SPPMI <- SPPMI_calc(co, q, d, len, n)
  svd_result <- svd(SPPMI, nu = p, nv = p)
  L <- svd_result$d[1:p]
  PMI_lr <- svd_result$u %*% diag(L) %*% t(svd_result$v)
  return(PMI_lr)
}

# Set parameters
set.seed(1)
d <- 100                # Number of words
n <- 2000               # Number of patients
q <- 5                  # Window size
c_len <- 800            # Length of time series
f <- d                  # Number of top words to consider
T1 <- c_len * q - q^2

# Initialize result data frames
typeI_results <- data.frame()
typeI_patient_level_results <- data.frame()

# File names including 'n' and 'd'
filename_typeI <- paste0("typeI_results_n", n, "_d", d, ".csv")
filename_typeI_patient <- paste0("typeI_patient_level_results_n", n, "_d", d, ".csv")

# Loop over g from 1 to 10
for (g in 1:100) {
  cat("Processing g =", g, "\n")
  
  # Use tryCatch to handle unexpected errors and save progress
  tryCatch({
    # Re-initialize variables for each g
    CO_patient_level <- matrix(0, n, f)
    CO_patient_level_rowsum <- matrix(0, n, f)
    co_null <- matrix(0, d, d)
    
    # Generate data under the null distribution
    total_seq <- matrix(sample(1:d, n * c_len, replace = TRUE, prob = rep(1 / d, d)), n, c_len)
    
    # Compute co-occurrence matrices
    if (n == 1) {
      tmp_seq <- total_seq
      co_entry <- coocur_cal_row(tmp_seq, d, q, f)
    } else {
      for (j in 1:n) {
        tmp_seq <- total_seq[j, ]
        co_row <- coocur_cal_row(tmp_seq, d, q, f)
        co_full <- coocur_cal(tmp_seq, d, q)
        co_null <- co_null + co_full
        CO_patient_level[j, ] <- co_row[g, 1:f]
        CO_patient_level_rowsum[j, ] <- rowSums(co_row)
      }
    }
    
    # Prepare variables for variance calculations
    co_tmp_list <- vector("list", f)
    for (i in 1:f) {
      co_tmp_list[[i]] <- matrix(0, n, d)
    }
    
    co_entry_mean <- co_null[1:f, ] / n  # f x d matrix
    co_rowsum_mean <- rowSums(co_null) / n
    
    for (j in 1:n) {
      tmp_seq <- total_seq[j, ]
      co_full <- coocur_cal(tmp_seq, d, q)
      for (i in 1:f) {
        co_tmp_list[[i]][j, ] <- co_full[i, ] / co_entry_mean[i, ] - sum(co_full[i, ]) / co_rowsum_mean[i] - rowSums(co_full) / co_rowsum_mean
      }
    }
    
    # Compute standard deviations for PMI entries
    PMI_entry_sd_co <- rep(0, f)
    for (i in 1:f) {
      PMI_entry_sd_co[i] <- sd(co_tmp_list[[g]][, i]) / sqrt(n)
    }
    
    p_w <- rowSums(co_null) / (2 * n * T1)
    
    # Compute theoretical variances under null distribution
    var_g <- rep(0, f)
    for (i in 1:f) {
      if (i == g) {
        var_g[i] <- (1 - p_w[i])^2 / (n * T1 * p_w[i]^2)
      } else {
        var_g[i] <- (1 - p_w[i] - p_w[g] + 2 * p_w[i] * p_w[g]) / (2 * n * T1 * p_w[i] * p_w[g])
      }
    }
    
    # Monte Carlo method to estimate Type I error rates
    m <- 100   # Number of repetitions
    PMI_entry_list <- matrix(0, m, f)
    typeI_patient_level <- rep(0, f)
    typeI <- rep(0, f)
    
    for (k in 1:m) {
      total_seq <- matrix(sample(1:d, n * c_len, replace = TRUE, prob = rep(1 / d, d)), n, c_len)
      co_etr <- matrix(0, d, d)
      if (n == 1) {
        tmp_seq <- total_seq
        co_full <- coocur_cal(tmp_seq, d, q)
        co_etr <- co_full
      } else {
        for (j in 1:n) {
          tmp_seq <- total_seq[j, ]
          co_full <- coocur_cal(tmp_seq, d, q)
          co_etr <- co_etr + co_full
        }
      }
      
      if (0 %in% rowSums(co_etr)) {
        stop("Invalid co-occurrence matrix!")
      }
      
      co_etr_f <- co_etr[1:f, 1:d]
      SP_est <- SPPMI_calc_entry(co_etr_f, q, g, c_len, n)
      PMI_entry_list[k, ] <- SP_est
      
      for (j in 1:f) {
        # Compare with theoretical variance
        if (abs(PMI_entry_list[k, j]) > sqrt(var_g[j]) * qnorm(0.975)) {
          typeI[j] <- typeI[j] + 1
        }
        # Compare with patient-level variance
        if (abs(PMI_entry_list[k, j]) > PMI_entry_sd_co[j] * qnorm(0.975)) {
          typeI_patient_level[j] <- typeI_patient_level[j] + 1
        }
      }
    }
    
    # Compute Type I error rates
    typeI <- typeI / m
    typeI_patient_level <- typeI_patient_level / m
    
    # Collect results with 'g' as a column
    typeI_results <- rbind(typeI_results, data.frame(g = g, i = 1:f, TypeI = typeI))
    typeI_patient_level_results <- rbind(typeI_patient_level_results, data.frame(g = g, i = 1:f, TypeI_Patient_Level = typeI_patient_level))
    
    # Save results every time g %% 5 == 0
    if (g %% 5 == 0) {
      write.csv(typeI_results, file = filename_typeI, row.names = FALSE)
      write.csv(typeI_patient_level_results, file = filename_typeI_patient, row.names = FALSE)
      cat("Results saved for g =", g, "\n")
    }
  }, error = function(e) {
    cat("Error at g =", g, ":", conditionMessage(e), "\n")
    # Save results before exiting
    write.csv(typeI_results, file = filename_typeI, row.names = FALSE)
    write.csv(typeI_patient_level_results, file = filename_typeI_patient, row.names = FALSE)
    stop("Terminating due to error.")
  })
}


write.csv(typeI_results, file = filename_typeI, row.names = FALSE)
write.csv(typeI_patient_level_results, file = filename_typeI_patient, row.names = FALSE)
cat("Final results saved.\n")
