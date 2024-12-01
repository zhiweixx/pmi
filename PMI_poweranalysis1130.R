arg = commandArgs(T)

d = as.numeric(arg[1])
ord = as.numeric(arg[2])
c_len = as.numeric(arg[3])

p <- floor(log(d)^2) + 1                                           #dimension of word vector
q <- 5                                                             #window size
f = d
alpha <- 1 - log(d)/(p^2)
kpa <- d^{-1} * ord 
g = 1
T1 = c_len*q - q^2

library(pracma)
library(MASS)
data_gene <- function(tm,p,alpha,V){
  #V: d*p, p=[log(d)^2]
  z <- rep(0,p)
  to_seq <- rep(0,tm)
  r <- matrix(rnorm(tm*p, sd = 1/sqrt(p)), tm, p)
  for(i in 1:tm){
    z <- sqrt(alpha)*z + sqrt(1-alpha)*r[i,]             #p*1  
    se <- exp(V %*% as.matrix(z))  #d*1
    to_seq[i] <- sample(1:d, 1, replace = T, prob = se/sum(se))
  }
  return(to_seq)
}
coocur_cal <- function(wd_seq, d ,q = 1){ #input one patient's time series data, output the co-occurrence matrix
  co_m <- matrix(0,d,d)
  l <- length(wd_seq)
  for (i in 1:(l-q)) {
    for (j in 1:min(q,l-i)) {
      a = wd_seq[i]
      b = wd_seq[i+j]
      if(a == b){
        co_m[a,a] <- co_m[a,a] + 2 
      }else{
        co_m[a,b] <- co_m[a,b] + 1
        co_m[b,a] <- co_m[b,a] + 1
      }
    }
  }
  return(co_m)
}
coocur_cal_row <- function(wd_seq, d ,q = 1,f){ #input one patient's EHR, output the first f rows of co-occurrence matrix
  co_m <- matrix(0,f,d)
  l <- length(wd_seq)
  for (i in 1:(l-q)) {
    for (j in 1:min(q,l-i)) {
      a = wd_seq[i]
      b = wd_seq[i+j]
      if ((a > f) & (b > f)) next
      if(a == b){
        co_m[a,a] <- co_m[a,a] + 2 
      }else{
        a_p = min(a,b)
        b_p = max(a,b)
        co_m[a_p,b_p] <- co_m[a_p,b_p] + 1
        if(b_p <= f) co_m[b_p, a_p] <- co_m[b_p, a_p] + 1
      }
    }
  }
  return(co_m)
}
SPPMI_calc_lr <- function(co, q = 1, d, len, n, p){ #compute PMI low rank estimator
  r <- rowSums(co)
  c <- (2*q*len - 2*q^2)*n #total co-occurrence
  SPPMI <- matrix(0,d,d)
  for (i in 1:d) {
    for(j in 1:i){
      SPPMI[i,j] = log(max(1e-3, c*co[i,j]/(r[i]*r[j])))
      SPPMI[j,i] = SPPMI[i,j]
    }
  }
  re <- svd(SPPMI, nu = p, nv = p)
  L = re$d[1:p]
  PMI_lr = re$u %*% diag(L) %*% t(re$v)
  return(PMI_lr)
}

SPPMI_calc <- function(co, q = 1, d, len, n){ #compute PMI low rank estimator
  r <- rowSums(co)
  c <- (2*q*len - 2*q^2)*n #total co-occurrence
  SPPMI <- matrix(0,d,d)
  for (i in 1:d) {
    for(j in 1:i){
      SPPMI[i,j] = log(max(1e-3, c*co[i,j]/(r[i]*r[j])))
      SPPMI[j,i] = SPPMI[i,j]
    }
  }
  return(SPPMI)
}

SPPMI_calc_entry <- function(co, q = 1, g, len, n){ #compute PMI[g, 1:f], g \leq f
  #co : f*d
  r <- rowSums(co)
  f <- length(r)
  if (g > f) cat("g should less than f!\n")
  c <- (2*q*len - 2*q^2)*n #total co-occurrence
  SPPMI <- rep(0,f)
  for (j in 1:f) {
    SPPMI[j] = log(max(1e-3, c*co[g,j]/(r[g]*r[j])))
  }
  return(SPPMI)
}

set.seed(1)
pi = rep(floor(p/3),3)
if(p%%3!=0) pi[1:(p%%3)] = pi[1:(p%%3)]+1
di = rep(floor(d/3),3)
if(d%%3!=0) di[1:(d%%3)] = di[1:(d%%3)]+1
di0 = cumsum(di) - di
di1 = d - cumsum(di)
U <- lapply(1:3, function(i){
  tmp = orth(randn(n=di[i],m=pi[i]+2))[,1:pi[i],drop=FALSE]
  return(rbind(matrix(0, nrow = di0[i], ncol = pi[i]),
               tmp,
               matrix(0, nrow = di1[i], ncol = pi[i])))
})
U = do.call("cbind", U)
L <- kpa*diag(p)
V <- U %*% L                                                       #V is d*p matrix
#estimate true PMI
long_len <- max(d^2,400^2)*1e3
EHR = data_gene(tm=long_len,p,alpha,V)
CO = coocur_cal(EHR, d ,q = q)
pmi_true = SPPMI_calc(co = CO, q, d, long_len, n=1)
pmi_lr_true = SPPMI_calc_lr(co = CO, q = q, d = d, len = long_len, n = 1, p=p)
pmi_true = matrix(0, d, d)
vec_p = apply(V, 1, function(x) exp(sum(x^2)/(2*p)))
vec_p = vec_p/sum(vec_p)
alpha_p = 2*sum(sapply(1:q, function(u) (c_len - u)*alpha^{u}))/(p*(2*q*c_len - 2*q^2))
V_tilde = sqrt(alpha_p) * (V - as.matrix(rep(1,d)) %*% (vec_p %*% V))
VVT = V_tilde %*% t(V_tilde)
rm(list = setdiff(ls(), c('V','pmi_true','pmi_lr_true','d','ord','n','c_len','arg')))

for(n in c(200,400,600,800,1000)){
  ##############################Generate n patients' word sequence and validation of theorem 3.5 ###########################################
  #########################estimate variance of PMI[g, 1:j] by cooccurence matrix & data splitting using non-low rank and low rank estimator##########################
  ############# use cooccurrence matrix to estimate the variance of PMI[g,1:f] without data splitting ############
  co_tmp_list <- vector("list",f)
  for(i in 1:f){
    co_tmp_list[[i]] <- matrix(0,n,d)
  }
  PMI_entry_sd_lr <- rep(0,f)                #use sqrt(n) PMI_lr to estimate the variance of PMI_lr[g,1:f]
  PMI_entry_sd_lr_est <- rep(0,f)            #use eq (2.6) in the paper to estimate the variance of PMI_lr[g,1:f]
  
  ############### data generation #################
  total_seq <- matrix(0, nrow = n, ncol = c_len)
  for(l in 1:n){                                                   
    total_seq[l,] <- data_gene(c_len,p,alpha,V)
  }
  co <- matrix(0,d,d)
  for (j in 1:n) {
    tmp2 <- total_seq[j,] 
    tmp3 <- coocur_cal(tmp2, d, q)
    co <-  co + tmp3
  }
  
  co_entry_mean <- co[1:f,]/n #f*d matrix
  co_rowsum_mean <- rowSums(co)/n
  SP_est <- SPPMI_calc(co, q, d, c_len, n)
  r_SP <- svd(SP_est, nu = p, nv = p)
  L_SP = r_SP$d[1:p]
  V_tmp = eigen(SP_est, symmetric = T)
  V_hat = sqrt(alpha_p)*V_tmp$vectors[,1:p] %*% diag(sqrt((V_tmp$values[1:p] * as.numeric(V_tmp$values[1:p]>0))))
  rm(V_tmp)
  
  
  ###################
  PMI_mean_est <- SP_est
  PMI_lr_mean_est <- SPPMI_calc_lr(co, q, d, c_len, n,p)
  
  P_SP = r_SP$u %*% diag(sign(L_SP)) %*% t(r_SP$v)
  tmp_P_g <- P_SP[g,]
  for(j in 1:n){
    tmp2 <- total_seq[j,] 
    tmp3 <- coocur_cal(tmp2, d, q)
    for(i in 1:f){
      co_tmp_list[[i]][j,] <- tmp3[i,]/co_entry_mean[i,] - sum(tmp3[i,])/co_rowsum_mean[i] - rowSums(tmp3)/co_rowsum_mean
    }
  }
  PMI_entry_sd_co <- rep(0,f)
  for(i in 1:f){
    PMI_entry_sd_co[i] <- sd(co_tmp_list[[g]][,i])/sqrt(n)
  }
  PMI_entry_sd_lr_co <- rep(0,f)
  tmp_g <- apply(co_tmp_list[[g]], MARGIN = 2, function(x) x - mean(x))
  for(i in 1:f){
    tmp_i <- apply(co_tmp_list[[i]], MARGIN = 2, function(x) x - mean(x))
    tmp_P_i <- P_SP[i,]
    PMI_entry_sd_lr_co[i] <- (tmp_P_g %*% (t(tmp_i) %*% tmp_i) %*% as.matrix(tmp_P_g) + tmp_P_i %*% (t(tmp_g) %*% tmp_g) %*% as.matrix(tmp_P_i) + 2* tmp_P_g %*% (t(tmp_i) %*% tmp_g) %*% as.matrix(tmp_P_i))/(n-1)
  }
  PMI_entry_sd_lr_co <- sqrt(PMI_entry_sd_lr_co/n)
  rm(tmp_i)
  rm(tmp_g)
  
  

  
  #####################    simulation on theorem 3.6 & 3.7   #######################
  ########## normality validation of low rank estimator ##############
  m = 100   #repeat m times, get m PMI matrices based on n EHRs
  PMI_entry_list_36 = matrix(0,m,f)
  PMI_entry_list_37 = matrix(0,m,f)
  cover_prob <- rep(0,f)
  cover_prob_lr <- rep(0,f)
  cover_prob_lr_2 <- rep(0,f)
  cover_prob_3 <- rep(0,f)
  cover_prob_lr_3 <- rep(0,f)
  cover_prob_4 <- rep(0,f)
  cover_prob_lr_4 <- rep(0,f)
  for(k in 1:m){
    total_seq <- matrix(0, nrow = n, ncol = c_len)
    for(l in 1:n){                                                   
      total_seq[l,] <- data_gene(c_len,p,alpha,V)
    }
    co_etr <- matrix(0,d,d)
    if(n == 1){
      tmp2 <- total_seq
      tmp3 <- coocur_cal(tmp2, d, q)
      co_etr <-  tmp3
    }else{
      for (j in 1:n) {
        tmp2 <- total_seq[j,] 
        tmp3 <- coocur_cal(tmp2, d, q)
        co_etr <-  co_etr + tmp3
      }
    }
    rm(tmp2)
    rm(tmp3)
    rm(total_seq)
    if(0 %in% rowSums(co_etr)){
      stop("Invalid co-occurrence matrix!")
    }
    SP_est0 <- SPPMI_calc(co = co_etr, q, d, c_len, n)                     #SP_est0 is the merged SPPMI matrix
    PMI_entry_list_36[k,1:f] <- SP_est0[g,1:f]#-PMI_mean_est[g,1:f]         #/PMI_entry_sd_co
    SP_est_lr <- SPPMI_calc_lr(co_etr, q, d, c_len, n,p)                     #SP_est_lr is the low rank estimator of SPPMI matrix
    PMI_entry_list_37[k,1:f] <- SP_est_lr[g, 1:f]#-PMI_lr_mean_est[g,1:f]   #/PMI_entry_sd_lr_co
    for(j in 1:f){
      if (abs(PMI_entry_list_36[k,j] - pmi_lr_true[g,j]) <= PMI_entry_sd_co[j]*qnorm(0.975)){
        cover_prob[j] <- cover_prob[j]+1
      }
      if (abs(PMI_entry_list_37[k,j] - pmi_lr_true[g,j]) <= PMI_entry_sd_lr_co[j]*qnorm(0.975)){
        cover_prob_lr[j] <- cover_prob_lr[j]+1
      }
      if (abs(PMI_entry_list_36[k,j] - pmi_true[g,j]) <= PMI_entry_sd_co[j]*qnorm(0.975)){
        cover_prob_3[j] <- cover_prob_3[j]+1
      }
      if (abs(PMI_entry_list_37[k,j] - pmi_true[g,j]) <= PMI_entry_sd_lr_co[j]*qnorm(0.975)){
        cover_prob_lr_3[j] <- cover_prob_lr_3[j]+1
      }
      if (abs(PMI_entry_list_37[k,j] - pmi_lr_true[g,j]) <= PMI_entry_sd_co[j]*qnorm(0.975)){
        cover_prob_4[j] <- cover_prob_4[j]+1   #no lr var with lr est
      }
      if (abs(PMI_entry_list_37[k,j] - pmi_lr_true[g,j]) <= PMI_entry_sd_lr_co[j]*qnorm(0.975)){
        cover_prob_lr_4[j] <- cover_prob_lr_4[j]+1
      }
    }
    rm(SP_est0)
    rm(SP_est_lr)
  }
  PMI_lr_mean_est <- PMI_lr_mean_est[g,1:f]
  PMI_mean_est <- PMI_mean_est[g,1:f]
  cover_prob <- cover_prob/m
  cover_prob_lr <- cover_prob_lr/m
  cover_prob_3 <- cover_prob_3/m
  cover_prob_lr_3 <- cover_prob_lr_3/m
  cover_prob_4 <- cover_prob_4/m
  cover_prob_lr_4 <- cover_prob_lr_4/m
  
  
  ################save results###############
  save(V, VVT, PMI_mean_est, PMI_entry_sd_co, PMI_entry_sd_lr_co, PMI_entry_list_36, PMI_entry_list_37,cover_prob, cover_prob_lr, cover_prob_3, cover_prob_lr_3, cover_prob_4, cover_prob_lr_4, pmi_true, pmi_lr_true, file=paste0('/n/data1/hsph/biostat/celehs/lab/zig728/PMI_inference/2407/1102CheckNormal_d',d,'_n',n,'_clen',c_len,'_ord',abs(ord),'.Rdata'))
}
