setwd("/home/users/y2564li/graph_proj")

library(readr)
library(Matrix)
library(dplyr)
library(tidyr)

## ################################## ##
# read data ####
pcor_1200tp_sess1 <- readRDS("data/corr/fun_pcor_ord1_1200tp_cor_sess1.RDS")
pcor_1200tp_sess2 <- readRDS("data/corr/fun_pcor_ord1_1200tp_cor_sess2.RDS")
pcor_1200tp_sess3 <- readRDS("data/corr/fun_pcor_ord1_1200tp_cor_sess3.RDS")
pcor_1200tp_sess4 <- readRDS("data/corr/fun_pcor_ord1_1200tp_cor_sess4.RDS")

## ################################## ##
# get min and cond k ####
# input: a list of 957 elements, each element containing a 68*68*68 array
compute_pcor_low_order_min_per_sess <- function(pcor_res_sess) {
  res_min <- list()
  res_k <- list()
  for (subj in 1:957) {
    arr <- pcor_res_sess[[subj]]
    arr_min <- matrix(ncol=68,nrow=68)
    arr_k <- matrix(ncol=68,nrow=68)
    for (i in 1:67) {
      for (j in (i+1):68) {
        val <- arr[i,j,]
        val <- ifelse(is.na(val),Inf,abs(val))
        arr_min[i,j] <- min(val)
        arr_min[j,i] <- arr_min[i,j]
        arr_k[i,j] <- which.min(val)
        arr_k[j,i] <- arr_k[i,j]
      }
    }
    res_min[[subj]] <- arr_min
    res_k[[subj]] <- arr_k
    print(paste0("subject",subj," done!"))
  }
  list(min=res_min,k=res_k)
}


## ################################## ##
# network ####
compute_ft_pval <- function(cor_m,n,d=0){
  cor_ft <- atanh(cor_m)
  se <- 1/sqrt(n-d-3)
  z_val <- cor_ft/se
  cor_ft_p <- 2 * pnorm(-abs(z_val))+diag(rep(1, 68))
  cor_ft_p
}

val_to_nw <- function(mat, n_edges=136,select_small=F,p=68) {
  select_large <- (!select_small)
  thre <- sort(mat[upper.tri(diag(1:p))],decreasing=select_large)[n_edges]
  if (select_small==T) {
    res <- ifelse(mat<=thre,1,0)
    return(res)
  }
  res <- ifelse(mat>=thre,1,0)
  return(res)
}

corr_to_nw <- function(cor_mat,type,n_edges=136,n_obs=NULL,d=NULL,alpha=NULL) {
  if (type=="mag") {
    res <- val_to_nw(abs(cor_mat), n_edges=n_edges)
    diag(res) <- 0
    res
  } else if (type=="sig") {
    p_mat <- compute_ft_pval(cor_m=cor_mat,n=n_obs,d=d)
    res <- ifelse(p_mat<alpha, 1, 0)
    diag(res) <- 0
    res
  } else { print("wrong type") }
}

compute_pcor_low_order_nw_per_sess <- function(pcor_min,n_obs,n_cond,alpha) {
  res_nw_mag <- list()
  res_nw_sig <- list()
  res_nw_mag <- lapply(1:957,function(subj) { corr_to_nw(pcor_min[[subj]],type="mag") })
  res_nw_sig <- lapply(1:957,function(subj) { 
    corr_to_nw(pcor_min[[subj]],type="sig",n_obs=n_obs,d=n_cond,alpha=alpha) })
  list(mag=res_nw_mag, sig=res_nw_sig)
} # [[type]][[sess_tp]][[subj]]


## ################################## ##
# compute for each sess ####
# sess1
pcor_1200tp_res <- compute_pcor_low_order_min_per_sess(pcor_1200tp_sess1)
pcor_1200tp_min <- pcor_1200tp_res[["min"]]
pcor_1200tp_k <- pcor_1200tp_res[["k"]]
saveRDS(pcor_1200tp_min, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_min_sess1.rds")
saveRDS(pcor_1200tp_k, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_k_sess1.rds")
pcor_1200tp_nw <- compute_pcor_low_order_nw_per_sess(pcor_1200tp_min,n_obs=1200,n_cond=1,alpha=.05/2278)
saveRDS(pcor_1200tp_nw, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess1.rds") 
print("#################### session 1 done! ####################")
# sess2
pcor_1200tp_res <- compute_pcor_low_order_min_per_sess(pcor_1200tp_sess2)
pcor_1200tp_min <- pcor_1200tp_res[["min"]]
pcor_1200tp_k <- pcor_1200tp_res[["k"]]
saveRDS(pcor_1200tp_min, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_min_sess2.rds")
saveRDS(pcor_1200tp_k, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_k_sess2.rds")
pcor_1200tp_nw <- compute_pcor_low_order_nw_per_sess(pcor_1200tp_min,n_obs=1200,n_cond=1,alpha=.05/2278)
saveRDS(pcor_1200tp_nw, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess2.rds") 
print("#################### session 2 done! ####################")

# sess3
pcor_1200tp_res <- compute_pcor_low_order_min_per_sess(pcor_1200tp_sess3)
pcor_1200tp_min <- pcor_1200tp_res[["min"]]
pcor_1200tp_k <- pcor_1200tp_res[["k"]]
saveRDS(pcor_1200tp_min, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_min_sess3.rds")
saveRDS(pcor_1200tp_k, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_k_sess3.rds")
pcor_1200tp_nw <- compute_pcor_low_order_nw_per_sess(pcor_1200tp_min,n_obs=1200,n_cond=1,alpha=.05/2278)
saveRDS(pcor_1200tp_nw, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess3.rds") 
print("#################### session 3 done! ####################")

# sess4
pcor_1200tp_res <- compute_pcor_low_order_min_per_sess(pcor_1200tp_sess4)
pcor_1200tp_min <- pcor_1200tp_res[["min"]]
pcor_1200tp_k <- pcor_1200tp_res[["k"]]
saveRDS(pcor_1200tp_min, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_min_sess4.rds")
saveRDS(pcor_1200tp_k, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_k_sess4.rds")
pcor_1200tp_nw <- compute_pcor_low_order_nw_per_sess(pcor_1200tp_min,n_obs=1200,n_cond=1,alpha=.05/2278)
saveRDS(pcor_1200tp_nw, "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess4.rds") 
print("#################### session 4 done! ####################")

# combined list
pcor_1200tp_nw_sess1 <- readRDS("data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess1.rds") 
pcor_1200tp_nw_sess2 <- readRDS("data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess2.rds") 
pcor_1200tp_nw_sess3 <- readRDS("data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess3.rds") 
pcor_1200tp_nw_sess4 <- readRDS("data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw_sess4.rds") 
saveRDS(list(pcor_1200tp_nw_sess1, pcor_1200tp_nw_sess2, pcor_1200tp_nw_sess3, pcor_1200tp_nw_sess4),
        "data/corr/pcor_ord1/fun_pcor_ord1_1200tp_nw.rds")











