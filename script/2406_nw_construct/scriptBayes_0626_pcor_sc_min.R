setwd("/home/users/y2564li/graph_proj")

library(readr)
library(Matrix)
library(dplyr)
library(tidyr)

## ################################## ##
# read and merge data ####
check_merge_index <- function(data) {
  n <- length(data)
  idx <- sample(1:n)[1:5]
  for (i in idx) { if ( as.numeric(names(data)[i]) != i) stop("index does not match") }
}
## each RDS is a list of length 2: "pcor_sc_1200tp" "pcor_sc_200tp"
## each element in RDS is a list of 100/200/.../957 arrays of dimension 68*68*68
merge_pcor_sc_data_per_sess <- function(sess) {
  data_merged_1200tp <- list()
  data_merged_200tp <- list()
  labels <- c("1to100","101to200","201to300","301to400","401to500","501to600",
              "601to700","701to800","801to900","901to957")
  for (i in 1:length(labels)) {
    lab <- labels[i]
    data <- readRDS(paste0("data/corr/pcor_sc/fun_pcor_sc_cor_sess",sess,"_",lab,".RDS"))
    i_start <- (i-1)*100+1
    if ( i_start != as.numeric(sub("to.*","",lab)) ) stop("start index error")
    i_end <- -1
    if (i==length(labels)) { i_end <- 957 } 
    else { i_end <- i_start+99 }
    if ( i_end != as.numeric(sub(".*to","",lab)) ) stop("end index error")
    data[["pcor_sc_1200tp"]] <- data[["pcor_sc_1200tp"]][i_start:i_end] # data is a list of length i_end
    data[["pcor_sc_200tp"]] <- data[["pcor_sc_200tp"]][i_start:i_end] 
    names(data[["pcor_sc_1200tp"]]) <- as.character(i_start:i_end)
    data_merged_1200tp <- append(data_merged_1200tp, data[["pcor_sc_1200tp"]])
    names(data[["pcor_sc_200tp"]]) <- as.character(i_start:i_end)
    data_merged_200tp <- append(data_merged_200tp, data[["pcor_sc_200tp"]])
    check_merge_index(data_merged_1200tp)
    check_merge_index(data_merged_200tp)
  }
  list(pcor_sc_1200tp=data_merged_1200tp, pcor_sc_200tp=data_merged_200tp)
}

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
        # if either i or j is 1 or 35, leave min value and the region being conditioned on as NA
        if ( intersect(c(i,j),c(1,35))%>%length == 0 ) { 
          val <- arr[i,j,]
          val <- ifelse(is.na(val),Inf,abs(val)) # I may had some coding inconsistency between NA and Inf
          arr_min[i,j] <- min(val)
          arr_min[j,i] <- arr_min[i,j]
          arr_k[i,j] <- which.min(val)
          arr_k[j,i] <- arr_k[i,j]
        }
      }
    }
    res_min[[subj]] <- arr_min
    res_k[[subj]] <- arr_k
    if (subj %% 100 == 0) print(paste0("Finding min pcor: subject",subj," done!"))
  }
  list(min=res_min,k=res_k)
}


## ################################## ##
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
    res <- ifelse(is.na(res),0,res) # region pairs involving 1 and 35 are set as NA and cannot be compared to numbers
    diag(res) <- 0
    res
  } else if (type=="sig") {
    p_mat <- compute_ft_pval(cor_m=cor_mat,n=n_obs,d=d)
    res <- ifelse(p_mat<alpha, 1, 0)
    res <- ifelse(is.na(res),0,res)
    diag(res) <- 0
    res
  } else { stop("wrong type") }
}

# input: 
## pcor_min: list of 957 minimum correlation matrices of dimension 68*68
## n_obs: number of observations, n_cond: number of nodes being conditional on, 
## alpha: significance level
compute_pcor_low_order_nw_per_sess <- function(pcor_min,n_obs,n_cond,alpha) {
  res_nw_mag <- list()
  res_nw_sig <- list()
  res_nw_mag <- lapply(1:957,function(subj) { corr_to_nw(pcor_min[[subj]],type="mag") })
  res_nw_sig <- lapply(1:957,function(subj) { 
    corr_to_nw(pcor_min[[subj]],type="sig",n_obs=n_obs,d=n_cond,alpha=alpha) })
  list(mag=res_nw_mag, sig=res_nw_sig)
} # [[type]][[subj]]


# computation ####
merge_pcor_sc_data_saveRDS <- function() {
  res_1200tp <- list()
  res_200tp <- list()
  for (sess in 1:4) {
    data <- merge_pcor_sc_data_per_sess(sess=sess)
    res_1200tp <- append(res_1200tp, data[["pcor_sc_1200tp"]])
    res_200tp <- append(res_200tp, data[["pcor_sc_200tp"]])
    saveRDS(data[["pcor_sc_1200tp"]], paste0("data/corr/pcor_sc/fun_pcor_sc_1200tp_cor_sess",sess,".RDS"))
    saveRDS(data[["pcor_sc_200tp"]], paste0("data/corr/pcor_sc/fun_pcor_sc_200tp_cor_sess",sess,".RDS"))
  }
}

merge_pcor_sc_data_saveRDS()
print("########### Merge data and saveRDS done! ###########")

compute_nw_and_saveRDS <- function() {
  for (tp in c("1200tp","200tp")) {
    n_obs <- sub("*tp", "", tp)%>%as.numeric()
    for (sess in 1:4) {
      pcor_res_sess <- readRDS(paste0("data/corr/pcor_sc/fun_pcor_sc_",tp,"_cor_sess",sess,".RDS"))
      res <- compute_pcor_low_order_min_per_sess(pcor_res_sess=pcor_res_sess)
      res_min <- res[["min"]] # list of 957
      res_k <- res[["k"]] # list of 957
      saveRDS(res_min, paste0("data/corr/pcor_sc/fun_pcor_sc_",tp,"_min_sess",sess,".rds"))
      saveRDS(res_k, paste0("data/corr/pcor_sc/fun_pcor_sc_",tp,"_k_sess",sess,".rds"))
      res_nw <- compute_pcor_low_order_nw_per_sess(pcor_min=res_min,n_obs=n_obs,n_cond=3,alpha=.05/2278)
      saveRDS(res_nw, paste0("data/corr/pcor_sc/fun_pcor_sc_",tp,"_nw_sess",sess,".rds")) 
      print(paste0("########### Compute network for ",tp," session",sess," done! ###########"))
    }
  }
}
compute_nw_and_saveRDS()

####################
s1 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_1200tp_nw_sess1.rds")
s2 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_1200tp_nw_sess2.rds")
s3 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_1200tp_nw_sess3.rds")
s4 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_1200tp_nw_sess4.rds")
saveRDS(list(s1,s2,s3,s4),"data/corr/pcor_sc/fun_pcor_sc_1200tp_nw.rds")
saveRDS(list(s1[["mag"]],s2[["mag"]],s3[["mag"]],s4[["mag"]]),
        "data/corr/pcor_sc/fun_pcor_sc_1200tp_nw_mag.rds")
saveRDS(list(s1[["sig"]],s2[["sig"]],s3[["sig"]],s4[["sig"]]),
        "data/corr/pcor_sc/fun_pcor_sc_1200tp_nw_sig.rds")


s1 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_200tp_nw_sess1.rds")
s2 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_200tp_nw_sess2.rds")
s3 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_200tp_nw_sess3.rds")
s4 <- readRDS("data/corr/pcor_sc/fun_pcor_sc_200tp_nw_sess4.rds")
saveRDS(list(s1,s2,s3,s4),"data/corr/pcor_sc/fun_pcor_sc_200tp_nw.rds")
saveRDS(list(s1[["mag"]],s2[["mag"]],s3[["mag"]],s4[["mag"]]),
        "data/corr/pcor_sc/fun_pcor_sc_200tp_nw_mag.rds")
saveRDS(list(s1[["sig"]],s2[["sig"]],s3[["sig"]],s4[["sig"]]),
        "data/corr/pcor_sc/fun_pcor_sc_200tp_nw_sig.rds")



