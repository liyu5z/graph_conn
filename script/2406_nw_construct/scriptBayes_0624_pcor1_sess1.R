setwd("/home/users/y2564li/graph_proj")

library(readr)
library(Matrix)
library(dplyr)
library(tidyr)
library(glasso)

conn_list = readRDS("/home/users/y2564li/graph_proj/data/fun_str_common.rds") 
fun_ts1_common = conn_list[[1]]
fun_ts2_common = conn_list[[2]]
fun_ts3_common = conn_list[[3]]
fun_ts4_common = conn_list[[4]] 
str_conn_common = conn_list[[5]]

ts = list(fun_ts1_common,fun_ts2_common,fun_ts3_common,fun_ts4_common) ## ts[[subj]][[sess]] subj=1,...,957, sess=1,2,3,4

ts_sub6 <- list()
for (sess in 1:4) {
  ts_sub6[[sess]] <- list()
  for (subj in 1:957) { ts_sub6[[sess]][[subj]] <- ts[[sess]][[subj]][seq(1,1200,by=6),] }
}

compute_pcor_ord1 <- function(i,j,k,data) {
  cor_data <- cor(data)
  num <- cor_data[i,j] - cor_data[i,k]*cor_data[j,k]
  den <- sqrt( (1-(cor_data[i,k])^2) * (1-(cor_data[j,k])^2) )
  num/den
}

count_edges <- function(conn,p=68) {
  conn <- ifelse(conn==0,0,1)
  return(sum(conn[upper.tri(diag(1:p))]))
}

sess <- 1
pcor_order1 <- list()
pcor_order1_sub6 <- list()
# initialize: pcor_order1[i,j,k]=cor(i,j|k)
for (subj in 1:957) {
  pcor_order1[[subj]] <- array(rep(NA,68^3),dim=c(68,68,68)) 
  pcor_order1_sub6[[subj]] <- array(rep(NA,68^3),dim=c(68,68,68)) 
}
t_start = Sys.time()
for (subj in 1:957) {
  for (i in 1:68) {
    pcor_order1[[subj]][i,,] <- t(sapply(X=1:68,FUN=compute_pcor_ord1,i=i,k=1:68,data=ts[[sess]][[subj]]))
    pcor_order1_sub6[[subj]][i,,] <- t(sapply(X=1:68,FUN=compute_pcor_ord1,i=i,k=1:68,data=ts_sub6[[sess]][[subj]]))
  } # end i
  print(paste0("subject",subj," done!"))
} # end subj
t_end = Sys.time()
print(t_end-t_start)

saveRDS(pcor_order1, "data/corr/fun_pcor_ord1_1200tp_cor_sess1.RDS")
saveRDS(pcor_order1_sub6, "data/corr/fun_pcor_ord1_200tp_cor_sess1.RDS")

