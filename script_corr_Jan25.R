# load libraries
library(readr)
library(ggplot2)
library(glmnet)
library(hdi)
library(magrittr)
library(Matrix)
library(ppcor)
library(tidyr)

## read structural connection data ####
path_data = '~/git_data/GraphMixed/Data/100307_aparc_connectome/'
id_subj = "100307"
communities_tb = read.table('~/git_data/GraphMixed/Data/100307_aparc_connectome/communities_aparc.csv', header=T, sep=',')

idx_ctx = c(1:3,5:35,36:38,40:70)
str_conn_raw = read.table(paste0(path_data, id_subj, '_structural', '.csv'), header=F, sep=' ')
str_conn = as.matrix(sparseMatrix(str_conn_raw$V1, str_conn_raw$V2, x = str_conn_raw$V3, symmetric = TRUE))[idx_ctx,idx_ctx]

(sum(str_conn==0)-68)/2 # =320, pairs of regions that are not connected

# connection represented by a binary matrix, 1=connected, 0=otherwise
str_conn.m <- ifelse(str_conn>0,1,0) + diag(rep(1,68)) 

# ***************************************** #
## read function data ####
ts = list()   # List of time-series for the 4 available MRI scan sessions
for (id_sess in 1:4) {
  ts[[id_sess]] = read.table(paste0(path_data, id_subj, '_', id_sess, '.csv'), header=F, sep=',')
}
  ## ts: list of 4 sessions
  ## ts[[i]]: list of 68 regions in ith time session with 1200 time points
# subsampled data
ts.sub6 <- list()
ts.sub6[[1]] <-ts[[1]][seq(1,1200,by=6),] # 1200=nrow(each session)
ts.sub6[[2]] <-ts[[2]][seq(1,1200,by=6),]
ts.sub6[[3]] <-ts[[3]][seq(1,1200,by=6),]
ts.sub6[[4]] <-ts[[4]][seq(1,1200,by=6),]

# ***************************************** #
# 1. Covariance ####
## full sample 
cov.avg <- ( cov(ts[[1]])+cov(ts[[2]])+cov(ts[[3]])+cov(ts[[4]]) )/4
image(cov.avg[communities_tb$id_ROI,communities_tb$id_ROI])

## subsample 
cov.sub6.avg <- ( cov(ts.sub6[[1]])+cov(ts.sub6[[2]])+cov(ts.sub6[[3]])+cov(ts.sub6[[4]]) )/4
image(cov.sub6.avg[communities_tb$id_ROI,communities_tb$id_ROI])

# ***************************************** #
# 2. Full partial correlation ####
## full sample
pcor.avg <- ( pcor(ts[[1]])$estimate+pcor(ts[[2]])$estimate+pcor(ts[[3]])$estimate+pcor(ts[[4]])$estimate )/4
image( pcor.avg ) # plot
image( pcor.avg[communities_tb$id_ROI,communities_tb$id_ROI] ) # plot by communities

## subsample
pcor.sub6.avg <- (pcor(ts.sub6[[1]])$estimate+pcor(ts.sub6[[2]])$estimate+
                    pcor(ts.sub6[[3]])$estimate+pcor(ts.sub6[[4]])$estimate)/4
image( pcor.sub6.avg ) # plot
image( pcor.sub6.avg[communities_tb$id_ROI,communities_tb$id_ROI] ) # plot by communities


## Fisher transformation ####
pnorm(abs(atanh( pcor.avg ))/ (1/sqrt(1200-3)), lower.tail=F)  # p-values using entire sample 
pnorm(abs(atanh( pcor.sub6.avg ))/ (1/sqrt(200-3)), lower.tail=F)  # p-values using subsample

# ***************************************** #
# 3. Low-order conditional dependence ####
## function to compute partial correaltions conditional on one variable
cond_cor_1 <- function(i,j,k,data) {
  num <- cor(data)[i,j] - cor(data)[i,k]*cor(data)[j,k]
  den <- sqrt( (1-(cor(data)[i,k])^2) * (1-(cor(data)[j,k])^2) )
  num/den
}
## full sample
  # initialize arrays to store correlations, (i,j,k)=rho_ij|k
pcor.order1 <- list()
pcor.order1[[1]] <- array(rep(NA,68^3),dim=c(68,68,68))
pcor.order1[[2]] <- array(rep(NA,68^3),dim=c(68,68,68))
pcor.order1[[3]] <- array(rep(NA,68^3),dim=c(68,68,68))
pcor.order1[[4]] <- array(rep(NA,68^3),dim=c(68,68,68))
 # compute correlations
for (i in 1:68) {
  pcor.order1[[1]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts[[1]]))
  pcor.order1[[2]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts[[2]]))
  pcor.order1[[3]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts[[3]]))
  pcor.order1[[4]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts[[4]]))
  print(paste0(i," done"))
}
  # average over 4 sessions
pcor.order1.avg <- (pcor.order1[[1]]+pcor.order1[[2]]+pcor.order1[[3]]+pcor.order1[[4]])/4
  # find the minimal correlation
pcor.order1.avg.min <- list()
pcor.order1.avg.min$m <- matrix(nrow=68,ncol=68)
pcor.order1.avg.min$k <- matrix(nrow=68,ncol=68)
for (i in 1:67) {
  for (j in (i+1):68) {
    val <- pcor.order1.avg[i,j,]
    val[is.na(val)] <- 100
    k <- sort(abs(val),index.return=T)$ix[1]
    pcor.order1.avg.min$k[i,j] <- k
    pcor.order1.avg.min$k[j,i] <- k
    pcor.order1.avg.min$m[i,j] <- val[k]
    pcor.order1.avg.min$m[j,i] <- val[k]
  }
}

## subsample
  # initialize arrays to store correlations, (i,j,k)=rho_ij|k
pcor.sub6.order1 <- list()
pcor.sub6.order1[[1]] <- array(rep(NA,68^3),dim=c(68,68,68))
pcor.sub6.order1[[2]] <- array(rep(NA,68^3),dim=c(68,68,68))
pcor.sub6.order1[[3]] <- array(rep(NA,68^3),dim=c(68,68,68))
pcor.sub6.order1[[4]] <- array(rep(NA,68^3),dim=c(68,68,68))
  # compute correlations
for (i in 1:68) {
  pcor.sub6.order1[[1]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts.sub6[[1]]))
  pcor.sub6.order1[[2]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts.sub6[[2]]))
  pcor.sub6.order1[[3]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts.sub6[[3]]))
  pcor.sub6.order1[[4]][i,,] <- t(sapply(X=1:68,FUN=cond_cor_1,i=i,k=1:68,data=ts.sub6[[4]]))
  print(paste0(i," done"))
}
  # average over 4 sessions
pcor.sub6.order1.avg <- (pcor.sub6.order1[[1]]+pcor.sub6.order1[[2]]+
                           pcor.sub6.order1[[3]]+pcor.sub6.order1[[4]])/4
  # find minimal correlations
pcor.sub6.order1.avg.min <- list()
pcor.sub6.order1.avg.min$m <- matrix(nrow=68,ncol=68)
pcor.sub6.order1.avg.min$k <- matrix(nrow=68,ncol=68)
for (i in 1:67) {
  for (j in (i+1):68) {
    val <- pcor.sub6.order1.avg[i,j,]
    val[is.na(val)] <- 100
    k <- sort(abs(val),index.return=T)$ix[1]
    pcor.sub6.order1.avg.min$k[i,j] <- k
    pcor.sub6.order1.avg.min$k[j,i] <- k
    pcor.sub6.order1.avg.min$m[i,j] <- val[k]
    pcor.sub6.order1.avg.min$m[j,i] <- val[k]
  }
}


# ***************************************** #
## 4. Summaries ####
### 4.1 Find pairs of regions with largest correlations ####
# Helper function: convert an index in a flattened/vectorized matrix to (row, col) index
## input: i=the index in a vectorized matrix, p=nrow (or ncol) of the squared matrix,
##    sym=F if the matrix is asymmetric and the function returns [row,col], 
##       =T if symmetric and the function returns row&col in an increasing order (note [row,col] and [col,row] equivalent in this case)
## output: index of row, index of column
get_sqrmat_index <- function(i,p,sym=FALSE) {
  t(sapply(X=i, FUN=function(ind) {
    i.r <- row(diag(1:p))[ind]
    i.c <- col(diag(1:p))[ind]
    if (sym==T) { c(min(i.r,i.c), max(i.r,i.c) ) } 
    else {c(i.r, i.c)}
  }))
}
# Helper function: find n values with largest magnitudes in a matrix
##  input: a correlation/symmetric matrix, and 
##    specify the number of greatest values (in absolute values) to find
##  output: n values of greatest magnitude
get_large_cor <- function(m.cor,n) {
  sort(abs(m.cor),decreasing=T)[(1:n)*2]
}
# Helper function: find the row and col indexes in a matrix corresponding to the n largest values (in magnitude)
##  input: mat=a correlation matrix, n=the number of greatest values (in absolute values) to find
##    p=nrow or ncol of the matrix where we are going to get (row, col) indexes
##    sym for the use in get_sqrmat_index, tol=a small value to allow errors in computation precisions
##  output: matrix indexes in the form of (row,col) which represent the regions of top n values
get_top_conn <- function(mat,n,p=68,sym=T,tol=1e-6) {
  if (sum(is.na(diag(mat)))==0) { mat <- mat-diag(diag(mat)) } # remove diagonal elements which are cor/cov with itself
  val.top <- get_large_cor(mat,n) # largest correlation values
  val.index <- sapply(val.top, FUN=function(val) {which(abs(abs(mat)-val)<=tol)})[1,] # get indexes
  get_sqrmat_index(val.index,p,sym)
}

## full sample - 20 pairs of regions with strongest correlations
cov.top20.conn <- get_top_conn(cov.avg, n=20)
cor.top20.conn <- get_top_conn(cor.avg, n=20)
pcor.top20.conn <- get_top_conn(pcor.avg,n=20)
pcor.order1.top20.conn <- get_top_conn(pcor.order1.avg.min$m,n=20)
## subsample - 20 pairs of regions with strongest correlations
cov.sub6.top20.conn <- get_top_conn(cov.sub6.avg,n=20)
cor.sub6.top20.conn <- get_top_conn(cor.sub6.avg,n=20)
pcor.sub6.top20.conn <- get_top_conn(pcor.sub6.avg,n=20)
pcor.sub6.order1.top20.conn <- get_top_conn(pcor.sub6.order1.avg.min$m,n=20)


### 4.2 Find common pairs ####

# Helper function: count the number of same pairs
##  input: two n pairs of regions
##  output: number of identical pairs
count_common_pairs <- function(conn1,conn2) {
  sum(apply(conn1,MARGIN=1,function(row1){
    apply(conn2,MARGIN=1,function(row2) {identical(row1,row2)}) }))
}
  # common pairs - entire sample vs subsample
count_common_pairs(cov.top20.conn,cov.sub6.top20.conn) # 19 pairs
count_common_pairs(cor.top20.conn,cor.sub6.top20.conn) # 16 pairs
count_common_pairs(pcor.top20.conn,pcor.sub6.top20.conn) # 17 pairs
count_common_pairs(pcor.order1.top20.conn,pcor.sub6.order1.top20.conn) # 19 pairs
# may also compare differnt measures based on the same dataset, e.g.
count_common_pairs(cov.top20.conn,cor.top20.conn)

### 4.3 Region community ####
# Helper function: determine whether pairs of regions are from the same community
##  input: n pairs of regions
##  output: n rows, each containing a T/F indicating whether they are in the same community and the two communities
if_same_comm <- function(conn) {
  t(apply(conn,MARGIN=1,FUN=function(x) {
    cbind(identical(communities_tb[communities_tb$id_ROI==x[1],"name_community"],
                    communities_tb[communities_tb$id_ROI==x[2],"name_community"]) ,
          communities_tb[communities_tb$id_ROI==x[1],"name_community"],
          communities_tb[communities_tb$id_ROI==x[2],"name_community"] ) }))
}
# e.g.
if_same_comm(cov.sub6.top20.conn) # only 1 pair contains regions from different communities





