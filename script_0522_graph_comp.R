library(igraph)
library(graphkernels)
library(ggplot2)

compute_shortest_path_matrix <- function(mat,use_NA=FALSE) {
  G <- igraph::graph_from_adjacency_matrix(mat, mode='undirected', weighted = NULL)
  res <- igraph::distances(G, algorithm="dijkstra")
  if (use_NA) { res <- ifelse(res==Inf,NA,res)-diag(NA,68) }
  else { res <- ifelse(res==Inf,0,res) }
  return(res)
}
compute_node_degrees <- function(nw) {
  return( sapply(1:68,function(i){ sum(nw[i,])-1 }) )
}
compute_fp_bin <- function(m_str,m_fun) {
  res <- count_edges(ifelse(m_str==1,0,m_fun)>0)
  return(res)
}
compute_fp_list_bin <- function(fun_nw,str_nw,
                               fnames=c("mag_all","sig_all","mag_sub6","sig_sub6"),
                               snames=c("sum_and_thre","thre_and_sum")) {
  df <- data.frame()
  for (fname in fnames) {
    for (sname in snames) {
      m_str <- ifelse(str_nw[[sname]]==0,0,1)
      m_fun <- ifelse(fun_nw[[fname]]==0,0,1)
      fp <- compute_fp(m_str,m_fun)
      if (wrsgraph::edit_dist_graph(m_fun,m_str)$fp != fp) print("different fp!")
      df <- rbind(df,c(fname,sname,fp))
    }
  }
  colnames(df) <- c("fun_nw","str_nw","false_positive")
  return(df)
}

compute_fp <- function(m_str,m_fun) { 
  res <- count_edges(ifelse(m_str>0,0,m_fun)) 
  return(res)
}
compute_ed_fp_list <- function(fun_nw,str_nw,
                               fnames=c("mag_all","sig_all","mag_sub6","sig_sub6"),
                               snames=c("sum_and_thre","thre_and_sum")) {
  df <- data.frame()
  for (fname in fnames) {
    for (sname in snames) {
      ed <- sum(str_nw[[sname]]!=fun_nw[[fname]])/2
      fp <- compute_fp(str_nw[[sname]],fun_nw[[fname]])
      if (wrsgraph::edit_dist_graph(fun_nw[[fname]],str_nw[[sname]])$ed != ed) print("different ed!")
      if (wrsgraph::edit_dist_graph(fun_nw[[fname]],str_nw[[sname]])$fp != fp) print("different fp!")
      df <- rbind(df,c(fname,sname,ed,fp))
    }
  }
  colnames(df) <- c("fun_nw","str_nw","edit_dist","false_positive")
  return(df)
}



for (comm in unique(communities_tb$id_community)) {
  id_roi <- subset(communities_tb,id_community==comm)$id_ROI
  print(c(comm,mean(str_gt_node_degrees[[1]][id_roi])))
}

# 1. Structural networks ####
communities_tb = read.table('Data_YL/communities_aparc.csv', header=T, sep=',')
## 1.1 Global thresholding ####
str_nw_gt <- readRDS("../rds_nw/str_nw_gt.rds")
str_gt_shortest_path <- lapply(str_nw_gt[1:2],compute_shortest_path_matrix)
rbind("bin_conn"=sapply(str_nw_gt[1:2],count_edges),
      "path_conn"=sapply(str_gt_shortest_path,function(nw) {count_edges(nw>0)}))

str_gt_node_degrees <- lapply(str_nw_gt[1:2],compute_node_degrees)

## 1.2 Local thresholding ####
str_nw_lt <- readRDS("../rds_nw/str_nw_lt.rds")
str_lt_shortest_path <- lapply(str_nw_lt[1:2],compute_shortest_path_matrix)
rbind("bin_conn"=sapply(str_nw_lt[1:2],count_edges),
      "path_conn"=sapply(str_lt_shortest_path,function(nw) {count_edges(nw>0)}))

str_lt_node_degrees <- lapply(str_nw_lt[1:2],compute_node_degrees)

## 1.3 Sign test ####
str_nw_st <- readRDS("../rds_nw/str_nw_st.rds")
str_st_shortest_path <- lapply(str_nw_st[c(1,3)],compute_shortest_path_matrix)
rbind("bin_conn"=sapply(str_nw_st[c(1,3)],count_edges),
      "path_conn"=sapply(str_st_shortest_path,function(nw) {count_edges(nw>0)}))

str_st_node_degrees <- lapply(str_nw_st[c(1,3)],compute_node_degrees)

## 1.4 summary-path ####
str_pl_summary <- cbind(cbind("bin_conn"=sapply(str_nw_gt[1:2],count_edges),
                              "path_conn"=sapply(str_gt_shortest_path,function(nw) {count_edges(nw>0)})),
                        cbind(sapply(str_nw_lt[1:2],count_edges),
                              sapply(str_lt_shortest_path,function(nw) {count_edges(nw>0)})),
                        cbind(sapply(str_nw_st[c(1,3)],count_edges),
                              sapply(str_st_shortest_path,function(nw) {count_edges(nw>0)})) )
str_pl_summary <- rbind(rep(c("bin_conn","path_conn"),3),str_pl_summary)
colnames(str_pl_summary) <- rep(c("global_thre","local_thre","sign_test"),each=2)
data.frame(str_pl_summary)

identical(str_gt_shortest_path[[1]]>0,str_gt_shortest_path[[2]]>0) # TRUE
### 2 path length networks given by gt have same places of connection, though path lengths differ
identical(str_lt_shortest_path[[1]]>0,str_lt_shortest_path[[2]]>0) # TRUE since both are guaranteed to be connected graphs

### barplot ####
par(mfrow=c(1,2))
barplot(table(str_gt_shortest_path[[1]][upper.tri(diag(68))]),col="lightgreen",
     main="global thre path len (sum_and_thre)",xlab="path length")
barplot(table(str_gt_shortest_path[[2]][upper.tri(diag(68))]),col="lightgreen",
     main="global thre path len (thre_and_sum)",xlab="path length")
par(mfrow=c(1,2))
barplot(table(str_lt_shortest_path[[1]][upper.tri(diag(68))]),col="lightgreen",
     main="local thre path len (sum_and_thre)",xlab="path length")
barplot(table(str_lt_shortest_path[[2]][upper.tri(diag(68))]),col="lightgreen",
     main="local thre path len (thre_and_sum)",xlab="path length")
par(mfrow=c(1,2))
barplot(table(str_st_shortest_path[[1]][upper.tri(diag(68))]),col="lightgreen",
     main="sign test path len (thre_count_p027)",xlab="path length")
barplot(table(str_st_shortest_path[[2]][upper.tri(diag(68))]),col="lightgreen",
     main="sign test path len (thre_count_p027)",xlab="path length")



# ################################################################ #
# ################################################################ #
# 2. Functional networks ####
## 2.1 Marginal correlation ####
fun_nw_mcor <- readRDS("../rds_nw/fun_nw_mcor.rds")
fun_mcor_shortest_path <- lapply(fun_nw_mcor[1:4],compute_shortest_path_matrix)
fun_mcor_node_degrees <- lapply(fun_nw_mcor[1:4],compute_node_degrees)
# mag_all: (62,28) prop==1, cutoff prop==0.3437827
# mag_sub6: (62,28) prop==1, cutoff prop==0.3322884
# sig_all: 34 pairs prop==1, cutoff prop=0.9955590; 62-57-28, no 62-28
# sig_sub6: 5 pairs prop==1, curoff prop==0.9195402

## 2.2 Partial correlation ####
fun_nw_pcor <- readRDS("../rds_nw/fun_nw_pcor.rds")
fun_pcor_shortest_path <- lapply(fun_nw_pcor[1:4],compute_shortest_path_matrix)
fun_pcor_node_degrees <- lapply(fun_nw_pcor[1:4],compute_node_degrees)

sort(fun_nw_pcor[["mag_sub6_prop"]][upper.tri(diag(68))],decreasing = T)[1:136]
arrayInd(which(fun_nw_pcor[["mag_sub6_prop"]]-diag(68)==1),c(68,68))
# mag_all: 2 pairs prop==1 ((58,24),(62,28)), cutoff prop==0.2813480
# mag_sub6: no prop==1, cutoff prop==0.1361024
# sig_all: 2 pairs prop==1 ((58,24),(62,28)), cutoff prop=0.2071578
# sig_sub6: no pairs prop==1, curoff prop==0.002873563


## 2.3 Low-order partial correlation ####
fun_nw_pcor_order1 <- readRDS("../rds_nw/fun_nw_pcor_order1.rds") 
fun_pcor_order1_shortest_path <- lapply(fun_nw_pcor_order1[1:4],compute_shortest_path_matrix)
fun_pcor_order1_node_degrees <- lapply(fun_nw_pcor_order1[1:4],compute_node_degrees)


## 2.4 Graphical lasso ####
fun_nw_gl <- readRDS("../rds_nw/fun_nw_gl.rds") 
fun_gl_shortest_path <- lapply(fun_nw_gl[1:2],compute_shortest_path_matrix)
fun_gl_node_degrees <- lapply(fun_nw_gl[1:2],compute_node_degrees)

## 2.5 summary-path ####
fun_pl_summary <- cbind(cbind(sapply(fun_nw_mcor[1:4],count_edges),
                              sapply(fun_mcor_shortest_path,function(nw) {count_edges(nw>0)})),
                        cbind(sapply(fun_nw_pcor[1:4],count_edges),
                              sapply(fun_pcor_shortest_path,function(nw) {count_edges(nw>0)})),
                        cbind(sapply(fun_nw_pcor_order1[1:4],count_edges),
                              sapply(fun_pcor_order1_shortest_path,function(nw) {count_edges(nw>0)})),
                        cbind(c(sapply(fun_nw_gl[1:2],count_edges),NA,NA),
                              c(sapply(fun_gl_shortest_path,function(nw) {count_edges(nw>0)}),NA,NA)))
fun_pl_summary <- rbind(rep(c("bin_conn","path_conn"),4),fun_pl_summary)
colnames(fun_pl_summary) <- rep(c("marg_cor","par_cor","pcor_order1","glasso"),each=2)
fun_pl_summary

## 2.5 plot-path ####
par(mfrow=c(2,2))
sapply(1:4,function(i) {
  nw <- fun_mcor_shortest_path[[i]]
  barplot(table(nw[upper.tri(diag(68))]),col="lightblue",xlab="path length",
          main=paste0("marginal corr path len (",names(fun_mcor_shortest_path)[i],")")) })
par(mfrow=c(2,2))
sapply(1:4,function(i) {
  nw <- fun_pcor_shortest_path[[i]]
  barplot(table(nw[upper.tri(diag(68))]),col="lightblue",xlab="path length",
          main=paste0("partial corr path len (",names(fun_pcor_shortest_path)[i],")")) })
par(mfrow=c(2,2))
sapply(1:4,function(i) {
  nw <- fun_pcor_order1_shortest_path[[i]]
  barplot(table(nw[upper.tri(diag(68))]),col="lightblue",xlab="path length",
          main=paste0("first-order pcorr path len (",names(fun_pcor_order1_shortest_path)[i],")")) })
par(mfrow=c(1,2))
sapply(1:2,function(i) {
  nw <- fun_gl_shortest_path[[i]]
  barplot(table(nw[upper.tri(diag(68))]),col="lightblue",xlab="path length",
          main=paste0("glasso path len (", names(fun_gl_shortest_path)[i],")")) })

## 2.6 prop summary ####
compute_prop_summary <- function(fun_nw,fun_method) {
  df <- data.frame()
  prop1 <- list()
  fnames <- names(fun_nw)
  fnames <- fnames[grepl("prop$", fnames)]
  for (fname in fnames) {
    cutoff <- sort(fun_nw[[fname]][upper.tri(diag(68))],decreasing=T)[136]
    df_cur <- c(fun_method, fname,cutoff)
    df <- rbind(df,df_cur)
    prop1_cur <- matrix(arrayInd(which(fun_nw[[fname]]-diag(68)==1),c(68,68)),ncol=2)
    prop1[[fname]] <- prop1_cur
  }
  colnames(df) <- c("fun_method","fun_nw","prop_cutoff")
  return(list(df=df,prop1=prop1))
}
compute_prop_summary(fun_nw_mcor,"mcor")$prop1
compute_prop_summary(fun_nw_pcor,"pcor")$prop1
compute_prop_summary(fun_nw_pcor_order1,"pcor1")$prop1
compute_prop_summary(fun_nw_gl,"gl")$prop1
df_fun_prop_cutoff <- rbind(compute_prop_summary(fun_nw_mcor,"mcor")$df,
                            compute_prop_summary(fun_nw_pcor,"pcor")$df,
                            compute_prop_summary(fun_nw_pcor_order1,"pcor1")$df,
                            compute_prop_summary(fun_nw_gl,"gl")$df)
df_fun_prop_cutoff$prop_cutoff <- as.numeric(df_fun_prop_cutoff$prop_cutoff)
df_fun_prop_cutoff$prop_cutoff <- round(df_fun_prop_cutoff$prop_cutoff,4)

# N_edges
c("N_edges",sapply(append(str_nw_gt[1:2],append(str_nw_lt[1:2],str_nw_st[c(1,3)])),count_edges))
sapply(append(fun_nw_mcor[1:4],append(fun_nw_pcor[1:4],fun_nw_pcor_order1[1:4])),count_edges)
sapply(fun_nw_gl[1:2],count_edges)

# ################################################################ #
# ################################################################ #
# 3. Comparison ####
## 3.0 plot FP changes ####
compare_conn <- function(m_fun,m_str,m_str_path) {
  df <- data.frame()
  df <- rbind(df, c(count_edges(ifelse(m_fun==1&m_str==1,1,0)),
                    count_edges(ifelse(m_fun==1&m_str_path>=1,1,0)),
                    count_edges(ifelse(m_fun==1&m_str==0,1,0)),
                    count_edges(ifelse(m_fun==1&m_str_path==0,1,0))))
  colnames(df) <- c("tp","tp_path","fp","fp_path")
  return(df)
}
compare_conn_list <- function(fun_nw,str_nw,str_path,fun_method,str_method,
                              fnames=c("mag_all","mag_sub6","sig_all","sig_sub6")) {
  snames <- names(str_path)
  df <- data.frame()
  for (fname in fnames) {
    for (sname in snames) {
      df_cur <- compare_conn(fun_nw[[fname]],str_nw[[sname]],str_path[[sname]])
      df_cur$fun_nw <- fname
      df_cur$str_nw <- sname
      df <- rbind(df,df_cur)
    }
  }
  df$fun_method <- fun_method
  df$str_method <- str_method
  df$tp <- as.numeric(df$tp)
  df$fp <- as.numeric(df$fp)
  df$tp_path <- as.numeric(df$tp_path)
  df$fp_path <- as.numeric(df$fp_path)
  return(df)
}
compare_conn(fun_nw_mcor[[1]],str_nw_gt[[1]],str_gt_shortest_path[[1]])
compare_conn_list(fun_nw_mcor,str_nw_gt,str_gt_shortest_path,fun_method="mcor",str_method="gt")

df_fp_change <- rbind(compare_conn_list(fun_nw_mcor,str_nw_gt,str_gt_shortest_path,fun_method="mcor",str_method="gt"),
                      compare_conn_list(fun_nw_mcor,str_nw_lt,str_lt_shortest_path,fun_method="mcor",str_method="lt"),
                      compare_conn_list(fun_nw_mcor,str_nw_st,str_st_shortest_path,fun_method="mcor",str_method="st"),
                      compare_conn_list(fun_nw_pcor,str_nw_gt,str_gt_shortest_path,fun_method="pcor",str_method="gt"),
                      compare_conn_list(fun_nw_pcor,str_nw_lt,str_lt_shortest_path,fun_method="pcor",str_method="lt"),
                      compare_conn_list(fun_nw_pcor,str_nw_st,str_st_shortest_path,fun_method="pcor",str_method="st"),
                      compare_conn_list(fun_nw_pcor_order1,str_nw_gt,str_gt_shortest_path,fun_method="pcor1",str_method="gt"),
                      compare_conn_list(fun_nw_pcor_order1,str_nw_lt,str_lt_shortest_path,fun_method="pcor1",str_method="lt"),
                      compare_conn_list(fun_nw_pcor_order1,str_nw_st,str_st_shortest_path,fun_method="pcor1",str_method="st"),
                      compare_conn_list(fun_nw_gl,str_nw_gt,str_gt_shortest_path,fun_method="gl",str_method="gt",fnames=names(fun_gl_shortest_path)),
                      compare_conn_list(fun_nw_gl,str_nw_lt,str_lt_shortest_path,fun_method="gl",str_method="lt",fnames=names(fun_gl_shortest_path)),
                      compare_conn_list(fun_nw_gl,str_nw_st,str_st_shortest_path,fun_method="gl",str_method="st",fnames=names(fun_gl_shortest_path)))
df_fp_change[which(df_fp_change$fun_method=="gl"),"fun_nw"] <- paste0("mag_",df_fp_change[which(df_fp_change$fun_method=="gl"),"fun_nw"])


ggplot(subset(df_fp_change, str_method=="gt"&str_nw=="sum_and_thre"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "Global threshold (sum_and_thre) - number of false positive edges") +
  theme_minimal() + facet_wrap(~ fun_method)
ggplot(subset(df_fp_change, str_method=="gt"&str_nw=="thre_and_sum"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "Global threshold (thre_and_sum) - number of false positive edges") +
  theme_minimal() + facet_wrap(~ fun_method)
ggplot(subset(df_fp_change, str_method=="lt"&str_nw=="sum_and_thre"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "Local threshold (sum_and_thre) - number of false positive edges") +
  theme_minimal() + facet_wrap(~ fun_method)
ggplot(subset(df_fp_change, str_method=="lt"&str_nw=="thre_and_sum"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "Local threshold (thre_and_sum) - number of false positive edges") +
  theme_minimal() + facet_wrap(~ fun_method)
ggplot(subset(df_fp_change, str_method=="st"&str_nw=="thre_count_p027"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "Sign test (thre_count_p027) - number of false positive edges") +
  theme_minimal() + facet_wrap(~ fun_method)
ggplot(subset(df_fp_change, str_method=="st"&str_nw=="thre_dens_p035"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "Sign test (thre_dens_p035) - number of false positive edges") +
  theme_minimal() + facet_wrap(~ fun_method)

ggplot(subset(df_fp_change, fun_method=="mcor"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "marginal cor - number of false positive edges") +
  theme_minimal() + facet_wrap(~ str_method)
ggplot(subset(df_fp_change, fun_method=="pcor"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "partial cor - number of false positive edges") +
  theme_minimal() + facet_wrap(~ str_method)
ggplot(subset(df_fp_change, fun_method=="pcor1"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "order1 pcor - number of false positive edges") +
  theme_minimal() + facet_wrap(~ str_method)
ggplot(subset(df_fp_change, fun_method=="gl"), 
       aes(x = factor(1), y = fp, group = fun_nw, color = fun_nw)) +
  geom_point(size = 1.5) + geom_point(aes(x = factor(2), y = fp_path), size = 1.5) +
  geom_segment(aes(x = factor(1), xend = factor(2), y = fp, yend = fp_path)) +
  scale_x_discrete(labels = c("fp", "fp_path")) +
  labs(x="",y="",title = "glasso - number of false positive edges") +
  theme_minimal() + facet_wrap(~ str_method)


## 3.1 ED&FP 136-edge networks ####
### mcor ####
compute_ed_fp_list(fun_nw_mcor,str_nw_gt)
compute_ed_fp_list(fun_nw_mcor,str_nw_lt)
compute_ed_fp_list(fun_nw_mcor,str_nw_st, snames=c("thre_count_p027","thre_dens_p035"))

### pcor ####
compute_ed_fp_list(fun_nw_pcor,str_nw_gt)
compute_ed_fp_list(fun_nw_pcor,str_nw_lt)
compute_ed_fp_list(fun_nw_pcor,str_nw_st, snames=c("thre_count_p027","thre_dens_p035"))

### pcor_order1 ####
compute_ed_fp_list(fun_nw_pcor_order1,str_nw_gt)
compute_ed_fp_list(fun_nw_pcor_order1,str_nw_lt)
compute_ed_fp_list(fun_nw_pcor_order1,str_nw_st, snames=c("thre_count_p027","thre_dens_p035"))
### gl ####
compute_ed_fp_list(fun_nw_gl,str_nw_gt,fnames=c("all","sub6"))
compute_ed_fp_list(fun_nw_gl,str_nw_lt,fnames=c("all","sub6"))
compute_ed_fp_list(fun_nw_gl,str_nw_st,fnames=c("all","sub6"), snames=c("thre_count_p027","thre_dens_p035"))

### pheatmap ####
ed1200 <- matrix(c(
  248, 246, 262, 256, 246, 248,
  243, 241, 255, 249, 241, 241,
  246, 244, 252, 246, 244, 244,
  246, 244, 250, 244, 244, 244,
  248, 246, 250, 246, 246, 246,
  248, 246, 254, 250, 246, 246,
  246, 244, 260, 254, 244, 246
), nrow = 7, byrow = TRUE)
row_names <- c("mcor_mag", "mcor_sig", "pcor_mag", "pcor_sig", "pcor1_mag", "pcor1_sig", "glasso")
col_names <- c("glob_sum_thre", "glob_thre_sum", "loc_sum_thre", "loc_thre_sum", "test_count_p027", "test_dens_p035")
dimnames(ed1200) <- list(row_names, col_names)
ed200 <- matrix(c(
  248, 246, 262, 256, 246, 248,
  248, 246, 258, 252, 246, 246,
  244, 242, 250, 244, 242, 242,
  248, 246, 252, 248, 246, 246,
  248, 246, 250, 246, 246, 246,
  250, 248, 252, 248, 248, 248,
  246, 244, 260, 254, 244, 246
), nrow = 7, byrow = TRUE)
dimnames(ed200) <- list(row_names, col_names)
pheatmap::pheatmap(ed1200, Rowv = NA, Colv = NA, scale = "none", main="Edit dist, 1200 time points",
                   col = heat.colors(256), margins = c(6, 6),cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(ed200, Rowv = NA, Colv = NA, scale = "none", main="Edit dist, 200 time points",
                   col = heat.colors(256), margins = c(6, 6),cluster_rows = F,cluster_cols = F)


## 3.2 FP path connection ####
### 3.2.1 dataframe ####
### mcor
df_pathfp_mcor <- rbind(compute_fp_list_bin(fun_nw_mcor,str_gt_shortest_path),
      compute_fp_list_bin(fun_nw_mcor,str_lt_shortest_path), # not meaningful
      compute_fp_list_bin(fun_nw_mcor,str_st_shortest_path, snames=c("thre_count_p027","thre_dens_p035")))
df_pathfp_mcor$fun_method <- "mcor"
df_pathfp_mcor$str_method <- rep(c("gt","lt","st"),each=8)

### pcor
df_pathfp_pcor <- rbind(compute_fp_list_bin(fun_nw_pcor,str_gt_shortest_path),
                        compute_fp_list_bin(fun_nw_pcor,str_lt_shortest_path), # not meaningful
                        compute_fp_list_bin(fun_nw_pcor,str_st_shortest_path, 
                                            snames=c("thre_count_p027","thre_dens_p035")))
df_pathfp_pcor$fun_method <- "pcor"
df_pathfp_pcor$str_method <- rep(c("gt","lt","st"),each=8)
### pcor_order1
df_pathfp_pcor_order1 <- rbind(compute_fp_list_bin(fun_nw_pcor_order1,str_gt_shortest_path),
                               compute_fp_list_bin(fun_nw_pcor_order1,str_lt_shortest_path),
                               compute_fp_list_bin(fun_nw_pcor_order1,str_st_shortest_path,
                                                   snames=c("thre_count_p027","thre_dens_p035")))
df_pathfp_pcor_order1$fun_method <- "pcor_order1"
df_pathfp_pcor_order1$str_method <- rep(c("gt","lt","st"),each=8)
### gl 
df_pathfp_gl <- rbind(compute_fp_list_bin(fun_nw_gl,str_gt_shortest_path,fnames=c("all","sub6")),
                      compute_fp_list_bin(fun_nw_gl,str_lt_shortest_path,fnames=c("all","sub6")),
                      compute_fp_list_bin(fun_nw_gl,str_st_shortest_path,fnames=c("all","sub6"),
                                          snames=c("thre_count_p027","thre_dens_p035")))
df_pathfp_gl$fun_nw <- paste0("mag_",df_pathfp_gl$fun_nw)
df_pathfp_gl$fun_method <- "gl"
df_pathfp_gl$str_method <- rep(c("gt","lt","st"),each=4)

df_pathfp <- rbind(subset(df_pathfp_mcor,str_method!="lt"),
                   subset(df_pathfp_pcor,str_method!="lt"),
                   subset(df_pathfp_pcor_order1,str_method!="lt"),
                   subset(df_pathfp_gl,str_method!="lt"))
for (fun in c("sig_all","sig_sub6")) {
  for (str in c("sum_and_thre","sum_and_thre","thre_count_p027","thre_dens_p035")) {
    df_pathfp <- rbind(df_pathfp,c(fun,str,0,"gl","gt"))
    df_pathfp <- rbind(df_pathfp,c(fun,str,0,"gl","st"))
  }
}
df_pathfp$fun_method <- factor(df_pathfp$fun_method)
df_pathfp$fun_nw <- factor(df_pathfp$fun_nw)
df_pathfp$str_method <- factor(df_pathfp$str_method)
df_pathfp$str_nw <- factor(df_pathfp$str_nw)
df_pathfp$false_positive <- as.numeric(df_pathfp$false_positive)

### 3.2.2 plot ####
## gt
ggplot(subset(df_pathfp,str_method=="gt"&str_nw=="sum_and_thre"), 
       aes(x =fun_method, y = false_positive, fill = fun_nw)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  ggtitle("False positives - FC vs global thre-based SC (path connection)")
## st - thre_count_p027
ggplot(subset(df_pathfp,str_method=="st"&str_nw=="thre_count_p027"), 
       aes(x = fun_method, y = false_positive, fill = fun_nw)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  theme_minimal() +
  ggtitle("False positives - FC vs sign test-based SC (thre_count_p027; path connection)")
## st - thre_dens_p035
ggplot(subset(df_pathfp,str_method=="st"&str_nw=="thre_dens_p035"), 
       aes(x = fun_method, y = false_positive, fill = fun_nw)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() + 
  ggtitle("False positives - FC vs sign test-based SC (thre_dens_p035; path connection)")


## 3.3 FP barplots, path length ####

get_fp_path_length_nonzero <- function(m_str,m_fun) {
  m_fp <- ifelse(m_str>0,0,m_fun) # filter out func connection at which str is disconnected
  val_fp <- m_fp[upper.tri(diag(68))]
  return(val_fp[val_fp>0])
}

get_fp_path_length_df <- function(fun_nw,str_nw) {
  fnames <- names(fun_nw)
  snames <- names(str_nw)
  df <- data.frame()
  for (sname in snames) {
    tabs <- lapply(fun_nw,function(fc) { 
      table(get_fp_path_length_nonzero(m_str=str_nw[[sname]],m_fun=fc))
      })
    max_len <- max(sapply(tabs,function(t) { max(as.numeric(names(t))) } ))
    for (fname in fnames) {
      df_cur <- data.frame(tabs[[fname]])
      if (nrow(df_cur) < max_len) {
        count_miss <- setdiff(1:max_len,as.numeric(names(tabs[[fname]])))
        add_row <- data.frame(cbind(as.character(count_miss),0))
        colnames(add_row) <- c("Var1","Freq")
        df_cur <- rbind(df_cur,add_row)
      }
      df_cur$Var1 <- as.numeric(df_cur$Var1)
      df_cur$Freq <- as.numeric(df_cur$Freq)
      df_cur$FunMethod <- fname
      df_cur$StrMethod <- sname
      df_cur <- df_cur[order(df_cur$Var1),]
      df <- rbind(df,df_cur)
    }
  }
  colnames(df) <- c("PathLength","Freq","FunMethod","StrMethod")
  return(df)
}

get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,str_nw=str_gt_shortest_path)

get_fp_pathlength_freq_barplot <- function(data,title) {
  ggplot(data, aes(x = PathLength, y = Freq, fill = FunMethod)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    theme_minimal() + scale_x_continuous(breaks = 1:7) + ggtitle(title)
}

### FP-path conn ####
#### mcor
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_gt_shortest_path),
                                           StrMethod=="sum_and_thre"),
                               title="False positive edge freq - FC (mar corr) vs SC (global thre, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(mar corr) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(mar corr) vs SC(sign test, thre_dens_p035, path connection)")
#### pcor
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_gt_shortest_path),
                                           StrMethod=="sum_and_thre"),
                               title="False pos edge freq - FC(pcor) vs SC(global thre, path connection)")

get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(pcor) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(pcor) vs SC(sign test, thre_dens_p035, path connection)")

#### pcor_order1
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_gt_shortest_path),
                                           StrMethod=="sum_and_thre"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(global thre, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(sign test, thre_dens_p035, path connection)")

#### glasso
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_gt_shortest_path),
                                           StrMethod=="sum_and_thre"),
                               title="False pos edge freq - FC(glasso) vs SC(global thre, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(glasso) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_st_shortest_path),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(glasso) vs SC(sign test, thre_dens_p035, path connection)")


### FP-136edge conn ####
#### mcor vs gt (sum_and_thre and thre_and_sum give the same results)
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="sum_and_thre"),
                               title="False positive edge freq - FC (mar corr) vs SC (global thre, sum_and_thre, original 136-edge binary connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="thre_and_sum"),
                               title="False positive edge freq - FC (mar corr) vs SC (global thre, thre_and_sum, original 136-edge binary connection)")
### this is very close to the sum_and_thre one, different at locations: 1  5  8 12 15 22

#### mcor vs lt 
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,str_nw=str_nw_lt),
                                           StrMethod=="sum_and_thre"),
                               title="False positive edge freq - FC (mar corr) vs SC (local thre, sum_and_thre, original 136-edge binary connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,str_nw=str_nw_lt),
                                           StrMethod=="thre_and_sum"),
                               title="False positive edge freq - FC (mar corr) vs SC (local thre, thre_and_sum, original 136-edge binary connection)")

#### mcor vs st  
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_nw_st[c(1,3)]),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(mar corr) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_mcor_shortest_path,
                                                                 str_nw=str_nw_st[c(1,3)]),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(mar corr) vs SC(sign test, thre_dens_p035, path connection)")


#### pcor vs gt 
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="sum_and_thre"),
                               title="False pos edge freq - FC(pcor) vs SC(global thre, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="thre_and_sum"),
                               title="False pos edge freq - FC(pcor) vs SC(global thre, path connection)")

#### pcor vs lt 
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,str_nw=str_nw_lt),
                                           StrMethod=="sum_and_thre"),
                               title="False positive edge freq - FC (pcor) vs SC (local thre, sum_and_thre, original 136-edge binary connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,str_nw=str_nw_lt),
                                           StrMethod=="thre_and_sum"),
                               title="False positive edge freq - FC (pcor) vs SC (local thre, thre_and_sum, original 136-edge binary connection)")

#### pcor vs st 
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_nw_st),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(pcor) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_shortest_path,
                                                                 str_nw=str_nw_st),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(pcor) vs SC(sign test, thre_dens_p035, path connection)")

#### pcor_order1 vs gt
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="sum_and_thre"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(global thre, 136-edge connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="thre_and_sum"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(global thre, 136-edge connection)")
#### pcor_order1 vs lt 
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_nw_lt),
                                           StrMethod=="sum_and_thre"),
                               title="False positive edge freq - FC(pcor_order1) vs SC(local thre, sum_and_thre, 136-edge connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_nw_lt),
                                           StrMethod=="thre_and_sum"),
                               title="False positive edge freq - FC(pcor_order1) vs SC(local thre, thre_and_sum, 136-edge connection)")
#### pcor_order1 vs st
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_nw_st),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_pcor_order1_shortest_path,
                                                                 str_nw=str_nw_st),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(pcor_order1) vs SC(sign test, thre_dens_p035, path connection)")

#### gl vs gt
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="sum_and_thre"),
                               title="False pos edge freq - FC(glasso) vs SC(global thre, 136-edge connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_nw_gt),
                                           StrMethod=="thre_and_sum"),
                               title="False pos edge freq - FC(glasso) vs SC(global thre, 136-edge connection)")

#### gl vs lt
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_nw_lt),
                                           StrMethod=="sum_and_thre"),
                               title="False positive edge freq - FC(glasso) vs SC(local thre, sum_and_thre, 136-edge connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_nw_lt),
                                           StrMethod=="thre_and_sum"),
                               title="False positive edge freq - FC(glasso) vs SC(local thre, thre_and_sum, 136-edge connection)")
#### gl vs st
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_nw_st),
                                           StrMethod=="thre_count_p027"),
                               title="False pos edge freq - FC(glasso) vs SC(sign test, thre_count_p027, path connection)")
get_fp_pathlength_freq_barplot(data=subset(get_fp_path_length_df(fun_nw=fun_gl_shortest_path,
                                                                 str_nw=str_nw_st),
                                           StrMethod=="thre_dens_p035"),
                               title="False pos edge freq - FC(glasso) vs SC(sign test, thre_dens_p035, path connection)")


## 3.4 FP mean path length ####
get_fp_path_length_nonzero(str_gt_shortest_path[["sum_and_thre"]],fun_mcor_shortest_path[["mag_all"]])

compute_mean_fp_pathlen <- function(fun_nw,str_nw, snames=c("sum_and_thre","thre_and_sum")) {
  fnames <- names(fun_nw)
  df <- data.frame()
  for (fname in fnames) {
    for (sname in snames) {
      df_cur <- c(fname,sname,mean(get_fp_path_length_nonzero(m_str=str_nw[[sname]],m_fun=fun_nw[[fname]])))
      df <- rbind(df,df_cur)
    }
  }
  colnames(df) <- c("fun_nw","str_nw","mean")
  df$mean <- round(as.numeric(df$mean),3)
  df$fun_nw <- factor(df$fun_nw)
  df$str_nw <- factor(df$str_nw)
  return(df)
}

### plot, path-conn ####
df_fp_pathlen_pathconn <- rbind(compute_mean_fp_pathlen(fun_mcor_shortest_path,str_gt_shortest_path),
                               compute_mean_fp_pathlen(fun_mcor_shortest_path,str_st_shortest_path,
                                                       snames=names(str_st_shortest_path)),
                               compute_mean_fp_pathlen(fun_pcor_shortest_path,str_gt_shortest_path),
                               compute_mean_fp_pathlen(fun_pcor_shortest_path,str_st_shortest_path,
                                                       snames=names(str_st_shortest_path)),
                               compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_gt_shortest_path),
                               compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_st_shortest_path,
                                                       snames=names(str_st_shortest_path)),
                               compute_mean_fp_pathlen(fun_gl_shortest_path,str_gt_shortest_path),
                               compute_mean_fp_pathlen(fun_gl_shortest_path,str_st_shortest_path,
                                                       snames=names(str_st_shortest_path)) )
df_fp_pathlen_pathconn$fun_method <- c(rep(c("mcor","pcor","pcor_order1"),each=16),rep("glasso",8))
df_fp_pathlen_pathconn$str_method <- c(rep(rep(c("gt","st"),each=8),3),rep(c("gt","st"),each=4))
df_fp_pathlen_pathconn[which(df_fp_pathlen_pathconn$fun_method=="glasso"),"fun_nw"] <- paste0("mag_",df_fp_pathlen_pathconn[which(df_fp_pathlen_pathconn$fun_method=="glasso"),"fun_nw"])

## FPmean_pathconn_gt
ggplot(subset(df_fp_pathlen_pathconn,str_method=="gt"), aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(global thre, path conn)")
ggplot(subset(df_fp_pathlen_pathconn,str_method=="st"&str_nw=="thre_count_p027"), 
       aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(sign test, thre_count_p027, path conn)")
ggplot(subset(df_fp_pathlen_pathconn,str_method=="st"&str_nw=="thre_dens_p035"), 
       aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(sign test, thre_dens_p035, path conn)")


compute_mean_fp_pathlen(fun_mcor_shortest_path,str_gt_shortest_path)
compute_mean_fp_pathlen(fun_mcor_shortest_path,str_st_shortest_path,snames=names(str_st_shortest_path))
compute_mean_fp_pathlen(fun_pcor_shortest_path,str_gt_shortest_path)
compute_mean_fp_pathlen(fun_pcor_shortest_path,str_st_shortest_path,snames=names(str_st_shortest_path))
compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_gt_shortest_path)
compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_st_shortest_path,snames=names(str_st_shortest_path))
compute_mean_fp_pathlen(fun_gl_shortest_path,str_gt_shortest_path)
compute_mean_fp_pathlen(fun_gl_shortest_path,str_st_shortest_path,snames=names(str_st_shortest_path))


### 136edge-conn ####
df_fp_pathlen_136conn <- rbind(compute_mean_fp_pathlen(fun_mcor_shortest_path,str_nw_gt),
                               compute_mean_fp_pathlen(fun_mcor_shortest_path,str_nw_lt),
                               compute_mean_fp_pathlen(fun_mcor_shortest_path,str_nw_st,
                                                       snames=names(str_st_shortest_path)),
                               compute_mean_fp_pathlen(fun_pcor_shortest_path,str_nw_gt),
                               compute_mean_fp_pathlen(fun_pcor_shortest_path,str_nw_lt),
                               compute_mean_fp_pathlen(fun_pcor_shortest_path,str_nw_st,
                                                       snames=names(str_st_shortest_path)),
                               compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_nw_gt),
                               compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_nw_lt),
                               compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_nw_st,
                                                       snames=names(str_st_shortest_path)),
                               compute_mean_fp_pathlen(fun_gl_shortest_path,str_nw_gt),
                               compute_mean_fp_pathlen(fun_gl_shortest_path,str_nw_lt),
                               compute_mean_fp_pathlen(fun_gl_shortest_path,str_nw_st,
                                                       snames=names(str_st_shortest_path)) )
df_fp_pathlen_136conn$fun_method <- c(rep(c("mcor","pcor","pcor_order1"),each=24),rep("glasso",12))
df_fp_pathlen_136conn$str_method <- c(rep(rep(c("gt","lt","st"),each=8),3),rep(c("gt","lt","st"),each=4))
df_fp_pathlen_136conn[which(df_fp_pathlen_136conn$fun_method=="glasso"),"fun_nw"] <- paste0("mag_",df_fp_pathlen_136conn[which(df_fp_pathlen_136conn$fun_method=="glasso"),"fun_nw"])


compute_mean_fp_pathlen(fun_mcor_shortest_path,str_nw_gt)
compute_mean_fp_pathlen(fun_mcor_shortest_path,str_nw_lt)
compute_mean_fp_pathlen(fun_mcor_shortest_path,str_nw_st,snames=names(str_st_shortest_path))
compute_mean_fp_pathlen(fun_pcor_shortest_path,str_nw_gt)
compute_mean_fp_pathlen(fun_pcor_shortest_path,str_nw_lt)
compute_mean_fp_pathlen(fun_pcor_shortest_path,str_nw_st,snames=names(str_st_shortest_path))
compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_nw_gt)
compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_nw_lt)
compute_mean_fp_pathlen(fun_pcor_order1_shortest_path,str_nw_st,snames=names(str_st_shortest_path))
compute_mean_fp_pathlen(fun_gl_shortest_path,str_nw_gt)
compute_mean_fp_pathlen(fun_gl_shortest_path,str_nw_lt)
compute_mean_fp_pathlen(fun_gl_shortest_path,str_nw_st,snames=names(str_st_shortest_path))

# FPmean_136conn_gt
ggplot(subset(df_fp_pathlen_136conn,str_method=="gt"), aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(global thre, 136-edge conn)")
ggplot(subset(df_fp_pathlen_136conn,str_method=="lt"), aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(local thre, 136-edge conn)")
ggplot(subset(df_fp_pathlen_136conn,str_method=="st"&str_nw=="thre_count_p027"), aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(sign test, thre_count_p027, 136-edge conn)")
ggplot(subset(df_fp_pathlen_136conn,str_method=="st"&str_nw=="thre_dens_p035"), aes(x = fun_method, y = mean, fill =fun_nw )) +
  geom_bar(stat = "identity", position = position_dodge()) + theme_minimal() +
  ggtitle("False positive mean shortest path length, FC vs SC(sign test, thre_dens_p035, 136-edge conn)")















