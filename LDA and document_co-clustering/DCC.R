###########################################################
###########################################################
#################### DCC Implementation ###################
###########################################################
###########################################################

library(dplyr)
library(janeaustenr)
library(tidytext)
library(R.matlab)
library(mclust)
library(Matrix)

#########################################################
# Copyright (C) 2016  Aghiles Salah. All rights reserved.
# License: Apache 2.0 License.
#########################################################

# Normalized Mutual Information
NMI = function(row_cluster, true_labels, n){
  N_kl = as.matrix(table(row_cluster, true_labels))
  N_k = rowSums(N_kl)
  N_l = colSums(N_kl)
  Num = 0
  denum_k = 0
  for (i in 1 : nrow(N_kl)) {
    denum_k = denum_k + (N_k[i] / n) * log((N_k[i] / n))
    for (j in 1 : ncol(N_kl)) {
      if (N_kl[i, j] != 0) {
        num = log(n) + log(N_kl[i, j])
        den = log(N_k[i]) + log(N_l[j])
        Num = Num + (N_kl[i, j] / n) * (num - den)
      }
    }
  }
  denum_l = 0
  for (j in 1 : ncol(N_kl)) {
    denum_l = denum_l + (N_l[j] / n) * log((N_l[j] / n))
  }
  resnmi = Num / (sqrt(denum_k * denum_l))
  resnmi
}

# Partition generator
par_gen = function(n, k, nb_par=1){
  
  par = as.integer(sample(as.numeric(1 : k), n, replace = TRUE))
  if (nb_par > 1) {
    for (i in 2 : nb_par) {
      par = rbind(par, as.integer(sample(as.numeric(1 : k), n, replace = TRUE)))
    }
    par = as.matrix(par)
  }
  par
}

# TF-IDF data representation
tf_idf = function(sparse_mat, l2_norm = FALSE){
  
  bin_mat = replace(sparse_mat, sparse_mat > 0, 1)
  tfidf_mat = sparse_mat + sparse_mat * log(1 + nrow(sparse_mat)) - t(t(sparse_mat) * log((1 + colSums(bin_mat))))
  
  if (l2_norm)
    tfidf_mat = tfidf_mat / sqrt(rowSums(tfidf_mat * tfidf_mat))
  
  tfidf_sp_mat = as(tfidf_mat, "dgCMatrix")
  tfidf_sp_mat
}

##########################################################
# Copyright (C) 2016  Aghiles Salah. All rights reserved.
# License: Apache 2.0 License.
##########################################################

dcc_sample = function(prob){
  sample(1 : length(prob), size = 1, replace = FALSE, prob = prob)
}

dcc = function(X, k, iter.max=100, stoch_iter.max=70, row_init=NULL, col_init=NULL, n_init=5, tol=1e-6){
  
  # Coherence tests
  
  if (stoch_iter.max >= iter.max) {
    stop("Erorr stoch_iter.max must be less than iter.max")
  }
  
  if (k >= min(dim(X))) {
    stop("More clusters than distinc rows/columns")
  }
  
  if (n_init > 1) {
    # coherence tests for row_init
    if (! is.null(row_init)) {
      if (is.matrix(row_init)) {
        if (dim(row_init)[1] < n_init)
          stop("Less row partitions than n_init")
        if (dim(row_init)[2] != nrow(X))
          stop("The length of the row partitions is different from the number of objects")
      }
      else
        stop("Error o_O', row_init is not a matrix")
    }
    
    # coherence tests for col_init
    if (! is.null(col_init)) {
      if (is.matrix(col_init)) {
        if (dim(col_init)[1] < n_init)
          stop("Less column partitions than n_init")
        if (dim(col_init)[2] != ncol(X))
          stop("The length of the col partitions is different from the number of columns")
      }
      else
        stop("Error o_O', col_init is not a matrix")
    }
  }
  else {
    if (! is.null(row_init)) {
      if (length(row_init) != nrow(X))
        stop("The length of the row partition must be equal to the number of objects")
    }
    if (! is.null(col_init)) {
      if (length(col_init) != ncol(X))
        stop("The length of the column partition must be equal to the number of columns")
    }
  }
  
  # normalize rows to have unit L2 norm
  X = as(X, "dgCMatrix")
  X = X / sqrt(rowSums(X * X))
  n = nrow(X)
  p = ncol(X)
  vtw = c(rep(0, iter.max))
  
  # perform one run
  do_one <- function(row_c, col_c){
    
    vtw = c(rep(0, iter.max))
    nbIter = iter.max
    
    # Build initial binary row-cluster indicator matrix
    Z = matrix(0, n, k)
    Z = as(Z, "dgCMatrix")
    Z[cbind(seq_along(row_c), row_c)] = 1
    
    # Build initial binary column-cluster indicator matrix
    W = matrix(0, p, k)
    W = as(W, "dgCMatrix")
    W[cbind(seq_along(col_c), col_c)] = 1
    
    
    # Compute initial row centroids
    MU_w = diag(1 / sqrt(table(col_c)))
    MU_w = as(MU_w, "dgCMatrix")
    
    # Compute initial column centroids
    MU_z = diag(1 / sqrt(table(row_c)))
    MU_z = as(MU_z, "dgCMatrix")
    
    #DCC alternating optimization
    for (iter in 1 : iter.max) {
      
      # row partitionning
      Zt = (X %*% W) %*% (MU_w %*% MU_z)
      Zt = as(Zt, "dgCMatrix")
      row_partition = apply(Zt, 1, which.max)
      Z = matrix(0, n, k)
      Z = as(Z, "dgCMatrix")
      Z[cbind(seq_along(row_partition), row_partition)] = 1
      
      # Update column centroids MU^z
      MU_z = diag(1 / sqrt(table(row_partition)))
      MU_z = as(MU_z, "dgCMatrix")
      
      # Column partitionning
      Wt = crossprod(X, (Z %*% (MU_z %*% MU_w)))
      Wt = as(Wt, "dgCMatrix")
      
      if (iter <= stoch_iter.max) {
        # Perform stochastic column assignement to avoid bad local solutions
        col_partition = apply(Wt, 1, dcc_sample)
      }
      else {
        col_partition = apply(Wt, 1, which.max)
      }
      
      W = matrix(0, p, k)
      W = as(W, "dgCMatrix")
      W[cbind(seq_along(col_partition), col_partition)] = 1
      
      # Update row centroids MU^w
      MU_w = diag(1 / sqrt(table(col_partition)))
      MU_w = as(MU_w, "dgCMatrix")
      
      
      # evaluate the criterion
      vtw[iter] = sum(Z * ((X %*% W) %*% (MU_w %*% MU_z)))
      
      if (iter > 1) {
        if (abs(vtw[iter] - vtw[iter - 1]) < tol) {
          nbIter = iter
          break
        }
      }
    }
    
    
    structure(list(rowcluster = row_partition, colcluster = col_partition, ll = vtw, iter = nbIter))
    
    #End of do_one
  }
  
  # Preparing initial row partition(s)
  if (is.null(row_init)) {
    row_c = as.integer(sample(as.numeric(1 : k), n, replace = TRUE))
  }
  else {
    if (is.matrix(row_init)) {
      row_c = as.integer(row_init[1,])
    }
    else {
      row_c = row_init
    }
  }
  
  # Preparing initial column partition(s)
  if (is.null(col_init)) {
    col_c = as.integer(sample(as.numeric(1 : k), p, replace = TRUE))
  }
  else {
    if (is.matrix(col_init)) {
      col_c = as.integer(col_init[1,])
    }
    else {
      col_c = col_init
    }
  }
  
  ###
  Run = do_one(row_c, col_c)
  best = Run$ll[Run$iter]
  index_best = 1;
  if (n_init >= 2) {
    bool = TRUE
    #Coherence tests
    if (is.matrix(row_init)) {
      if (dim(row_init)[1] != n_init) {
        bool = FALSE
      }
    }
    if (is.matrix(col_init)) {
      if (dim(col_init)[1] != n_init) {
        bool = FALSE
      }
    }
    if (bool) {
      for (i in 2 : n_init) {
        # row partition preparation
        if (is.null(row_init)) {
          row_c = as.integer(sample(as.numeric(1 : k), n, replace = TRUE))
        }
        else {
          if (is.matrix(row_init)) {
            row_c = as.integer(row_init[i,])
          }
          else {
            stop("Error o_O', row_init is not a matrix")
          }
        }
        # column partition preparation
        if (is.null(col_init)) {
          col_c = as.integer(sample(as.numeric(1 : k), p, replace = TRUE))
        }
        else {
          if (is.matrix(col_init)) {
            col_c = as.integer(col_init[i,])
          }
          else {
            stop("Error o_O', col_init is not a matrix")
          }
        }
        
        RRun <- do_one(row_c, col_c)
        if (! is.na(RRun$ll[RRun$iter])) {
          if (RRun$ll[RRun$iter] > best) {
            Run <- RRun
            best <- Run$ll[Run$iter]
            index_best = i
          }
        }
      }
    }
    else {
      stop("Error o_O', the number of row/column partitions do not equal n_init ")
    }
  }
  Run$index_best = index_best
  Run
}

##############################################################################

# Import Data
setwd("C:/Users/jueju/Desktop/LDA and Co-clustering algorithms with data")
sparse_matrix <- read.csv("doc_word_matrix_stemmingf.csv", header = FALSE)
b = as.matrix(sparse_matrix)
colnames(b) = NULL
B <- Matrix(b, sparse = TRUE)

# Collect co-clustering results & obtain running time 
set.seed(88)

k_range  = 1:18
DCC = data.frame(1:nrow(b))
Running_Time = vector()

for (i in 1:length(k_range)){
  num_cluster = as.integer(k_range[i]) + 2
  
  # Running time
  start <- proc.time()
  dcc_result_raw <- dcc(X=B, k=num_cluster, iter.max=100, n_init=5, stoch_iter.max=70)
  time = proc.time() - start
  Running_Time[i] =as.numeric(time[3])
  
  # Store results
  DCC_results = dcc_result_raw$rowcluster
  nn = data.frame(as.numeric(DCC_results))
  DCC = cbind(DCC, nn)
}

Running_Time = as.data.frame(Running_Time)
names(DCC) = c("x", 3:20)

#write.csv(DCC, "DCC_results_stemming.csv")
#write.csv(Running_Time, "DCC_running_time_stemming.csv")



