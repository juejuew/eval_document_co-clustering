#######################################################
#######################################################
############# Model-based Implementation ##############
#######################################################
#######################################################

library(blockcluster)

# Import Data
setwd("C:/Users/jueju/Desktop/LDA and Co-clustering algorithms with data")
sparse_matrix <- read.csv("doc_word_matrix_stemmingf.csv", header = FALSE)
b = as.matrix(sparse_matrix)
colnames(b) = NULL
dim(b)

# Collect co-clustering results & obtain running time
set.seed(88)

k_range  = 1:2
MB = data.frame(1:nrow(b))
Running_Time = vector()

for (i in 1:length(k_range)){
  num_cluster = as.integer(k_range[i]) + 2
  
  # Running time
  start <- proc.time()
  strategy = coclusterStrategy( nbinititerations = 2, nbxem = 2, nbiterations_int = 2
                                , nbiterationsxem = 5, nbiterationsXEM = 50, epsilonXEM=1e-5)
  out<-coclusterContingency(b, nbcocluster=c(num_cluster,4), strategy = strategy)
  time = proc.time() - start
  Running_Time[i] = as.numeric(time[3])
  
  # Store results
  nn = data.frame(as.numeric(out@rowclass))
  MB = cbind(MB, nn)
}

Running_Time = as.data.frame(Running_Time)
names(MB) = c("x", 3:20)

#write.csv(MB, "vx_model_based_stemming.csv")
#write.csv(Running_Time, "model_based_stemming_time.csv")

