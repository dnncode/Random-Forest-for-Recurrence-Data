# ---------------------------------------------------------------
# ---------------------------------------------------------------
# The code is based on the paper:
#"Analysis of Large Heterogeneous Repairable System Reliability Data with Static System Attributes and Dynamic Sensor Measurement in Big Data Environment"
# Available at ArXiv: ***

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# RF-R: Random Forests for Repairable System Reliability Analysis
# This version of code handles both the static system attributes and dynamic sensor measurement
# (RF-R-MCF.R is the code that handles only the static sytem attributes)
# Date: 12/11/2018
# ---------------------------------------------------------------
# Dependencies: SNOW (Simple Network of Workstations) package in R for simple parallel computing
# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------

rm(list=ls())
t.start = proc.time() # record the start time
wd= "C:/Users/xl027/Desktop/24. SurvivalTree/data"
setwd(wd)

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# User specified parameters
# ---------------------------------------------------------------
# ---------------------------------------------------------------
C.index = TRUE 
# If TRUE, the program computes OOB prediction C-index 
# the default and recommended setting is TRUE

n.tree = 500 # number of trees to be grown
B.ratio = 0.75 # percentage of samples to be included in a bootstrap sample
n.OOB.sample = 100 
# when evaluating the OOB prediction C-index, sometimes it is possible
# to perform the calculation based on a randomly selected subset of the
# OOB samples (n.OOB.sample specifies the size of the random subset).
# This reduces the runtime of the code, but it is not always necessary.

n.node = 32 # Number of threads on a single computing node for parallel computing
            # If necessary, run multiple procesess in parallel on a computing cluster (HPC is highly recommended for large data sets)
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------
# ---------------------------------------------------------------
CaseStudy1 = TRUE  
# The Case Study data is used as an illustration
# Hence, CaseStudy1 is set to TRUE
if (CaseStudy1){
  data = read.csv("CaseStudy1/data.csv",header=TRUE,stringsAsFactors = FALSE)
  # data: a matrix that contains the failure data;
  # nrow(data) = number of systems
  # In each row, the first value is the installation time (default value:0)
  # In each row, the last non-negative value is the censoring time
  # In each row, "-1" is just a place holder
  X = read.csv("CaseStudy1/x.table.csv",header=TRUE,stringsAsFactors = FALSE)
  # X: a matrix that contains static system attributes
  # dim(X): number of systems X number of system attributes
  load("CaseStudy1/z.list.RData")  # z.list is loaded
  # z.list is a list that contains the dynamic sensor measurement
  # dim(z.list): number of systems
  # Each element of z.list is a matrix
  # The dimension of the matrix above is: number of measurement by each sensor X the number of sensor types 
} 
# Note: all X and Z have been standardized on [0,1]

# Remove systems with incomplete system attributes information
# In general, it does not make sense to impute system attributes
case = which(!is.na(rowSums(X)))
data=data[case,]; X = X[case,]; z.list=z.list[case]

n = nrow(data)  # number of systems available
#dimension = ncol(z.list[[1]])+1+ncol(X)  # dimension of the log-linear rate function; not used

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Load sub-routines
# ---------------------------------------------------------------
# ---------------------------------------------------------------
source("../code/nhpp_function_2.R")  # perform regularized (Lasso) MLE to fit the NHPP model
source("../code/nhpp_c_2.R") # calculate the distance between the fitted NHPP between two daughter nodes
source("../code/tree_split_nhpp_2.R") # for the Greedy algorithm, calculates the distance measure between two daughter nodes
source("../code/tree_grow_nhpp.R") 
# based on the results returned from "tree_split_nhpp_2.R", perform node split and update the tree structure
source("../code/terminal_test.R") # tree node split termination test
source("../code/tree_predict_nhpp.R") # for a given tree in the forest, perform the prediction for given input
source("../code/tree_predict_Cindex_nhpp_2.R") # compute the OOB prediction C-index
source("../code/tree_predict_ensemble_nhpp_2.R") # perform the ensemble prediction
source("../code/terminal_mle.R") 
# for each terminal node of a tree, perform regularized (Lasso) MLE to fit the NHPP model

library(snow) # use the SNOW (Simple Network of Workstations) package in R for simple parallel computing
cl <- makeCluster(n.node, type = "SOCK") # create computing cluster
# load all sub-routines on all nodes
clusterEvalQ(cl, source("../code/nhpp_function_2.R"))
clusterEvalQ(cl, source("../code/nhpp_c_2.R"))
clusterEvalQ(cl, source("../code/tree_split_nhpp_2.R"))
clusterEvalQ(cl, source("../code/tree_grow_nhpp.R"))
clusterEvalQ(cl, source("../code/terminal_test.R"))
clusterEvalQ(cl, source("../code/tree_predict_nhpp.R"))
clusterEvalQ(cl, source("../code/tree_predict_Cindex_nhpp_2.R"))
clusterEvalQ(cl, source("../code/tree_predict_ensemble_nhpp_2.R"))
clusterEvalQ(cl, source("../code/terminal_mle.R"))

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Basic settings for RF
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# number of covariates randomly selected for each split of a tree node
# the users can certainly change the current settings
if (ncol(X)<=2){
  X.ratio = ncol(X)
}else{
  X.ratio = max( max(floor( ncol(X)/3 ), 1), 2 )
}
if (X.ratio>4){X.ratio=4}

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Running the RF-R algorithm
# ---------------------------------------------------------------
# ---------------------------------------------------------------
tree.list=vector("list", length = n.tree)
# a list that contains all information for a tree in the forests
# length(tree.list) = number of trees (n.tree)
oob.error.matrix = array(0/0, dim=c(   min(n.OOB.sample, n-floor(n*B.ratio)),n.tree))
# a matrix that contains the OOB predictions for all OOB samples and for all trees.
# each column contains the OOB predictions, based on a single tree, for all OOB samples
oob.error = array(0/0, dim=c(n.tree,1))
# a vector that contains the OOB error (in terms of C-index) for each tree.
oob.error.ensemble = array(0/0, dim=c(n.tree,1))
# a vector that contains the cumulative OOB error (in terms of C-index) for ensembles trees.
permute.error = array(0/0, dim=c(n.tree,ncol(X)))
# for variables importance calculation: the OOB error when a particular system attributes are permuted
terminal.mle.list =vector("list", length = n.tree)
# a list that contains the MLE output for all terminals of a tree
# length(terminal.mle.list) = number of trees (n.tree)
B.history = array(0/0, dim=c(floor(n*B.ratio),n.tree))
# record the system IDs which are included in each bootstrap sample

# Now, start the for loop to grow a forest; each loop grows a single tree
for (i.tree in 1:n.tree){
  time1 = proc.time() # record the start time of each loop
  
  # Step 1: Create a bootstrap sample
  B = sample(c(1:n),floor(n*B.ratio),replace = FALSE)
  data.B = data[B,]
  X.B = X[B,]
  Z.B = z.list[B]
  B.history[,i.tree] = B
  
  # Step 2: Grow a tree using the Greedy Algorithm
  
  # Perform the first split
  initial.node = c(-999,-999,0,0)
  # this is how a tree structure is defined in our code:
  # the first entry specifies the split variable (Here, -999 is just a place holder)
  # the second entry specifies the split point (-999 is just a place holder)
  # the third entry specifies the tree branches (1: go to left; 2: go to right; 0: place holder)
  # the fourth entry specifies if this node is a terminal node (1: yes, terminal node; 0: no)
  initial.id = B # IDs of data included in the bootstrap sample
  x.sample = sample(1:ncol(X), X.ratio, replace = FALSE) # generate candidate covariates for splitting
  Data_id = c(1:length(B)) # sequence of data on node
  
  # Calculate the distance between two nodes given candidate split variable and split point
  #time1 = proc.time()
  obj.list = lapply(1:length(x.sample), tree.sub.fun2, Data_id_sub2=Data_id,
                    X.B_sub2=X.B, Z.B_sub2=Z.B, data.B_sub2=data.B, X.sample = x.sample)
  #proc.time()-time1
  # the lapply above can be parallelized. length(x.sample) is determined by X.ratio, which is not bigger than 4 specified by line 117
  # obj.list = parLapply(cl, 1:length(x.sample), tree.sub.fun2, Data_id_sub2=Data_id,
          # X.B_sub2=X.B, Z.B_sub2=Z.B, data.B_sub2=data.B, X.sample = x.sample) 
  # parallel computing using "snow" package; mcf_c is the computing bottleneck. 
  obj = do.call(cbind,obj.list)
  
  # Select the optimum split variable and split point --> Split the node
  x.id.opt = which(apply(obj,2,max,na.rm=TRUE)== max(apply(obj,2,max,na.rm=TRUE),na.rm=TRUE))[1]
  x.opt = x.sample[x.id.opt]
  x.target = X.B[Data_id,
                 x.opt]
  v.min = min(x.target,na.rm=TRUE)
  v.max = max(x.target,na.rm=TRUE)  
  x.delta = (v.max-v.min)/20
  x.seq = seq(v.min+3*x.delta,v.max-3*x.delta,x.delta)
  thres.opt = mean( x.seq[which(obj[,x.id.opt]==max(apply(obj,2,max,na.rm=TRUE),na.rm=TRUE))],na.rm=TRUE )
  daughter.l.id = Data_id[ which(X.B[Data_id,x.opt]<=thres.opt) ]
  daughter.r.id = Data_id[ which(X.B[Data_id,x.opt]>thres.opt) ]
  # Test if the left and right daughter nodes are terminal nodes
  terminal.l = terminal.test(data.B[daughter.l.id,])
  terminal.r = terminal.test(data.B[daughter.r.id,])
  # Update the tree structure
  terminal.list = list()
  n.terminal = 2 # number of nodes
  tmp = list()
  tmp[[1]] = rbind(initial.node, c(x.opt,thres.opt,1,terminal.l))
  tmp[[2]] = daughter.l.id
  terminal.list[[1]] = tmp
  tmp[[1]] = rbind(initial.node, c(x.opt,thres.opt,2,terminal.r))
  tmp[[2]] = daughter.r.id
  terminal.list[[2]] = tmp
  print(c(i.tree,1))
  
  # After the first node split, start the splitting iteration
  stop.rule = 0
  tree.level = 1 # current depth of a tree
  while (stop.rule == 0){ # stop.rule == 1 if all nodes are terminal nodes
    print("start lapply_tree.grow") # current depth of a tree
    tmp.list = lapply(1:n.terminal, tree.grow, 
                      Data_input=data.B, X.B_input = X.B, Z.B_input = Z.B,Terminal.list_input=terminal.list)
    #tmp.list = parLapply(cl, 1:n.terminal, tree.grow, 
    #                 Data_input=data.B, X.B_input = X.B, Z.B_input = Z.B,
    #			Terminal.list_input=terminal.list) 
    # cannot use parLapply here: tree.grow calls tree.sub.fun2 which replies on snow.
    print("complete lapply_tree.grow")
    terminal.list = do.call(c, tmp.list)
    n.terminal = length(terminal.list)
    
    sum = 0 
    for (i.node in 1:n.terminal){
      sum = sum + terminal.list[[i.node]][[1]][ nrow( terminal.list[[i.node]][[1]] ) , 4 ] 
      #print(terminal.list[[i.node]][[1]][ nrow( terminal.list[[i.node]][[1]] ) , 4 ] )
    }
    if (sum >= n.terminal){stop.rule=1}
    tree.level = tree.level + 1
    print(c(i.tree,tree.level))
  }
  
  tree.list[[i.tree]] = terminal.list # terminal.list records the tree structure of the current tree
  terminal.mle.tree = terminal_mle(terminal_mle_input=terminal.list, 
                                   X.train=X.B, Z.train=Z.B, data.train=data.B)
  # perform MLE at each terminal node
  terminal.mle.list[[i.tree]] = terminal.mle.tree
  
  
  # Step 3: generate OOB prediction C-index and evaluate covariate importance
  # (Steps 1 and 2 complete the tree growing process)
  # variable selection
  data.OOB = data[-B,]
  X.OOB = X[-B,]
  Z.OOB = z.list[-B]
  if (nrow(data.OOB)>n.OOB.sample){ 
    OOB.sample = sample(1:nrow(data.OOB),n.OOB.sample,replace=FALSE) 
    X.OOB =  X.OOB[OOB.sample, ]
    data.OOB = data.OOB[OOB.sample, ]
    Z.OOB = Z.OOB[OOB.sample]
  }
  
  if (!C.index){  # the current version always set C.index to TRUE, hence, !C.index=FALSE
    
    #oob.error.matrix[,i.tree] = sapply(1:nrow(X.OOB),tree_predict, 
    #                                   X.pred = X.OOB, tree.trained= terminal.list, 
    #							data.train= data.B, data.pred = data.OOB,
    # 							Z.pred = Z.OOB, Z.train=Z.B)
    oob.error.matrix[,i.tree] = parSapply(cl, 1:nrow(X.OOB), tree_predict, 
                                          X.pred = X.OOB, tree.trained= terminal.list, 
                                          data.train= data.B, data.pred = data.OOB,
                                          Z.pred = Z.OOB, Z.train=Z.B) 
    case = which(oob.error.matrix[,i.tree]==Inf)
    oob.error.matrix[case,i.tree] = 0/0
    oob.error[i.tree,1] = mean(oob.error.matrix[,i.tree],na.rm=TRUE)
  }else{
    #tmp = lapply(1:nrow(X.OOB),tree_predict_C, 
    #                                 X.pred = X.OOB, tree.trained= terminal.list, 
    #						data.train= data.B, data.pred = data.OOB,
    #						Z.pred = Z.OOB, Z.train=Z.B)
    tmp = parLapply(cl, 1:nrow(X.OOB), tree_predict_C, 
                    X.train=X.B,
                    X.pred = X.OOB, tree.trained= terminal.list, 
                    data.train= data.B, data.pred = data.OOB,
                    Z.pred = Z.OOB, Z.train=Z.B, tree.mle=terminal.mle.tree )
    tmp = do.call(rbind, tmp)
    pairs = t( combn(c(1:nrow(X.OOB)), 2) ) # create all pairs from OOB sample
    C.tmp = parSapply(cl, 1:nrow(pairs), tree_predict_C_index, 
                      Pairs = pairs, Data.input=tmp) # compute C-index
    #C.tmp = sapply(1:nrow(pairs), tree_predict_C_index, 
    #                                   Pairs = pairs, Data.input=tmp)
    oob.error[i.tree,1] = 1-mean(C.tmp, na.rm=TRUE) # error: 1-C-index
  }
  
  # Compute C-index after permuting the values for each covariate (for variable selectin purposes)
  for (j in 1:ncol(X)){
    X.OOB.permute = X.OOB
    X.OOB.permute[,j] = sample(X.OOB[,j], nrow(X.OOB), replace = FALSE)
    if (!C.index){ # the current version always set C.index to TRUE, hence, !C.index=FALSE
      #tmp = sapply(1:nrow(X.OOB),tree_predict, 
      #             X.pred = X.OOB.permute, tree.trained= terminal.list, data.train= data.B, data.pred = data.OOB,
      #  	     Z.pred = Z.OOB, Z.train=Z.B)
      tmp =  parSapply(cl, 1:nrow(X.OOB), tree_predict, 
                       X.pred = X.OOB.permute, tree.trained= terminal.list, data.train= data.B, data.pred = data.OOB,
                       Z.pred = Z.OOB, Z.train=Z.B) 
      case = which(tmp==Inf)
      tmp[case] = 0/0
      permute.error[i.tree,j] = mean(tmp,na.rm=TRUE)-oob.error[i.tree,1]
    }else{
      tmp = parLapply(cl, c(1:nrow(X.OOB)), tree_predict_C, X.train=X.B,
                      X.pred = X.OOB.permute, tree.trained= terminal.list, 
                      data.train= data.B, data.pred = data.OOB,Z.pred = Z.OOB, Z.train=Z.B, tree.mle=terminal.mle.tree)
      #tmp = lapply(1:nrow(X.OOB), tree_predict_C, 
      #                            X.pred = X.OOB.permute, tree.trained= terminal.list, 
      #				     data.train= data.B, data.pred = data.OOB,Z.pred = Z.OOB, Z.train=Z.B)
      tmp = do.call(rbind, tmp)
      pairs = t( combn(c(1:nrow(X.OOB)), 2) )
      C.tmp = parSapply(cl, 1:nrow(pairs), tree_predict_C_index, 
                        Pairs = pairs, Data.input=tmp)
      #C.tmp = sapply(1:nrow(pairs), tree_predict_C_index, 
      #                             Pairs = pairs, Data.input=tmp)
      permute.error[i.tree,j]   = (1-mean(C.tmp, na.rm=TRUE))- oob.error[i.tree,1] 
    }
    print(j)
  }
  
  # Track ensemble OOB error as the number of trees grows
  output =  parLapply(cl, 1:i.tree, tree_predict_ensemble,
                      X.train.ensemble = X.B,
                      X.pred.ensemble = X.OOB, 
                      tree.trained.ensemble= tree.list, 
                      data.input = data, 
                      Z.input = z.list,	
                      data.pred.ensemble = data.OOB,
                      Z.pred.ensemble = Z.OOB,
                      B.hist=B.history, tree.mle.list=terminal.mle.list)   # output gives the ensemble rate prediction 
  output = do.call(cbind, output)
  pred.ensemble = apply(output, 1, mean, na.rm=TRUE)
  
  obs= array()
  for (i in 1:nrow(data.OOB)){
    n.failure = length(which(data.OOB[i,]>0))-1
    t.length = data.OOB[i,n.failure+1]
    obs[i] = n.failure/t.length
  }
  if (!C.index){
    case.keep = which(obs>0)
    oob.error.ensemble[i.tree] = mean ( abs(obs[case.keep] - pred.ensemble[case.keep])/obs[case.keep], na.rm=TRUE) 
  }else{
    tmp = cbind(pred.ensemble, obs)
    pairs = t( combn(c(1:nrow(tmp)), 2) )
    C.tmp = parSapply(cl, 1:nrow(pairs), tree_predict_C_index, 
                      Pairs = pairs, Data.input=tmp)
    #C.tmp = sapply(1:nrow(pairs), tree_predict_C_index, 
    #                           Pairs = pairs, Data.input=tmp)
    oob.error.ensemble[i.tree]   = 1-mean(C.tmp, na.rm=TRUE)
  }
  
  
  print(oob.error.ensemble[1:i.tree]) # print the ensemble OOB error as the number of trees grows
  
  # plot "ensemble OOB error as the number of trees grows" and "variable importance" on-line
  tmp = apply(permute.error, 2, mean, na.rm=TRUE)*100
  par(mfrow=c(2,1))
  par(las=TRUE)
  par(mar=c(4,4,1,1))
  print(barplot(tmp, names.arg=paste("", ncol(X), sep=" "), 
                col="blue", ylab="variable importance"))
  print(plot(oob.error.ensemble[1:i.tree],type="s"))
  
  print(proc.time()- time1) # print the runtime for growing this tree
  
}
t.end = proc.time() - t.start # print the total runtime for growing this tree 

save.image("CaseStudy1/output_Cindex.RData") # save the output 
stopCluster(cl) # close the computing cluster


# --------------- END -----------------------







