# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R.R
# Description: for the Greedy algorithm, calculates the distance measure between two daughter nodes
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# split based on one variable
tree.sub.fun1 = function(i, thres.range, Data_on_Node, Z_on_Node, X_on_Node, covariate){
  thres = thres.range[i]
  group1 = which(covariate<=thres)
  group2 = which(covariate>thres)
  Data_group1 = Data_on_Node[group1, ]
  Data_group2 = Data_on_Node[group2, ]
  Z_group1 = Z_on_Node[group1]
  Z_group2 = Z_on_Node[group2]
  X_group1 = X_on_Node[group1,]
  X_group2 = X_on_Node[group2,]
  if ( (length(group1)>5)&(length(group2)>5) ){
    obj = nhpp_c(Data1=Data_group1, Data2= Data_group2, 
                 Z1=Z_group1, Z2=Z_group2,
                 X1=X_group1, X2=X_group2)
    #system.time(nhpp_c(Data1=Data_group1, Data2= Data_group2, Z1=Z_group1, Z2=Z_group2))
  }else{
    obj = -Inf
  }
  return(obj)
}

# parallel search for the best variable and best split; 
# call tree.sub.fun1
tree.sub.fun2 = function(i, Data_id_sub2, X.B_sub2, Z.B_sub2, data.B_sub2,X.sample){
  x.target = X.B_sub2[Data_id_sub2,X.sample[i]]
  v.min = min(x.target,na.rm=TRUE)
  v.max = max(x.target,na.rm=TRUE)  
  x.delta = (v.max-v.min)/20
  x.seq = seq(v.min+3*x.delta,v.max-3*x.delta,x.delta)
  
  #obj = sapply(1:length(x.seq), tree.sub.fun1, thres.range=x.seq, 
  #            Data_on_Node=data.B_sub2[Data_id_sub2,],
  #            Z_on_Node=Z.B_sub2[Data_id_sub2],
  #            covariate=x.target)
  
  obj = parSapply(cl,1:length(x.seq), tree.sub.fun1, thres.range=x.seq, 
                  Data_on_Node=data.B_sub2[Data_id_sub2,],
                  Z_on_Node=Z.B_sub2[Data_id_sub2],
                  X_on_Node=X.B_sub2[Data_id_sub2,], 
                  covariate=x.target)
  
  return(obj)
}


