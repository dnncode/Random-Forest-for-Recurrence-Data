# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R-MCF.R
# Description: compute the OOB prediction C-index
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------

tree_predict_C = function(i, X.pred, data.pred, tree.trained, data.train){
  # X.pred: input X
  # tree.trained: trained tree stucture
  # X.pred = X.OOB[1,]; tree.trained= terminal.list; data.train= data.B; data.pred = data.OOB[1,]
  
  n.terminal = length(tree.trained)
  output = sapply(1:n.terminal, tree_predict_sub, X.pred.0=X.pred[i,], tree.trained.0=tree.trained )
  terminal.select = tree.trained[[which(output==1)[1]]]
  
  terminal.struct = terminal.select[[1]]
  terminal.data.id = terminal.select[[2]]
  
  # actual
  n.failure = length(which(data.pred[i,]>0))-1
  t.length = data.pred[i,n.failure+1]
  rate = n.failure/t.length
  
  # pred:
  mcf.pred = mcf(data.train[terminal.data.id, ])
  case = which(t.length - mcf.pred$time>0)
  if (length(case)==0){case=1}
  case = max(case,na.rm=TRUE)
  rate.pred = mcf.pred$mcf[case]/mcf.pred$time[case]
  #plot(mcf.pred$time, mcf.pred$mcf)
  
  # oob error
  error = c(rate.pred, rate)
  return(error)
}


tree_predict_sub = function(i, X.pred.0, tree.trained.0 ){
  tree = data.frame(tree.trained.0[[i]][[1]][-1,1:3])
  tmp = 0
  for (j in 1:nrow(tree)){
    if (tree[j,3]==1){
      tmp = tmp + as.numeric( (X.pred.0[,tree[j,1]] - tree[j,2])<=0 )
    }else{
      tmp = tmp + as.numeric( (X.pred.0[,tree[j,1]] - tree[j,2])>0 )
    }
  }
  output = 0
  if (tmp == nrow(tree)){
    output = 1
  }
  return(output)
}

tree_predict_C_index = function(i, Pairs, Data.input){
  temp = Data.input[Pairs[i,],]
  judge = (temp[2,2] - temp[2,1] ) *  (temp[2,1] - temp[1,1] )
  Concordance = as.numeric((judge>=0))
  return(Concordance)
}














