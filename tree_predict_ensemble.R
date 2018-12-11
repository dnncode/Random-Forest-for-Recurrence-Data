# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R-MCF.R
# Description: perform the ensemble prediction
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------


tree_predict_ensemble = function(i, X.pred.ensemble, data.input, tree.trained.ensemble, data.pred.ensemble, B.hist){
  # X.pred: input X
  # tree.trained: trained tree stucture
  # X.pred.ensemble = X.OOB; tree.trained.ensemble= tree.list; data.input=data; data.pred.ensemble = data.OOB, B.hist=B.history
  
  tmp = tree.trained.ensemble[[i]]
  output = sapply(1:nrow(X.pred.ensemble), tree_predict_ensemble_sub,
                  X.pred=X.pred.ensemble, data.pred=data.pred.ensemble, 
                  tree.trained=tmp, data.train=data.input[B.hist[,i],])
  return(output)
  
}



tree_predict_ensemble_sub = function(i, X.pred, data.pred, tree.trained, data.train){
  # X.pred: input X
  # tree.trained: trained tree stucture
  # X.pred = X.OOB; tree.trained= terminal.list; data.train= data.B; data.pred = data.OOB
  
  n.terminal = length(tree.trained)
  output = sapply(1:n.terminal, tree_predict_sub, X.pred.0=X.pred[i,], tree.trained.0=tree.trained )
  terminal.select = tree.trained[[which(output==1)[1]]]
  
  terminal.struct = terminal.select[[1]]
  terminal.data.id = terminal.select[[2]]
  
  # pred:
  #mcf.pred = mcf(data.train[terminal.data.id, ])
  #n.failure = length(which(data.pred[i,]>0))-1
  #t.length = data.pred[i,n.failure+1]
  #t.length = 100
  #case = which(t.length - mcf.pred$time>0)
  #if (length(case)==0){
  #  MTBF = 0
  #}else{
  #  case = max(case,na.rm=TRUE) 
  #  MTBF = mcf.pred$mcf[case]
  #}
  
  mcf.pred = mcf(data.train[terminal.data.id, ])
  rate.pred = mcf.pred$mcf[nrow(mcf.pred)]/mcf.pred$time[nrow(mcf.pred)]

  return(rate.pred)
}


tree_predict_sub = function(i, X.pred.0, tree.trained.0 ){
  tree = data.frame(tree.trained.0[[i]][[1]])[-1,1:3]
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