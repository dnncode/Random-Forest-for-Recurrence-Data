# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R.R
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


tree_predict_C = function(i, X.train, X.pred, data.pred, tree.trained, data.train, Z.pred, Z.train, tree.mle){
  # X.pred: input X
  # tree.trained: trained tree stucture

  n.terminal = length(tree.trained)
  output = sapply(1:n.terminal, tree_predict_sub, X.pred.0=X.pred[i,], tree.trained.0=tree.trained )
  terminal.case = which(output==1)[1]
  terminal.select = tree.trained[[terminal.case]]
  

  terminal.struct = terminal.select[[1]]
  terminal.data.id = terminal.select[[2]]
  
  # actual
  n.failure = length(which(data.pred[i,]>0))-1
  t.length = data.pred[i,n.failure+1]
  rate = n.failure/t.length
  
  # pred:
  #nhpp.t = nhpp.beta(Data = data.train[terminal.data.id, ], Z = Z.train[terminal.data.id], X=X.train[terminal.data.id, ])
  nhpp.t = tree.mle[,terminal.case]
  #x.mat = t( matrix( rep(as.numeric(X.pred[i,]), nrow(Z.pred[[i]])), ncol=nrow(Z.pred[[i]])) )
  #tmp = cbind(array(1, dim=c(nrow(Z.pred[[i]]),1)), x.mat, Z.pred[[i]])
  tmp = cbind(array(1, dim=c(nrow(Z.pred[[i]]),1)), Z.pred[[i]])
  rate.pred = mean( exp(tmp %*% nhpp.t), na.rm=TRUE)

  # oob error
  error = c(rate.pred, rate)
  return(error)
}


tree_predict_sub = function(i, X.pred.0, tree.trained.0 ){
  tree = data.frame(   matrix(tree.trained.0[[i]][[1]][-1,1:3],ncol=3)    )
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














