# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R.R
# Description: for each terminal node of a tree, perform regularized (Lasso) MLE to fit the NHPP model
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------

terminal_mle = function(terminal_mle_input, X.train, Z.train, data.train){
  n.terminal = length(terminal_mle_input)
  terminal.mle = array(0, dim=c(ncol(Z.train[[1]])+1, n.terminal))
  for (i in 1:n.terminal){
    terminal.data.id = terminal_mle_input[[i]][[2]]
    # pred:
    nhpp.t = nhpp.beta(Data = data.train[terminal.data.id, ], 
                       Z = Z.train[terminal.data.id], X=X.train[terminal.data.id, ]) 
    terminal.mle[,i] = nhpp.t
  }
  return(terminal.mle)
}
