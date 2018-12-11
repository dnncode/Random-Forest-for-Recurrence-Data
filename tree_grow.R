# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R-MCF.R
# Description: based on the results returned from "tree_split_nhpp_2.R", perform node split and update the tree structure
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# -

tree.grow = function(i.node, Data_input, X.B_input,Terminal.list_input){
  # Data_input=data.B; X.B_input = X.B; Terminal.list_input = terminal.list
  
  Initial.node_input = Terminal.list_input[[i.node]][[1]]
  Data_id_input = Terminal.list_input[[i.node]][[2]] 
  
  if (Initial.node_input[nrow(Initial.node_input),4]==0){
    
    if (ncol(X.B_input)<=2){
  	X.ratio = ncol(X.B_input)
    }else{
  	X.ratio = max( max(floor( ncol(X.B_input)/3 ), 1), 2 )
    }
    if (X.ratio>4){X.ratio=4}

    x.sample = sample(1:ncol(X.B_input), X.ratio, replace = FALSE)
    obj.list = lapply(1:length(x.sample), tree.sub.fun2, Data_id_sub2=Data_id_input,
                      X.B_sub2=X.B_input, data.B_sub2=Data_input, X.sample = x.sample)
    obj = do.call(cbind,obj.list)
    

    if (TRUE){
      x.id.opt = which(apply(obj,2,max,na.rm=TRUE)== max(apply(obj,2,max,na.rm=TRUE),na.rm=TRUE))
      x.opt = x.sample[x.id.opt]
      x.target = X.B_input[Data_id_input,
                           x.opt]
      v.min = min(x.target,na.rm=TRUE)
      v.max = max(x.target,na.rm=TRUE)  
      x.delta = (v.max-v.min)/20
      x.seq = seq(v.min+3*x.delta,v.max-3*x.delta,x.delta)
      thres.opt = mean( x.seq[which(obj[,x.id.opt]==max(apply(obj,2,max,na.rm=TRUE),na.rm=TRUE))],na.rm=TRUE )
      daughter.l.id =Data_id_input[ which(X.B_input[Data_id_input,x.opt]<=thres.opt) ]
      daughter.r.id =Data_id_input[ which(X.B_input[Data_id_input,x.opt]>thres.opt) ]
      
      # test if this node is terminal
      terminal.l = terminal.test(Data_input[daughter.l.id,])
      terminal.r = terminal.test(Data_input[daughter.r.id,])
      
      # create terminal node list
      terminal.sub.list = list()
      tmp = list()
      tmp[[1]] = rbind(Initial.node_input, c(x.opt,thres.opt,1,terminal.l))
      tmp[[2]] = daughter.l.id
      terminal.sub.list[[1]] = tmp
      tmp[[1]] = rbind(Initial.node_input, c(x.opt,thres.opt,2,terminal.r))
      tmp[[2]] = daughter.r.id
      terminal.sub.list[[2]] = tmp
      
    }else{
      terminal.sub.list = list()
      tmp = list()
      tmp[[1]] = Initial.node_input
      tmp[[1]][nrow(tmp[[1]]),4] = 1
      tmp[[2]] = Data_id_input
      terminal.sub.list[[1]] = tmp
    }
    
  }else{
    #terminal.sub.list = Terminal.list_input[[i.node]]
    terminal.sub.list = list()
    tmp = list()
    tmp[[1]] = Initial.node_input
    tmp[[2]] = Data_id_input
    terminal.sub.list[[1]] = tmp
  }
  
  return(terminal.sub.list)

}


