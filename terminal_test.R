# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R.R
# Description: tree node split termination test
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------

terminal.test = function(test.data){
  
  # rule 1: at least d systems with failures. 
  d = 20 # the minimum number of systems with failures in a terminal node
  d.count = 0
  for (i in 1:nrow(test.data)){
    if ( length(  which(test.data[i,]==-1) ) < (ncol(test.data) - 2) ){
      d.count = d.count + 1
    }
  }
  output = 0
  if (d.count<d){
    output = 1
  }
  output
}
