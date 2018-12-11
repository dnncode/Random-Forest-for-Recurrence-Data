# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R.R
# Description: calculate the distance between to fitted NHPP between two daughter nodes
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------


nhpp_c = function(Data1, Data2, Z1, Z2, X1, X2){
  
  # Approach: Distance in two functions
 
    # to reduce the computation complexity; perform a random sample on Data1 and Data2
    if (nrow(Data1)>20){
      case = sample(1:nrow(Data1), 20, replace = FALSE)
      Data1 = Data1[case, ]
      Z1 = Z1[case]
      X1=X1[case,]
    }
    if (nrow(Data2)>20){
      case = sample(1:nrow(Data2), 20, replace = FALSE)
      Data2 = Data2[case, ]
      Z2 = Z2[case]
      X2=X2[case,]
    }
   
    nhpp1 = nhpp.beta(Data = Data1, Z = Z1, X=X1)
    nhpp2 = nhpp.beta(Data = Data2, Z = Z2, X=X2)
 
    t.c = min( max(Data1,na.rm=TRUE), max(Data2,na.rm=TRUE) )
   
    dist = f.dist(para1=nhpp1, para2=nhpp2)

    return(dist)
}


# Inner product of two rate functions
f.dist = function(para1, para2){
  #Example 1:
  dx = 0.1
  x = cbind(1,matrix( rep(seq(0,1,dx),length(para1)-1), ncol=length(para1)-1 ))
  Inner.D = sqrt( sum( ( (exp(x %*% para1) - exp(x %*% para2))^2 ) * dx ) )
  return(Inner.D)
}
#f.dist(para1=nhpp1, para2=nhpp2)




