# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R.R
# Description: perform regularized (Lasso) MLE to fit the NHPP model (X and Z in the loglinear model)
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------

loglike.nhpp.2 = function(i, para, F.input, Z.input, X.input){
  f = F.input[i, ]
  z = Z.input[[i]]
  x = X.input[i,]
  #x.mat = t( matrix( rep(as.numeric(x), nrow(z)), ncol=nrow(z)) )
  #tmp = cbind(array(1, dim=c(nrow(z),1)), x.mat, z)
  tmp = cbind(array(1, dim=c(nrow(z),1)), z)
  lambda = exp(tmp %*% para)  
  tmp = max(  which(f>0))
  if (tmp==2){
    l = - sum(lambda[-1,1])
  }
  if (tmp>2){
    case = f[2:(tmp-1)]
    l = - sum(lambda[-1,1]) + sum( log(lambda[as.numeric(case),1])  )  
  }
  return(l)
} #loglike.nhpp.2(1, para=Par0, F.input=F.pass, Z.input=Z.pass, X.input=X.pass)

loglike_total_nhpp.2 = function(Par, F.pass, Z.pass, X.pass){
  l_i = sapply(1:nrow(F.pass), loglike.nhpp.2, para=Par, F.input=F.pass, Z.input=Z.pass, X.input=X.pass)
  l = -sum(l_i, na.rm=TRUE)
  return(l)
}

grr <- function(Par, F.pass, Z.pass, X.pass) { 
  grad = array(0, dim=c(length(Par),nrow(F.pass)))
  for (i in 1:nrow(F.pass)){
    f = F.pass[i, ]
    z = Z.pass[[i]]
    x = X.pass[i,]
    #x.mat = t( matrix( rep(as.numeric(x), nrow(z)), ncol=nrow(z)) )
    #tmp = cbind(array(1, dim=c(nrow(z),1)), x.mat, z)
    tmp = cbind(array(1, dim=c(nrow(z),1)), z)
    lambda = exp(tmp %*% Par)
    case = f[2:(max(  which(f>0))-1)]
    for (ii in 1:length(Par)){
      term1 = lambda * tmp[,ii] 
      term2 = tmp[as.numeric(case),ii] 
      grad[ii,i] = sum(term1[-1])-sum( term2 )  
    }
  }
  grad = rowSums(grad)
  return(grad)
}
# grr(Par=Par0,F.pass=Data, Z.pass=Z, X.pass=X)

nhpp.beta = function(Data, Z, X){   
  #Dim = ncol(Z[[1]])+1+ncol(X)
  Dim = ncol(Z[[1]])+1
  Par0 = array(0,dim=c(Dim,1))
  Par0[1] = log(0.01)
  mle = optim(Par0, loglike_total_nhpp.2, gr=grr, method="BFGS", F.pass=Data, Z.pass=Z, X.pass=X)
  output = mle$par 
  return(output)
}
# nhpp.beta(Data = data, Z = z.list)



