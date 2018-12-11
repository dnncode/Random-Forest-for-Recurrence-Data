# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Sub-routine for RF-R-MCF.R
# Description: calculate the distance between the fitted MCF between two daughter nodes
# Date: 04/20/2018
# ---------------------------------------------------------------
# Author: 
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# ---------------------------------------------------------------
# ---------------------------------------------------------------


mcf_c = function(Data1, Data2){
  

  # Approach 1: formal rank test (computationally intenstive)
  if (FALSE){  

  # to reduce the computation complexity; perform a random sample on Data1 and Data2
  if (nrow(Data1)>100){
    Data1 = Data1[sample(1:nrow(Data1), 100, replace = FALSE), ]
  }
  if (nrow(Data2)>100){
    Data2 = Data2[sample(1:nrow(Data2), 100, replace = FALSE), ]
  }


  n.c = 5 # number of test points
  
  #Data1 = data.B[1:100,]; Data2 = data.B[101:200,]
  mcf1 = mcf(Data1)
  mcf2 = mcf(Data2)
  
  #Get a number of common time points. 
  t1 = mcf1[,1]
  n.t1 = length(t1)
  t2 = mcf2[,1]
  n.t2 = length(t2)
  t.end = floor( min( t1[length(t1)], t2[length(t2)]) )
  t.step = floor( t.end / n.c )  # number of points to compare
  t.c = seq( floor(t.end/t.step), t.end, t.step )  # t for comparison
  #n.c = length(t.c)
  n.mcf1 = (ncol(mcf1)-6)/2
  n.mcf2 = (ncol(mcf2)-6)/2
  
  #Difference
  t1.c = t2.c = var.c1 = var.c2 = case1 = case2 = array()
  for (i in 1:n.c){
    case1[i] = max(  which( (t.c[i] >= t1) ) )
    t1.c[i] = mcf1[case1[i],"mcf"]
    var.c1[i] = mcf1[case1[i],"variance"]
    case2[i] = max(  which( (t.c[i] >= t2) ) )
    t2.c[i] = mcf2[case2[i],"mcf"]
    var.c2[i] = mcf2[case2[i],"variance"]
  }
  z = t1.c - t2.c 
  

  # get the covariance of Z
  # pre-computed cov.d
  tmp = lapply(c(1:n.c), mcf_c_sub1,n.c.sub1=n.c)
  CASE1 = do.call(rbind, tmp)
  tmp = lapply(c(1:nrow(CASE1)), mcf_c_sub2, mcf.in=mcf1, n.mcf=n.mcf1,case.in=CASE1)
  tmp = do.call(c, tmp)
  cov.d.precomputed = array(0/0, dim=c(n.c, n.c))
  for (i in 1:nrow(CASE1)){
    cov.d.precomputed[CASE1[i,1], CASE1[i,2]] = cov.d.precomputed[CASE1[i,2], CASE1[i,1]] = tmp[i]
  }
  cov.mat1 = array(0, dim=c(n.c, n.c))
  for (i in 1:(n.c-1)){
    for (j in (i+1):n.c){
      cov.mat1[i,j] = cov.mat1[j,i] = sum( cov.d.precomputed[1:i,1:j] )
    }
  }
  cov.mat1 = cov.mat1 + diag(var.c1)
  

  tmp = lapply(c(1:nrow(CASE1)), mcf_c_sub2, mcf.in=mcf2, n.mcf=n.mcf2,case.in=CASE1)
  tmp = do.call(c, tmp)
  cov.d.precomputed = array(0/0, dim=c(n.c, n.c))
  for (i in 1:nrow(CASE1)){
    cov.d.precomputed[CASE1[i,1], CASE1[i,2]] = cov.d.precomputed[CASE1[i,2], CASE1[i,1]] = tmp[i]
  }
  cov.mat2 = array(0, dim=c(n.c, n.c))
  for (i in 1:(n.c-1)){
    for (j in (i+1):n.c){
      cov.mat2[i,j] = cov.mat2[j,i] = sum( cov.d.precomputed[1:i,1:j] )
    }
  }
  cov.mat2 = cov.mat2 + diag(var.c2)
  

  cov.mat = cov.mat1 + cov.mat2 
  # eigen(cov.mat)$value  
  
  # chi-square test:
  z.chi = matrix(z[-1],nrow=1)
  cov.mat.chi = cov.mat[-1,-1]
  tmp = try(solve(cov.mat.chi), silent =TRUE)
  if (class(tmp)=="try-error"){
    chi.test = -Inf
  }else{
    chi.test = z.chi %*% solve(cov.mat.chi) %*% t(z.chi)
  }
  chi.statistic = pchisq(chi.test, df=n.c-1)
  #qchisq(0.9, df=n.c-1)
  #max.z = max(abs(z))
  return(chi.test)
  #return(chi.statistic) #chi.statistic = 1-p.value

  }


  # Approach 2: Average distance in failure rate
  if (TRUE){  
    # to reduce the computation complexity; perform a random sample on Data1 and Data2
    if (nrow(Data1)>100){
      Data1 = Data1[sample(1:nrow(Data1), 100, replace = FALSE), ]
    }
    if (nrow(Data2)>100){
      Data2 = Data2[sample(1:nrow(Data2), 100, replace = FALSE), ]
    }
    #n.c = 5 # number of test points 

    mcf1 = mcf(Data1)
    mcf2 = mcf(Data2)
 
    #t1 = mcf1[,1]
    t1 = mcf1[,1]
    n.t1 = length(t1)
    #t2 = mcf2[,1]
    t2 = mcf2[,1]
    n.t2 = length(t2)
    #t.end = floor( min( t1[length(t1)], t2[length(t2)]) )
    #t.step = floor( t.end / n.c )  # number of points to compare
    #t.c = seq( floor(t.end/t.step), t.end, t.step )  # t for comparison
    t.c = sort( unique(c(round(t1), round(t2))) )[-1]
    n.c = length(t.c) 

    #Difference
    t1.c = t2.c = case1 = case2 = array(0, dim=c(n.c,1))
    for (i in 1:n.c){
      case1[i] = max(  which( (t.c[i] >= t1) ) )
      t1.c[i] = mcf1[case1[i],"mcf"]
      case2[i] = max(  which( (t.c[i] >= t2) ) )
      t2.c[i] = mcf2[case2[i],"mcf"]
    }
    z = t1.c - t2.c
    L2.test = sqrt( sum( z^2 * c(t.c[1],diff(t.c)), na.rm=TRUE ) )/max(t.c, na.rm=TRUE)
    return(L2.test)
  }



  if (FALSE){
    test1 = mcf(Data1)
    test2 = mcf(Data2)
    plot(test2$time, test2$mcf,type="s")
    lines(test1$time, test1$mcf,col="red")
  }
  

}



# subroutines of mcf_c
mcf_c_sub1 = function(i,n.c.sub1){ # create a sequence for computing cov.d
  output = cbind( rep(i,(n.c.sub1-i+1)), c(i:n.c.sub1)   )
  return(output)
}

mcf_c_sub2 = function(i, mcf.in, n.mcf,case.in){ # get elements of precomputed cov.d based on CASE
  output = cov.d(delta=mcf.in[,2:(n.mcf+1)], 
                 delta.dot=mcf.in[,2*n.mcf+2], 
                 d=mcf.in[,(n.mcf+2):(2*n.mcf+1)], 
                 d.bar=mcf.in[,2*n.mcf+4], case.t=case.in[i,1], case.m=case.in[i,2])
  return(output)
}


cov.d = function(delta, delta.dot, d, d.bar, case.t, case.m){
  # delta = delta1; delta.dot = delta.dot1; d=d1; d.bar=d.bar1
  if (case.t<case.m){
    tmp = delta[case.m,]/delta.dot[case.m]*(d[case.m,]-d.bar[case.m])*d[case.t,]
    cov.dtk.dtm1 = sum(tmp)/delta.dot[case.t]
  }
  if (case.t>case.m){
    tmp = delta[case.t,]/delta.dot[case.t]*(d[case.t,]-d.bar[case.t])*d[case.m,]
    cov.dtk.dtm1 = sum(tmp)/delta.dot[case.m]
  }
  if (case.t==case.m){
    tmp = delta[case.t,]/delta.dot[case.t]*(d[case.t,]-d.bar[case.t])^2
    cov.dtk.dtm1 = sum(tmp)/delta.dot[case.t]
  }
  return(cov.dtk.dtm1)
}

