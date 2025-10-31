############################################################
# Title: SCPCA+PPCA
#
# Description: SCPCA+PPCA code for analyzing microbial community differences
#
############################################################

library(PoissonPCA)
library(abind)

# This is the CPCA function. X is a list of datasets,
# n_g is a vector with the number of obs in each dataset, k is the usual k (leave at 0 and it'll
# set it to p), and iter is the # of optimization iterations
steppy=function (X, n_g, k = 0, iter = 30, ...) 
{
  p <- dim(X)[1]
  mcas <- dim(X)[3]
  if (k == 0) {
    k <- p
  }
  # iter <- 15
  n <- n_g/sum(n_g)
  D <- array(0, dim = c(p, mcas))
  CPC <- array(0, dim = c(p, p))
  Qw <- diag(1, p)
  s <- array(0, dim = c(p, p))
  for (m in 1:mcas) {
    s <- s + n[m] * X[, , m]
  }
  res <- eigen(s)
  q0 <- res$vectors
  d0 <- diag(res$values, p)
  if (d0[1, 1] < d0[p, p]) {
    q0 <- q0[, ncol(q0):1]
  }
  for (ncomp in 1:k) {
    q <- q0[, ncomp]
    d <- array(0, dim = c(1, mcas))
    for (m in 1:mcas) {
      d[, m] <- t(q) %*% X[, , m] %*% q
    }
    for (i in 1:iter) {
      s <- array(0, dim = c(p, p))
      for (m in 1:mcas) {
        s <- s + n_g[m] * X[, , m]/d[, m]
      }
      w <- s %*% q
      if (ncomp != 1) {
        w <- Qw %*% w
      }
      q <- w/as.numeric(sqrt((t(w) %*% w)))
      for (m in 1:mcas) {
        d[, m] <- t(q) %*% X[, , m] %*% q
      }
    }
    D[ncomp, ] <- d
    CPC[, ncomp] <- q
    Qw <- Qw - q %*% t(q)
  }
  out <- list(D = D[1:ncomp, ], CPC = CPC[, 1:ncomp], ncomp = ncomp)
  return(out)
}

# This is a modified version of one of the inner functions of PPCA, which calls Toby's C code
# to generate the scores
loglambdascores=function (X, V, d, k, mu) 
{
  p <- length(X)
  dstar <- d
  dstar[dstar < 0.01] <- 0.01
  dstar[seq_len(k)] <- 0
  M <- V %*% diag(dstar) %*% t(V)
  answer <- vector("double", k + p + 2)
  returnval <- .C("get_scores_log", answer, as.integer(p), 
                  as.integer(X), as.double(as.vector(V)), as.double(dstar), 
                  as.integer(k), as.double(mu), as.double(as.vector(M)))
  answer <- returnval[[1]]
  return(list(scores = answer[seq_len(k)], means = answer[seq_len(p) + 
                                                            k], convergence = answer[p + k + 1], accuracy = answer[p + 
                                                                                                                     k + 2]))
}

## eigvecs is one p x p mat, eigvals is mat with p rows and S columns, ind indicates which column of eigvals (s from 1:S)
# future molly: this is a modified version of another inner function of PPCA (which calls the function above), and it only works 
# if you used a sequencing depth correction with PPCA (some lines need to be slightly diff otherwise).
# i don't remember exactly what exactly was modified about this compared to the version in the PPCA pkg
computescorez_seq=function(dat,eigvecs,eigvals,ind){
  k=dim(dat)[2]-2
  Pscores=matrix(0,dim(dat)[1],k)
  Pmeans=dat
  transformation=makelogtransformation(3,6)
  gX=transformation$g(as.vector(as.matrix(dat)))
  dim(gX)=dim(dat)
  mu=colMeans(gX)
  p <- dim(dat)[2]
  V=as.matrix(eigvecs)
  proj <- t(V) %*% rep(1/sqrt(p), p)
  index <- which(abs(proj) == max(abs(proj)))
  V <- V - rep(1/sqrt(p), p) %*% t(proj)
  V <- cbind(rep(1/sqrt(p), p), V[, -index])
  for (i in seq_len(p - 1)[-1]) {
    V[, i] <- V[, i]/sqrt(sum(V[, i]^2))
    V[, (i + 1):p] <- V[, (i + 1):p] - V[, i] %*% 
      t(V[, i]) %*% V[, (i + 1):p]
  }
  V[, p] <- V[, p]/sqrt(sum(V[, p]^2))
  D <- c(1, eigvals[,ind][-index])
  for(i in seq_len(dim(dat)[1])) {
    gsl=loglambdascores(X=dat[i,],V=V,d=D,k=k+1,mu=mu)
    Pscores[i,]=gsl$scores[seq_len(k)+1]
    Pmeans[i,]=gsl$means
  }
  
  return(Pscores)
}

# This is the wrapper that takes the data and runs PPCA, CPCA, and gets the scores
poissoncpc=function(data){
  corn=colnames(data[[1]])
  for(i in 2:length(data)) {
    corn=intersect(colnames(data[[i]]),corn)
  }
  corn=corn[-length(corn)]
  newdata=lapply(data,function(x){subset(x,select=corn)})
  p=length(corn)
  obj=list()
  var=list()
  for(i in seq_len(length(newdata))){
    obj[[i]]=Poisson_Corrected_PCA(as.matrix(newdata[[i]]), k=p-2, transformation = "log",seqdepth="compositional")
    print(paste("Done PoissonPCA for dataset",i,"of",length(newdata)))
    var[[i]]=obj[[i]]$variance
  }
  cpcprep=abind(var,along=3)
  n_g=sapply(newdata,nrow)
  commonpcs=steppy(X=cpcprep,n_g=n_g)
  
  stepscores=list()
  for(i in seq_len(length(newdata))){
    ind=i
    stepscores[[i]]=computescorez_seq(dat=as.matrix(newdata[[i]]),eigvecs=commonpcs$CPC,eigvals=commonpcs$D,ind)
    print(paste("Done computing scores for dataset",i,"of",length(newdata)))
  }
  return(list(scores=stepscores,data=newdata))
}
