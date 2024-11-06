library(huge)
library(mvtnorm)
library(alphastable)
library(fMultivar)
library(Matrix)
require(ROCR)

#-----------------------------------------------------------------------------------------------------------
#  Generate subblock covariance matrix
#  Input:
#     p      ------ dimension
#  Output:
#    sigma   ------ covariance matrix
#-----------------------------------------------------------------------------------------------------------


AR_block = function(p, rho = 0.5){
  return(toeplitz(rho^(0:(p-1))))
}


banded_block <- function(p){
  theta <- matrix(0, p, p)
  for (i in 1: p){
    for (j in i: p){
      theta[i,j]<-max(0,(1-abs(i-j)/10))
      theta[j,i]<-theta[i,j]
    }
  }
  return(theta)
}


sparse_block = function(p, q, prob = 0.4){
  theta = matrix(0, p, q)
  for (i in 1:p) {
    for(j in 1:q){
      theta[i, j] = rbinom(1, 1, prob) * runif(1, 0.2, 0.6) * sign(runif(1, -1, 1))
    }
  }
  theta = theta - diag(diag(theta))
  return(theta)
}


equal_block = function(p, q, A = 0.2){
  theta = matrix(A, p, q)
  theta = theta - diag(diag(theta))
  return(theta)
}





#-----------------------------------------------------------------------------------------------------------
#  Generate full precision matrix
#  Input:
#     p      ------ full dimension
#     m      ------ size of each node
#  Output:
#   Sigma    ------ covariance matrix
#   Omega    ------ precision matrix
#-----------------------------------------------------------------------------------------------------------

Chain_model = function(p, m = floor(p/20), subblock = "sparse", sparse.prob = 0.4, equal.element = 0.2){
  n.groups = floor(p/m)
  group.size = c(rep(m, n.groups-1), m + p%%m)
  group.index = rep(1:n.groups, group.size)
  start.ind = c(1, 1+cumsum(group.size[-n.groups]))
  end.ind = cumsum(group.size)
  perm = sample.int(n.groups)
  
  tmp = matrix(0, n.groups, n.groups)
  theta = matrix(0, p, p)
  for(i in 1:(n.groups-1)){
    node1 = perm[i]
    node2 = perm[i+1]
    tmp[node1, node2] = 1
    tmp[node2, node1] = 1
    if(subblock == "sparse")
      theta[start.ind[node1]:end.ind[node1], start.ind[node2]:end.ind[node2]] = sparse_block(group.size[node1], group.size[node2], prob = sparse.prob)
    if(subblock == "equal")
      theta[start.ind[node1]:end.ind[node1], start.ind[node2]:end.ind[node2]] = matrix(equal.element, group.size[node1], group.size[node2])
  }
  theta = theta + t(theta)
  for(i in 1:n.groups){
    theta[start.ind[i]:end.ind[i], start.ind[i]:end.ind[i]] = AR_block(group.size[i], rho = 0.5)
  }
  
  gamma = 0.5 - min(eigen(theta, symmetric = T, only.values = T)$values)
  Omega = theta + gamma*diag(p)
  Sigma = solve(Omega)
  return(list(Sigma = Sigma, Omega = Omega, edge.mat = tmp, group.index = group.index))
}


Nearest_neighbor_model = function(p, m = floor(p/20), k = 4, subblock = "sparse", sparse.prob = 0.4, equal.element = 0.2){
  n.groups = floor(p/m)
  group.size = c(rep(m, n.groups-1), m + p%%m)
  group.index = rep(1:n.groups, group.size)
  start.ind = c(1, 1+cumsum(group.size[-n.groups]))
  end.ind = cumsum(group.size)
  
  tmp = matrix(0, n.groups, n.groups)
  unit.square = matrix(runif(2*n.groups, 0, 1), n.groups, 2)
  for(i in 1:n.groups){
    dis.i = apply(unit.square, 1, function(x){return(sum((x-unit.square[i,])^2))})
    neighbor.i = order(dis.i)[2:(k+1)]
    tmp[i, neighbor.i] = 1
  }
  tmp = tmp + t(tmp)
  tmp[tmp>0] = 1
  
  # Remove edges randomly
  for(i in 1:n.groups){
    neighbor.i = which(tmp[i,] > 0)
    if(length(neighbor.i) > k){
      remove.i = sample(neighbor.i, length(neighbor.i - k))
      tmp[i, remove.i] = 0
      tmp[remove.i, i] = 0
    }
  }
  edge.loc = which(tmp > 0, arr.ind = T)
  
  theta = matrix(0, p, p)
  for(i in 1:nrow(edge.loc)){
    node1 = edge.loc[i, 1]
    node2 = edge.loc[i, 2]
    if(node1 < node2){
      if(subblock == "sparse")
        theta[start.ind[node1]:end.ind[node1], start.ind[node2]:end.ind[node2]] = sparse_block(group.size[node1], group.size[node2], prob = sparse.prob)
      if(subblock == "equal")
        theta[start.ind[node1]:end.ind[node1], start.ind[node2]:end.ind[node2]] = matrix(equal.element, group.size[node1], group.size[node2])
    }
  }
  theta = theta + t(theta)
  for(i in 1:n.groups){
    theta[start.ind[i]:end.ind[i], start.ind[i]:end.ind[i]] = AR_block(group.size[i], rho = 0.5)
  }
  
  gamma = 0.5 - min(eigen(theta, symmetric = T, only.values = T)$values)
  Omega = theta + gamma*diag(p)
  Sigma = solve(Omega)
  return(list(Sigma = Sigma, Omega = Omega, edge.mat = tmp, group.index = group.index))
}




#-----------------------------------------------------------------------------------------------------------
#  Generate data matrix
#  Input:
#     n      ------ number of observations
#     p      ------ dimension size
#     m      ------ size of each node
#     prec.type  ------ type of covariance matrix: (1) Chain+sparse; (2) Chain+equal; (3) Neighbor+sparse; (4) Neighbor+equal
#     g.type ------ type of transform function: (1) Zou Hui JASA 2022; (2) Gaussian CDF; (3) Symmetric Power 
#  Output:
#     X      ------ n x p data matrix
#-----------------------------------------------------------------------------------------------------------
Generate_XZ = function(n, p, m = floor(p/20), k = 4, prec.type = 1, mu = NULL, g.type = 1, sparse.prob = 0.4, equal.element = 0.2,
                       mu.trans = 0.05, sigma.trans = 0.4, alpha.trans = 3){
  if(is.null(mu))
    mu = rep(0, p)
  
  Model = switch(prec.type,
                 Chain_model(p, m, subblock = "sparse", sparse.prob = sparse.prob),
                 Chain_model(p, m, subblock = "equal", equal.element = equal.element),
                 Nearest_neighbor_model(p, m, k, subblock = "sparse", sparse.prob = sparse.prob),
                 Nearest_neighbor_model(p, m, k, subblock = "equal", equal.element = equal.element))
  Z = matrix(rnorm(n*p),n,p) %*% chol(Model$Sigma) + matrix(mu, n, p)
  X = matrix(0, n, p)
  
  if(g.type == 1){
    for(j in 1:p){
      X[,j] = switch(j %% 5 + 1,
                     exp(Z[,j]),
                     Z[,j],
                     sign(Z[,j]) * sqrt(abs(Z[,j])),
                     pnorm(Z[,j]),
                     Z[,j]^3)
    }
  }
  
  if(g.type == 2){
    g0 = function(t){return(pnorm((t-mu.trans)/sigma.trans))}
    for(j in 1:p){
      mu.j = mu[j]
      sigma.j = Model$Sigma[j,j] %>% sqrt()
      # temp1 = function(t){return(g0(t) * dnorm((t-mu.j)/sigma.j))}
      # temp2 = integrate(temp1, -10, 10)$value
      # temp3 = function(y){return((g0(y) - temp2)^2 * dnorm((y-mu.j)/sigma.j))}
      # temp4 = integrate(temp3, -10, 10)$value %>% sqrt()
      # X[,j] = mu.j + sigma.j * (g0(Z[,j]) - temp2) / temp4
      MC.sample = g0(rnorm(1e5, mu.j, sigma.j))
      X[,j] = mu.j + sigma.j * (g0(Z[,j]) - mean(MC.sample)) / sd(MC.sample)
    }
  }
  
  if(g.type == 3){
    g0 = function(t){return(sign(t)*abs(t)^alpha.trans)}
    for(j in 1:p){
      mu.j = mu[j]
      sigma.j = Model$Sigma[j,j] %>% sqrt()
      # temp1 = function(t){return((g0(t-mu.j))^2 * dnorm((t-mu.j)/sigma.j))}
      # temp2 = integrate(temp1, -100, 100)$value %>% sqrt()
      # X[,j] = mu.j + sigma.j * g0(Z[,j] - mu.j) / temp2
      MC.sample = g0(rnorm(1e5, 0, sigma.j))
      X[,j] = mu.j + sigma.j * g0(Z[,j] - mu.j) / sd(MC.sample)
    }
  }
  
  return(list(X = X, Z = Z,
              Sigma = Model$Sigma, Omega = Model$Omega, 
              edge.mat = Model$edge.mat, group.index = Model$group.index))
}





#---------------------------------------------------------------------------------------------------------
#  Precision matrix evaluation
#  Input:
#     mat.est   ------ estimation of precision matrix.
#     mat.true  ------ true precision matrix.
#  Output:
#     
#---------------------------------------------------------------------------------------------------------
Precision.Evaluate = function(mat.est, mat.true, edge.est, edge.true){
  dmat = mat.est - mat.true
  d_M = norm(dmat, "M")
  d_S = norm(dmat, "2")
  d_F = norm(dmat, "F")
  
  pred <- prediction(abs(as.vector(mat.est)), abs(as.vector(mat.true)) > 1e-6)
  auc <- unlist(performance(pred, measure = "auc")@y.values)
  
  supp.est = which(abs(as.vector(edge.est)) > 1e-6)
  supp.true = which(abs(as.vector(edge.true)) > 1e-6)
  rec = length(intersect(supp.true, supp.est))/length(supp.true)
  prec = length(intersect(supp.true, supp.est))/max(length(supp.est),1)
  F1.score = ifelse(prec + rec == 0, 0, 2*prec*rec / (prec + rec))

  
  res = c(d_M, d_S, d_F, auc, prec, rec, F1.score)
  names(res) = c("Max", "Spectral", "Frobenius", "AUC", "Precision", "Recall", "F1-Score")
  return(res)
}

