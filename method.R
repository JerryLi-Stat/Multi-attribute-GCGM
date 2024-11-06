library(magrittr)
library(huge)

#----------------------------------------------------------------------------------------------------------------
#  Transform the data matrix X to the one following normal distribution 
#  with Winsorized estimator or normal score estimator.
#
#  Input:
#     X      ------ data matrix
#  Output:
#     X.trans   ------ the transformed data matrix
#----------------------------------------------------------------------------------------------------------------
Nonparanormal_trans_Winsorized = function(X){
  X = as.matrix(X)
  n = nrow(X); p = ncol(X)
  delta.n = 1/(4 * n^(1/4) * sqrt(pi*log(n)))
  X.mean = colMeans(X)
  X.sigma = apply(X, 2, sd) * sqrt((n-1)/n)
  X.trans = matrix(0, n, p)
  
  for(j in 1:p){
    X.j = X[,j] %>% as.vector()
    F.j.hat = rank(X.j)/n
    F.j.tilde = F.j.hat %>% pmax(., delta.n) %>% pmin(., 1-delta.n)
    X.trans[,j] = X.mean[j] + X.sigma[j] * qnorm(F.j.tilde)
  }
  # mu.trans = colMeans(X.trans)
  # Sigma.trans = cov(X.trans) * (n-1)/n
  
  return(X.trans)
  # return(list(X.trans = X.trans))
}


Nonparanormal_trans_NormalScore = function(X){
  X = as.matrix(X)
  n = nrow(X); p = ncol(X)
  X.mean = colMeans(X)
  X.sigma = apply(X, 2, sd) * sqrt((n-1)/n)
  X.trans = matrix(0, n, p)
  
  for(j in 1:p){
    X.j = X[,j] %>% as.vector()
    F.j.tilde = rank(X.j)/(n+1)
    X.trans[,j] = X.mean[j] + X.sigma[j] * qnorm(F.j.tilde)
  }
  # mu.trans = colMeans(X.trans)
  # Sigma.trans = cov(X.trans) * (n-1)/n
  
  return(X.trans)
}




#----------------------------------------------------------------------------------------------------------------
#  ADMM algorithm for sparse-group graphical Lasso
#
#  Input:
#     X             ------ n x p data matrix
#     group.index   ------ 1 x p vector, the index that indicates the group of each attribute
#     lambda, alpha ------ the tuning parameters
#     W.init        ------ Initial W. Used for warm start
#     Omega.ada     ------ The result precision matrix for the first time. Used for adaptive sparse group Lasso
#  Output:
#     Omega.hat     ------ p x p precision matrix
#     edge.hat      ------ n.groups x n.groups matrix, estimated edge-set
#----------------------------------------------------------------------------------------------------------------
Sparse_Group_Graphical_Lasso_ADMM = function(X, group.index, lambda, alpha, W.init = NULL, Omega.ada = NULL,
                                             rho0 = 2, mu = 10, tol.abs = 1e-4, tol.rel = 1e-4, iter.max = 1e3){
  X = as.matrix(X)
  n = nrow(X); p = ncol(X); n.groups = length(unique(group.index))
  S = cov(X) * (n-1)/n
  U.old = matrix(0, p, p)
  Omega.old = (diag(S))^(-1)
  if(is.null(W.init))
    W.init = matrix(0, p, p)
  if(is.null(Omega.ada))
    Omega.ada = matrix(1, p, p)
  W.old = W.init
  rho = rho0
  converged = F
  iter = 0
  
  while(!(converged) & (iter <= iter.max)){
    eigen.temp = eigen(S - rho * (W.old - U.old), symmetric = T)
    V.temp = eigen.temp$vectors
    d.temp = eigen.temp$values
    D.tilde = diag((-d.temp + sqrt(d.temp^2 + 4*rho)) / (2 * rho))
    Omega.new = V.temp %*% D.tilde %*% t(V.temp)
    
    W.new = matrix(0, p ,p)
    A.temp = Omega.new + U.old
    # off-diagnal subblocks of W, upper-triangle part
    for(i in 1:(n.groups-1)){
      for(j in (i+1):n.groups){
        attribute.i = which(group.index == i)
        attribute.j = which(group.index == j)
        B.temp = Soft_thre(A.temp[attribute.i, attribute.j], alpha*lambda/rho)
        W.new[attribute.i, attribute.j] = B.temp * pmax(0, 1 - lambda * (1-alpha) / (rho * norm(B.temp, "F")))
      }
    }
    W.new = W.new + t(W.new)    # complete the under-triangle part
    # diagnal subblocks of W
    for(i in 1:n.groups){
      attribute.i = which(group.index == i)
      A.ii = A.temp[attribute.i, attribute.i]
      OffDiagA.ii = A.ii - diag(diag(A.ii))
      W.new[attribute.i, attribute.i] = Soft_thre(OffDiagA.ii, alpha*lambda/rho) + diag(diag(A.ii))
    }
    
    U.new = U.old + Omega.new - W.new
    d.p = norm(Omega.new - W.new, "F")
    d.d = rho * norm(Omega.new - Omega.old, "F")
    tol.pri = p*tol.abs + tol.rel*max(norm(Omega.new, "F"), norm(W.new, "F"))
    tol.dual = p*tol.abs + tol.rel*norm(U.new, "F") * rho   # it should be ||U||_F * rho here
    if((d.p <= tol.pri) & (d.d <= tol.dual))
      converged = T
    if(d.p > mu*d.d){
      rho = 2*rho
      U.new = U.new / 2
    }
    if(d.d > mu*d.p){
      rho = rho/2
      U.new = 2 * U.new
    }
    
    W.old = W.new; Omega.old = Omega.new; U.old = U.new
    iter = iter + 1
  }
  
  Omega.hat = W.new
  edge.hat = Find_edge(Omega.hat, group.index)
  return(list(Omega.hat = Omega.hat,
              edge.hat = edge.hat))
}




Sparse_Group_Graphical_Lasso_BIC = function(X, group.index, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.1, 
                                            alpha0 = 0.1, alpha = seq(0.1, 0.9, by=0.1), 
                                            rho0 = 2, mu = 10, tol.abs = 1e-4, tol.rel = 1e-4, iter.max = 1e3,
                                            HowToTune = "OneByOne", use.corr = T){
  # Useful for test
  # lambda = NULL; nlambda = NULL; lambda.min.ratio = 0.1; alpha0 = 0.1; alpha = seq(0.1, 0.9, by=0.1)
  # rho0 = 2; mu = 10; tol.abs = 1e-4; tol.rel = 1e-4; iter.max = 1e3; use.corr = T
  
  tic = proc.time()
  X = as.matrix(X)
  n = nrow(X); p = ncol(X); n.groups = length(unique(group.index))
  group.size = unname(table(group.index))
  S = cov(X) * (n-1)/n
  if(use.corr){
    sigmas = sqrt(diag(S))
    S = diag(sigmas^(-1)) %*% S %*% diag(sigmas^(-1))
    X = scale(X) * sqrt(n/(n-1))
  }

  
  if(length(alpha) == 1)
    alpha0 = alpha
  if (!is.null(lambda)) 
    nlambda = length(lambda)
  if (is.null(lambda)) {
    if (is.null(nlambda)) 
      nlambda = 10
    if(HowToTune == "OneByOne")
      lambda = AdaLambda_SparseGroup(X, group.index, S, nlambda, alpha0, rho0, lambda.min.ratio)
    if(HowToTune == "ExpendGrid"){
      lambda = sapply(alpha, function(x){return(AdaLambda_SparseGroup(x, group.index, S, nlambda, x, rho0, lambda.min.ratio))})
    }
  }

  
  
  if(HowToTune == "OneByOne"){
    # Select lambda with alpha=alpha0
    Omega.list = list()
    Omega.list[[1]] = matrix(0, p, p)
    BIC.lambda = rep(0, nlambda)
    for(l in 1:nlambda){
      lambda.l = lambda[l]
      W.init = Omega.list[[max(l-1,1)]]   # warm start
      res.l = Sparse_Group_Graphical_Lasso_ADMM(X, group.index, lambda = lambda.l, alpha = alpha0, W.init = W.init,
                                                rho0 = rho0, mu = mu, tol.abs = tol.abs, tol.rel = tol.rel, iter.max = iter.max)
      Omega.list[[l]] = res.l$Omega.hat
      # BIC.lambda[l] = sum(diag(S %*% res.l$Omega.hat)) - log(det(res.l$Omega.hat)) +
      #                   sum(res.l$edge.hat * tcrossprod(group.size)) * log(n) / 2   # BIC in JMLR(2014)
      BIC.lambda[l] = sum(diag(S %*% res.l$Omega.hat)) - log(det(res.l$Omega.hat)) +
        sum(res.l$Omega.hat != 0) * log(n) / (2*n)                 # BIC in IEEE(2021)
    }
    BIC.min = min(BIC.lambda)
    lambda.fin = max(lambda[BIC.lambda <= BIC.min])    # there may be multiple lambdas that matches the minimum BIC, we choose the largest lambda
    
    # Select alpha with lambda=lambda.fin
    nalpha = length(alpha)
    Omega.list = list()
    Omega.list[[1]] = matrix(0, p, p)
    BIC.alpha = rep(0, nalpha)
    for(l in 1:nalpha){
      alpha.l = alpha[l]
      W.init = Omega.list[[max(l-1,1)]]   # warm start
      res.l = Sparse_Group_Graphical_Lasso_ADMM(X, group.index, lambda = lambda.fin, alpha = alpha.l, W.init = W.init,
                                                rho0 = rho0, mu = mu, tol.abs = tol.abs, tol.rel = tol.rel, iter.max = iter.max)
      Omega.list[[l]] = res.l$Omega.hat
      # BIC.alpha[l] = sum(diag(S %*% res.l$Omega.hat)) - log(det(res.l$Omega.hat)) +
      #                   sum(res.l$edge.hat * tcrossprod(group.size)) * log(n) / 2   # BIC in JMLR(2014)
      BIC.alpha[l] = sum(diag(S %*% res.l$Omega.hat)) - log(det(res.l$Omega.hat)) +
        sum(res.l$Omega.hat != 0) * log(n) / (2*n)                 # BIC in IEEE(2021)
    }
    BIC.min = min(BIC.alpha)
    alpha.fin = min(alpha[BIC.alpha <= BIC.min])    # there may be multiple alphas that matches the minimum BIC, we choose the smallest alpha
    id.fin = match(alpha.fin, alpha)
    Omega.fin = Omega.list[[id.fin]]
    if(use.corr){
      Omega.fin = diag(sigmas^(-1)) %*% Omega.fin %*% diag(sigmas^(-1))
    }
    edge.fin = Find_edge(Omega.fin, group.index)
  }
  

  
  
  
  if(HowToTune == "ExpendGrid"){
    nalpha = length(alpha)
    Omega.list = list()
    Omega.list[[1]] = matrix(0, p, p)
    BIC.lambda.alpha = matrix(0, nlambda, nalpha)


    for(i.l in 1:nlambda){
      for(i.a in 1:nalpha){
        alpha.i = alpha[i.a]
        lambda.i = lambda[i.l, i.a]
        W.init = Omega.list[[max((i.l-1)*nalpha + i.a - 1 - min(i.l-1, 1)*(nalpha-1), 1)]]
        res.i = Sparse_Group_Graphical_Lasso_ADMM(X, group.index, lambda = lambda.i, alpha = alpha.i, W.init = W.init,
                                                  rho0 = rho0, mu = mu, tol.abs = tol.abs, tol.rel = tol.rel, iter.max = iter.max)
        Omega.list[[(i.l-1)*nalpha + i.a]] = res.i$Omega.hat
        # BIC.lambda.alpha[i.l, i.a] = sum(diag(S %*% res.i$Omega.hat)) - log(det(res.i$Omega.hat)) +
        #                              sum(res.i$edge.hat * tcrossprod(group.size)) * log(n) / 2   # BIC in JMLR(2014)
        BIC.lambda.alpha[i.l, i.a] = sum(diag(S %*% res.i$Omega.hat)) - log(det(res.i$Omega.hat)) +
                                     sum(res.i$Omega.hat != 0) * log(n) / (2*n)                 # BIC in IEEE(2021)
      }
    }
    BIC.min = min(BIC.lambda.alpha)
    loc.min = which(BIC.lambda.alpha <= BIC.min, arr.ind = T)
    loc.lambda.min = min(loc.min[,1])
    loc.min = loc.min[loc.min[,1] == loc.lambda.min, ] %>% matrix(., ncol=2)
    loc.alpha.min = min(loc.min[,2])
    lambda.fin = lambda[loc.lambda.min]
    alpha.fin = alpha[loc.alpha.min]
    Omega.fin = Omega.list[[(loc.lambda.min-1)*nalpha + loc.alpha.min]]
    if(use.corr){
      Omega.fin = diag(sigmas^(-1)) %*% Omega.fin %*% diag(sigmas^(-1))
    }
    edge.fin = Find_edge(Omega.fin, group.index)
  }
  toc = proc.time() - tic
  

  
  return(list(Omega = Omega.fin, edge = edge.fin,
              lambda = lambda, lambda.fin = lambda.fin,
              alpha = alpha, alpha.fin = alpha.fin,
              time = unname(toc[3])))
}





### Compute the smallest lambda that gives a no-edge model
AdaLambda_SparseGroup = function(X, group.index, S, nlambda, alpha0 = 0.1, rho0 = 2, 
                                 lambda.min.ratio = 0.1, tau = 0.002){
  n.groups = length(unique(group.index))
  alpha0 = min(max(alpha0, 0.1), 0.9)  # we truncate alpha0 into [0.1, 0.9]
  eigen.S = eigen(S, symmetric = T)
  V.S = eigen.S$vectors
  d.S = eigen.S$values
  Omega.first.iter = V.S %*% diag((-d.S + sqrt(d.S^2 + 4*rho0)) / (2 * rho0)) %*% t(V.S)
  lambda.extrem = max(abs(Omega.first.iter - diag(diag(Omega.first.iter)))) * rho0 / alpha0
  lambda.max = lambda.extrem / 2
  lambda.min = lambda.min.ratio * lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  return(lambda)
}


Find_edge = function(Omega, group.index){
  n.groups = length(unique(group.index))
  edge.fin = matrix(0, n.groups, n.groups)
  for(i in 1:(n.groups-1)){
    for(j in (i+1):n.groups){
      attribute.i = which(group.index == i)
      attribute.j = which(group.index == j)
      edge.fin[i,j] = (norm(Omega[attribute.i, attribute.j], "F") > 0)
      edge.fin[j,i] = edge.fin[i,j]
    }
  }
  return(edge.fin)
}





#----------------------------------------------------------------------------------------------------------------
#  JMLR(2014)
#
#  Input:
#     X             ------ p x p data matrix
#     group.index   ------ 1 x p vector, the index that indicates the group of each attribute
#     lambda        ------ the tuning parameters
#     Omega.init    ------ Initial W. Used for warm start
#  Output:
#     Omega         ------ p x p precision matrix
#     edge          ------ n.groups x n.groups matrix, estimated edge-set
#----------------------------------------------------------------------------------------------------------------
Multi_Attribute_Graph = function(X, group.index, lambda, Omega.init = NULL, Sigma.init = NULL, tol = 1e-3, iter.max = 1e3){
  X = as.matrix(X)
  n = nrow(X); p = ncol(X); n.groups = length(unique(group.index))
  S = cov(X) * (n-1)/n
  if(is.null(Omega.init))
    Omega.init = diag(diag(S)^(-1))
  Omega = Omega.old = Omega.init
  if(is.null(Sigma.init))
    Sigma.init = diag(diag(S))
  Sigma = Sigma.init
  converged = F
  iter = 0
  
  while(!(converged) & (iter <= iter.max)){
    for(i in 1:n.groups){
      a.i = which(group.index == i)
      t = 1
      positive.define = F
      while(!(positive.define)){
        for(j in 1:n.groups){
          a.j = which(group.index == j)
          temp = Omega.old[a.i, a.j] + t*(Sigma[a.i, a.j] - S[a.i, a.j])
          if(norm(temp, "F") > t*lambda){
            Omega[a.i, a.j] = (1 - t*lambda/norm(temp, "F")) * temp
          }else{Omega[a.i, a.j] = 0}
        }
        Omega[-a.i, a.i] = t(Omega[a.i, -a.i])
        if(min(eigen(Omega, symmetric = T, only.values = T)$values) > 0)
          positive.define = T
        t = t/2
      } # End of while positive
      # Omega.old.minusi.inverse = Sigma[-a.i, -a.i] - Sigma[-a.i, a.i] %*% solve(Sigma[a.i, a.i]) %*% Sigma[a.i, -a.i]
      # Sigma[a.i, a.i] = solve(Omega[a.i, a.i] - Omega[a.i, -a.i] %*% Omega.old.minusi.inverse %*% Omega[-a.i, a.i])
      # Sigma[-a.i, -a.i] = Omega.old.minusi.inverse + Omega.old.minusi.inverse %*% Omega[-a.i, a.i] %*% Sigma[a.i, a.i] %*% Omega[a.i, -a.i] %*% Omega.old.minusi.inverse
      # Sigma[a.i, -a.i] = -Omega[a.i, a.i] %*% Omega[a.i, -a.i] %*% Sigma[-a.i, -a.i]
      # Sigma[-a.i, a.i] = t(Sigma[a.i, -a.i])
      # print(all.equal(Omega%*%Sigma, diag(p)))
      Sigma = solve(Omega)
    }# End of for i
    if(abs(sum(diag(S %*% Omega)) - log(det(Omega)) + lambda*Group_matrix_norm(Omega, group.index) - p - log(det(Sigma))) < tol)
      converged = T
    iter = iter + 1
    Omega.old = Omega
  }
  edge.hat = Find_edge(Omega, group.index)
  return(list(Omega = Omega, 
              edge = edge.hat,
              Sigma = Sigma))
}



Multi_Attribute_Graph_BIC = function(X, group.index, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.1, tol = 1e-3, iter.max = 1e3){
  ### For test
  # lambda = NULL; nlambda = NULL; lambda.min.ratio = 0.1; tol = 1e-3; iter.max = 1e3
  
  tic = proc.time()
  X = as.matrix(X)
  n = nrow(X); p = ncol(X); n.groups = length(unique(group.index))
  group.size = unname(table(group.index))
  S = cov(X) * (n-1)/n
  if (!is.null(lambda)) 
    nlambda = length(lambda)
  if (is.null(lambda)) {
    if (is.null(nlambda)) 
      nlambda = 10
    lambda = AdaLambda_Multi(X, group.index, nlambda, lambda.min.ratio)
  }
  
  
  Omega.list = list()
  Omega.list[[1]] = diag(diag(S)^(-1))
  Sigma.list = list()
  Sigma.list[[1]] = diag(diag(S))
  BIC.lambda = rep(0, nlambda)
  for(l in 1:nlambda){
    lambda.l = lambda[l]
    # Omega.init = Omega.list[[max(l-1,1)]]   # warm start
    # Sigma.init = Sigma.list[[max(l-1,1)]]
    Omega.init = diag(diag(S)^(-1))
    Sigma.init = diag(diag(S))
    res.l = Multi_Attribute_Graph(X, group.index, lambda = lambda.l, Omega.init = Omega.init, Sigma.init = Sigma.init, tol = tol, iter.max = iter.max)
    Omega.list[[l]] = res.l$Omega
    Sigma.list[[l]] = res.l$Sigma
    BIC.lambda[l] = sum(diag(S %*% res.l$Omega)) - log(det(res.l$Omega)) +
                      sum(res.l$edge * tcrossprod(group.size)) * log(n) / 2   # BIC in JMLR(2014)
    # BIC.lambda[l] = sum(diag(S %*% res.l$Omega)) - log(det(res.l$Omega)) +
    #   sum(res.l$Omega.hat != 0) * log(n) / (2*n)                 # BIC in IEEE(2021)
  }
  BIC.min = min(BIC.lambda)
  lambda.fin = max(lambda[BIC.lambda <= BIC.min])    # there may be multiple lambdas that matches the minimum BIC, we choose the largest lambda
  id.fin = match(lambda.fin, lambda)
  Omega.fin = Omega.list[[id.fin]]
  edge.fin = Find_edge(Omega.fin, group.index)
  
  toc = proc.time() - tic
  return(list(Omega = Omega.fin, edge = edge.fin,
              lambda = lambda, lambda.fin = lambda.fin,
              time = unname(toc[3])))
}


### Compute the smallest lambda that gives a no-edge model
AdaLambda_Multi = function(X, group.index, nlambda, lambda.min.ratio = 0.1){
  n.groups = length(unique(group.index))
  # Omega.init = diag(diag(S))
  # Sigma.init = diag(diag(S)^(-1))
  # temp = Omega.init + Sigma.init - S
  # lambda.max = 0
  # for(i in 1:n.groups){
  #   a.i = which(group.index == i)
  #   for(j in i:n.groups){
  #     a.j = which(group.index == j)
  #     lambda.max = max(lambda.max, norm(temp[a.i, a.j], "F"))
  #   }
  # }
  
  lambda.0 = 1
  finish = F
  temp = Multi_Attribute_Graph(X, group.index, lambda = lambda.0, iter.max = 1)
  if(sum(temp$edge) == 0){
    while(!(finish)){
      lambda.1 = lambda.0 * 0.8
      temp = Multi_Attribute_Graph(X, group.index, lambda = lambda.1, iter.max = 1)
      if(sum(temp$edge) > 0)
        finish = T
      lambda.0 = lambda.1
    }
    lambda.max = lambda.1 * 1.25
  }
  if(sum(temp$edge) > 0){
    while(!(finish)){
      lambda.1 = lambda.0 * 1.25
      temp = Multi_Attribute_Graph(X, group.index, lambda = lambda.1, iter.max = 1)
      if(sum(temp$edge) == 0)
        finish = T
      lambda.0 = lambda.1
    }
    lambda.max = lambda.1
  }

  lambda.min = lambda.min.ratio * lambda.max
  lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  return(lambda)
}


Group_matrix_norm <- function(mat, group.index, type = "F"){
  n.groups = length(unique(group.index))
  res = 0
  for(i in 1:n.groups){
    a.i = which(group.index == i)
    for(j in 1:n.groups){
      a.j = which(group.index == j)
      res = res + norm(mat[a.i, a.j], type = type)
    }
  }
  res
}

Soft_thre <- function(z, lambda){
  return(sign(z)*pmax(0, abs(z)-lambda))
}
