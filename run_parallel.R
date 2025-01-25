library(huge)
library(mvtnorm)
library(alphastable)
library(fMultivar)
library(Matrix)
library(ROCR)
library(magrittr)
library(pcaPP)
library(SGL)
library(BB)
library(clue)

rm(list = ls())





Process <- function(ns, n, p, m, k, prec.type, g.type, mu = NULL,
                    run.oracle = T, run.raw = T, run.lasso = T, run.grouplasso = T,
                    run.Winsorized = T, run.kendall = T, run.cmc = T){
  source("generate_func.R")
  source("cmc_method.R")
  
  ### Generate data
  seed = 20240320 + 100*ns
  set.seed(seed)
  Dat = Generate_XZ(n, p, m, k, prec.type, mu, g.type)
  X = Dat$X; Z = Dat$Z; group.index = Dat$group.index
  Sigma.star = Dat$Sigma; Omega.star = Dat$Omega; edge.star = Dat$edge.mat
  rm(Dat)
  
  
  Res = NULL; methods = NULL
  
  ### Our method
  tic = proc.time()
  X.NS = Nonparanormal_trans_NormalScore(X)
  res.NS = SGL_BIC(X.NS, group.index)
  toc = proc.time() - tic
  time = toc[3]
  
  eva.NS = Precision.Evaluate(res.NS$Omega, Omega.star, res.NS$Edge, edge.star)
  Res = unname(eva.NS)
  methods = "Our"
  metrics.names = names(eva.NS)
  
  
  
  ### Oracle
  if(run.oracle){
    tic = proc.time()
    res.oracle = SGL_BIC(Z, group.index)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.oracle = Precision.Evaluate(res.oracle$Omega, Omega.star, res.oracle$Edge, edge.star)
    Res = cbind(Res, unname(eva.oracle))
    methods = c(methods, "Oracle")
  }
  
  
  ### Raw X
  if(run.raw){
    tic = proc.time()
    res.raw = SGL_BIC(X, group.index, use.corr = F)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.raw = Precision.Evaluate(res.raw$Omega, Omega.star, res.raw$Edge, edge.star)
    Res = cbind(Res, unname(eva.raw))
    methods = c(methods, "RawX")
  }
  
  ### Only Lasso
  if(run.lasso){
    tic = proc.time()
    X.NS = Nonparanormal_trans_NormalScore(X)
    res.lasso = SGL_BIC(X.NS, group.index, alpha=1)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.lasso = Precision.Evaluate(res.lasso$Omega, Omega.star, res.lasso$Edge, edge.star)
    Res = cbind(Res, unname(eva.lasso))
    methods = c(methods, "OnlyLasso")
  }
  
  
  ### Only GroupLasso
  if(run.grouplasso){
    tic = proc.time()
    X.NS = Nonparanormal_trans_NormalScore(X)
    res.grouplasso = SGL_BIC(X.NS, group.index, alpha=0)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.grouplasso = Precision.Evaluate(res.grouplasso$Omega, Omega.star, res.grouplasso$Edge, edge.star)
    Res = cbind(Res, unname(eva.grouplasso))
    methods = c(methods, "OnlyGroupLasso")
  }
  
  
  ### Winsorized transform
  if(run.Winsorized){
    tic = proc.time()
    X.Winsorized = Nonparanormal_trans_Winsorized(X)
    res.Winsorized = SGL_BIC(X.Winsorized, group.index)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.Winsorized = Precision.Evaluate(res.Winsorized$Omega, Omega.star, res.Winsorized$Edge, edge.star)
    Res = cbind(Res, unname(eva.Winsorized))
    methods = c(methods, "Winsorized")
  }
  
  
  ### Kendall's tau
  if(run.kendall){
    tic = proc.time()
    res.Kendall = SGL_BIC(X, group.index, cov_method = 2)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.Kendall = Precision.Evaluate(res.Kendall$Omega, Omega.star, res.Kendall$Edge, edge.star)
    Res = cbind(Res, unname(eva.Kendall))
    methods = c(methods, "Kendall")
  }
  
  
  ### CMC
  if(run.cmc){
    tic = proc.time()
    res.CMC = CMC_BIC(X, group.index)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.CMC = Precision.Evaluate(res.CMC$Omega, Omega.star, res.CMC$Edge, edge.star)
    Res = cbind(Res, unname(eva.CMC))
    methods = c(methods, "CMC")
  }
  
  
  
  rownames(Res) = metrics.names
  Res = as.data.frame(t(Res))
  Res$Method = methods
  Res$time = time
  
  return(Res)
}



#---------------------------------------------------------------------------------------------------------#
# General framework of parallel --------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------#

library(doParallel)
library(foreach)
# cl_num <- detectCores()*2/3
cl_num <- 2
cl <- makeCluster(cl_num,outfile="")
# cl <- makeCluster(cl_num,outfile="",type="FORK")
registerDoParallel(cl)
simu = 50

k = 4
n.seq = 300
p.seq = c(100, 150, 200)
m = 2    # size of each node
prec.seq = 1:4
g.seq = 2:3
setting = expand.grid(n.seq, p.seq, prec.seq, g.seq)


filename = paste("Run_smallp_sparse0.8_equal0.9_0115.RData",sep="")


AllRes = list()
AllResMean = list()
AllResSD = list()
for(tt in 1:nrow(setting)){
  n = setting[tt, 1]
  p = setting[tt, 2]
  prec.type = setting[tt, 3]
  g.type = setting[tt, 4]
  
  
  tic = proc.time()
  RES0 = foreach(ns=1:simu,.packages = c("Matrix", "mvtnorm", "fMultivar", "alphastable", "huge", "ROCR", "magrittr", "pcaPP", "SGL", "BB", "clue")) %dopar%
    Process(ns, n, p, m, k, prec.type, g.type, mu = NULL,
            run.oracle = F, run.raw = F, run.lasso = F, run.grouplasso = F, run.Winsorized = T, run.kendall = T, run.cmc = T)
  toc=proc.time()-tic
  
  RES0 = do.call(base::rbind,RES0)
  ResMean = aggregate(subset(RES0,select=-Method),by=list(RES0$Method),mean)
  ResSD = aggregate(subset(RES0,select=-Method),by=list(RES0$Method),sd)
  
  RES0$n = n
  RES0$p = p
  RES0$m = m
  RES0$PrecType = prec.type
  RES0$TransType = g.type
  
  ResMean$n = n
  ResMean$p = p
  ResMean$m = m
  ResMean$PrecType = prec.type
  ResMean$TransType = g.type
  
  ResSD$n = n
  ResSD$p = p
  ResSD$m = m
  ResSD$PrecType = prec.type
  ResSD$TransType = g.type
  
  
  AllRes[[tt]] = RES0
  AllResMean[[tt]] = ResMean
  AllResSD[[tt]] = ResSD
  
  total_seconds <- toc[3]
  hours <- total_seconds %/% 3600
  minutes <- (total_seconds %% 3600) %/% 60
  seconds <- round(total_seconds %% 60, digits = 3)
  
  print(paste(paste(tt,"Cost Time:",sep="-"), paste(hours, "hours,", minutes, "minutes, and", seconds, "seconds")))
}

stopCluster(cl)

AllResMean = do.call(base::rbind,AllResMean)
AllResSD = do.call(base::rbind,AllResSD)
save(AllRes, AllResMean, AllResSD, file=filename)
