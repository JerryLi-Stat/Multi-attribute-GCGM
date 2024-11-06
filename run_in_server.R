library(huge)
library(mvtnorm)
library(alphastable)
library(fMultivar)
library(Matrix)
require(ROCR)
library(magrittr)

rm(list = ls())





Process <- function(ns, n, p, m, k, prec.type, g.type, mu = NULL,
                    run.oracle = T, run.raw = T, run.lasso = T, run.grouplasso = T){
  source("generate_func.R")
  source("method.R")
  
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
  res.NS = Sparse_Group_Graphical_Lasso_BIC(X.NS, group.index)
  toc = proc.time() - tic
  time = toc[3]

  eva.NS = Precision.Evaluate(res.NS$Omega, Omega.star, res.NS$edge, edge.star)
  Res = unname(eva.NS)
  methods = "NS_SGL"
  metrics.names = names(eva.NS)
  
  
  
  ### Oracle
  if(run.oracle){
    tic = proc.time()
    res.oracle = Sparse_Group_Graphical_Lasso_BIC(Z, group.index)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.oracle = Precision.Evaluate(res.oracle$Omega, Omega.star, res.oracle$edge, edge.star)
    Res = cbind(Res, unname(eva.oracle))
    methods = c(methods, "Oracle")
  }
  
  
  ### Raw X
  if(run.raw){
    tic = proc.time()
    res.raw = Sparse_Group_Graphical_Lasso_BIC(X, group.index, use.corr = F)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.raw = Precision.Evaluate(res.raw$Omega, Omega.star, res.raw$edge, edge.star)
    Res = cbind(Res, unname(eva.raw))
    methods = c(methods, "RawX")
  }

  ### Only Lasso
  if(run.lasso){
    tic = proc.time()
    X.NS = Nonparanormal_trans_NormalScore(X)
    res.lasso = Sparse_Group_Graphical_Lasso_BIC(X.NS, group.index, alpha=1)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.lasso = Precision.Evaluate(res.lasso$Omega, Omega.star, res.lasso$edge, edge.star)
    Res = cbind(Res, unname(eva.lasso))
    methods = c(methods, "OnlyLasso")
  }
  
  
  ### Only GroupLasso
  if(run.grouplasso){
    tic = proc.time()
    X.NS = Nonparanormal_trans_NormalScore(X)
    res.grouplasso = Sparse_Group_Graphical_Lasso_BIC(X.NS, group.index, alpha=0)
    toc = proc.time() - tic
    time = c(time, toc[3])
    
    eva.grouplasso = Precision.Evaluate(res.grouplasso$Omega, Omega.star, res.grouplasso$edge, edge.star)
    Res = cbind(Res, unname(eva.grouplasso))
    methods = c(methods, "OnlyGroupLasso")
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
cl_num <- 30
# cl <- makeCluster(cl_num,outfile="")
cl <- makeCluster(cl_num,outfile="",type="FORK")
registerDoParallel(cl)
simu = 100

k = 4
n.seq = c(500, 1000) 
p.seq = c(150, 300, 450)
prec.seq = 1:4
g.seq = 2:3
setting = expand.grid(n.seq, p.seq, prec.seq, g.seq)


filename = paste("RandomDesign_Change_n_p_0426.RData",sep="")


AllRes = list()
AllResMean = list()
AllResSD = list()
for(tt in 1:nrow(setting)){
  n = setting[tt, 1]
  p = setting[tt, 2]
  prec.type = setting[tt, 3]
  g.type = setting[tt, 4]
  m = 3
  
  tic = proc.time()
  RES0 = foreach(ns=1:simu,.packages = c("Matrix", "mvtnorm", "fMultivar", "alphastable", "huge", "ROCR", "magrittr")) %dopar%
    Process(ns, n, p, m, k, prec.type, g.type, mu = NULL,
            run.oracle = T, run.raw = T, run.lasso = F, run.grouplasso = T)
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




#-----------------------------------------------------------------------------#
# Transfer result matrix into csv file ---------------------------------------
#-----------------------------------------------------------------------------#
# library(openxlsx)
# rm(list = ls())
# n.seq = c(500, 1000) 
# p.seq = c(150, 300, 450)
# prec.seq = 1:4
# g.seq = 2:3
# sheet.setting = expand.grid(prec.seq, g.seq)
# sheets.id = 1:nrow(sheet.setting)
# number.of.sheets = length(sheets.id)
# prec.names = c("Chain+sparse",  "Chain+equal", "Neighbor+sparse", "Neighbor+equal")
# trans.names = c("Mixture in ZouHui", "Gaussian CDF", "Symmetric Power")
# 
# # # for(d in 1:number_of_sheets){
# # #   filename = paste0("RandomDesign_n100_dist", d, "_Change_p_sigma_1027.RData")
# # #   load(filename)
# # #   output[[d]] = Summary_sheet(AllResMean, AllResSD)
# # # }
# # 
# # # AllResMean[which(AllResMean$DistType == 1), ]
# # 
# load("~/Rcode/RandomDesign_Change_n_p_0418.RData")
# output = list()
# for(i in 1:number.of.sheets){
#   prec.id = sheet.setting[i, 1]
#   g.id = sheet.setting[i, 2]
#   sheet.rows = which((AllResMean$PrecType == prec.id) & (AllResMean$TransType == g.id))
#   sheet.mean = AllResMean[sheet.rows, ]
#   sheet.sd = AllResSD[sheet.rows, ]
#   temp = Summary_sheet(sheet.mean, sheet.sd)
#   temp = cbind(temp, prec.names[prec.id])
#   temp = cbind(temp, trans.names[g.id])
#   colnames(temp)[(ncol(temp)-1) : ncol(temp)] = c("PrecType", "TransType")
#   output[[i]] = temp
# }
# # names(output) = paste0("DistType", 1:number_of_sheets)
# write.xlsx(output, "SGL_0426.xlsx")
# # 
# # 
# # The columns of output matrix are c("n", "p", "our", "GL", "Oracle", "RawX", "index_name")
# Summary_sheet <- function(ResMean, ResSD){
#   methods = unique(ResMean$Group.1)
#   Evaluation_indexes = names(ResMean)[2:8]   # "Max", "Spectral", "Frobenius", "AUC", "Precision", "Recall", "F1-score"
#   number.of.settings = nrow(ResMean)/length(methods)
#   output = matrix(0, number.of.settings*length(Evaluation_indexes), length(methods)+3)
# 
#   combine_mean_and_sd = function(mean_vec, sd_vec){
#     mean_vec = sprintf("%.3f", mean_vec)
#     sd_vec = sprintf("%.3f", sd_vec)
#     L = length(mean_vec)
#     res = rep(0, L)
#     for (l in 1:L) {
#       res[l] = paste0(mean_vec[l], "(", sd_vec[l], ")")
#     }
#     return(res)
#   }
# 
#   pointer = 1
#   for (i in n.seq){ # loop for n
#     temp_mean = ResMean[ResMean$n==i,]
#     temp_sd = ResSD[ResSD$n==i,]
#     for (k in 2:8) { # loop for index
#       for (j in p.seq) { # loop for p
#         location = which(temp_mean$p==j)
#         output[pointer, 1] = i
#         output[pointer, 2] = j
#         output[pointer, 2+(1:length(methods))] =
#           combine_mean_and_sd(temp_mean[location, k], temp_sd[location, k])
#         output[pointer, 3+length(methods)] = names(ResMean)[k]
#         pointer = pointer + 1
#       }
#     }
#   }
# 
#   # sort_method = c(3, 2, 1)    # change the order of methods in output
#                                 # c("our", "rcec", "coat")
#   sort_method = 1:length(methods)
#   output[, 2+(1:length(methods))] = output[, 2+sort_method]
#   colnames(output) = c("n", "p", methods[sort_method], "index_name")
#   return(output)
# }
# 







#---------------------------------------------------------------------------#
# NCI-60 data -------------------------------------------------------------
#---------------------------------------------------------------------------#
# rm(list = ls())
# library(tidyverse)
# library(readxl)
# library(huge)
# library(cowplot)
# source("method.R")
# Protein <- read_xls("Protein__Antibody_Array_DTB_log2.xls", skip = 10) |>
#   tibble()
# Gene <- read_xls("RNA__Affy_HG_U95(A_E)_RMA.xls", skip = 10) |>
#   tibble()
# 
# # Check whether the cells are the same in these datasets
# all.equal(colnames(Protein)[7:66], colnames(Gene)[8:67])
# 
# # Select the common gene id's
# protein_id = Protein$`Entrez gene id e`
# protein_id_unique = unique(protein_id)
# gene_id = Gene$`Entrez gene id e`
# gene_id_unique = unique(gene_id)
# final_id = intersect(protein_id_unique, gene_id_unique)
# p = length(final_id)
# n = 60
# 
# # Filter the rows we need
# Protein_Dat <- Protein |>
#   filter(`Entrez gene id e` %in% final_id) |>
#   select(-c(1, 3:6)) 
# Gene_Dat <- Gene |>
#   filter(`Entrez gene id e` %in% final_id) |>
#   select(-c(1:2, 4:7))
# 
# # Construct our final dataset(60 x 180)
# # Note: we combine the rows with same id by their means
# GeneProtein <- matrix(0, n, p*2)
# for(j in 1:p){
#   GeneProtein[, 2*(j-1)+1] <- Gene_Dat |>
#     filter(`Entrez gene id e` == final_id[j]) |>
#     select(-1) |>
#     apply(2, mean)
#   GeneProtein[, 2*j] <- Protein_Dat |>
#     filter(`Entrez gene id e` == final_id[j]) |>
#     select(-1) |>
#     apply(2, mean)
# }
# group_id <- rep(1:p, each = 2)
# GeneProtein_NS = huge.npn(GeneProtein)
# 
# 
# 
# ### 5-fold CV to select tuning parameters
# Sparse_Group_Graphical_Lasso_CV <- function(X, group.index, nlambda = 10, lambda.min.ratio = 0.1, 
#                                             alpha0 = 0.1, alpha = seq(0.1, 0.9, by=0.1), rho0 = 2,
#                                             n.fold = 5, use.corr = T){
#   tic = proc.time()
#   n = nrow(X)
#   kp = ncol(X)
#   n.groups = length(unique(group.index))
#   group.size = unname(table(group.index))
#   S = cov(X) * (n-1)/n
#   if(use.corr){
#     sigmas = sqrt(diag(S))
#     S = diag(sigmas^(-1)) %*% S %*% diag(sigmas^(-1))
#     X = scale(X) * sqrt(n/(n-1))
#   }
#   
#   nalpha = length(alpha)
#   if(nalpha == 1)
#     alpha0 = alpha
#   fold.id <- rep(1:n.fold, length.out = n)[sample.int(n)]
#   lambda = AdaLambda_SparseGroup(X, group.index, S, nlambda, alpha0, rho0, lambda.min.ratio)
#   
#   # Select lambda with alpha=alpha0
#   ErrorMat.lambda = matrix(0, n.fold, nlambda)
#   Omega.list = list()
#   Omega.list[[1]] = matrix(0, kp, kp)
#   for(i in 1:n.fold){
#     X.train = X[fold.id != i, ]
#     X.valid = X[fold.id == i, ]
#     S.valid = cov(X.valid) #* (n-1)/n
#     for(l in 1:nlambda){
#       lambda.l = lambda[l]
#       W.init = Omega.list[[max(l-1, 1)]]
#       res.l = Sparse_Group_Graphical_Lasso_ADMM(X.train, group.index, lambda = lambda.l, alpha = alpha0, W.init = W.init, rho0 = rho0)
#       Omega.list[[l]] = res.l$Omega.hat
#       ErrorMat.lambda[i, l] = (n/2) * (sum(res.l$Omega.hat * S.valid) - determinant(res.l$Omega.hat, logarithm = TRUE)$modulus[1])
#     }
#   }
#   errors.lambda = apply(ErrorMat.lambda, 2, mean)
#   loc.lambda = which.min(errors.lambda)
#   lambda.fin = lambda[loc.lambda]
#   W.init = Omega.list[[loc.lambda]]
#   
#   # Select alpha with lambda = lam
#   if(nalpha == 1){
#     alpha.fin = alpha
#   }
#   else{
#     ErrorMat.alpha = matrix(0, n.fold, nalpha)
#     Omega.list = list()
#     Omega.list[[1]] = matrix(0, kp, kp)
#     for(i in 1:n.fold){
#       X.train = X[fold.id != i, ]
#       X.valid = X[fold.id == i, ]
#       S.valid = cov(X.valid) #* (n-1)/n
#       for(l in 1:nalpha){
#         alpha.l = alpha[l]
#         W.init = Omega.list[[max(l-1, 1)]]
#         res.l = Sparse_Group_Graphical_Lasso_ADMM(X.train, group.index, lambda = lambda.fin, alpha = alpha.l, W.init = W.init, rho0 = rho0)
#         Omega.list[[l]] = res.l$Omega.hat
#         ErrorMat.alpha[i, l] = (n/2) * (sum(res.l$Omega.hat * S.valid) - determinant(res.l$Omega.hat, logarithm = TRUE)$modulus[1])
#       }
#     }
#     errors.alpha = apply(ErrorMat.alpha, 2, mean)
#     loc.alpha = which.min(errors.alpha)
#     alpha.fin = alpha[loc.alpha]
#     W.init = Omega.list[[loc.alpha]]
#   }
#   
#   # Final computation
#   res.fin = Sparse_Group_Graphical_Lasso_ADMM(X, group.index, lambda = lambda.fin, alpha = alpha.fin, W.init = W.init, rho0 = rho0)
#   Omega.fin = res.fin$Omega.hat
#   if(use.corr){
#     Omega.fin = diag(sigmas^(-1)) %*% Omega.fin %*% diag(sigmas^(-1))
#   }
#   edge.fin = Find_edge(Omega.fin, group.index)
#   
#   toc = proc.time() - tic
#   return(list(Omega = Omega.fin, edge = edge.fin,
#               lambda = lambda, lambda.fin = lambda.fin,
#               alpha = alpha, alpha.fin = alpha.fin,
#               time = unname(toc[3])))
# }
# 
# 
# 
# SGL_StableSelection <- function(X, group.index, lambda, alpha, edge.full, B = 100){
#   X = as.matrix(X)
#   n = nrow(X)
#   kp = ncol(X)
#   n.groups = length(unique(group.index))
#   edge.percent = matrix(0, n.groups, n.groups)
#   reproduced.proportion <- 0
#   
#   for(t in 1:B){
#     X.t = X[sample.int(n, replace = T), ]
#     res.t = Sparse_Group_Graphical_Lasso_ADMM(X.t, group.index, lambda = lambda, alpha = alpha)
#     edge.percent = edge.percent + res.t$edge.hat
#     reproduced.proportion <- reproduced.proportion + sum(res.t$edge.hat + edge.full == 2) / sum(edge.full)
#     if(t %% 10 == 0){cat(t%/%10)}
#   }
#   edge.percent = edge.percent / B
#   reproduced.proportion <- reproduced.proportion / B
#   return(list(edge = edge.percent, prop = reproduced.proportion))
# }
# 
# 
# 
# res.NS.SGL = Sparse_Group_Graphical_Lasso_CV(GeneProtein_NS, group_id)
# res.NS.SGL.SS <- SGL_StableSelection(GeneProtein_NS, group_id, res.NS.SGL$lambda.fin, res.NS.SGL$alpha.fin, res.NS.SGL$edge)
# res.NS.SGL.SS$prop
# sum(res.NS.SGL.SS$edge >= 0.90)/2
# edge.NS.SGL.SS <- res.NS.SGL.SS$edge >= 0.90
# p.NS.SGL <- ggplot(melt(edge.NS.SGL.SS), aes(x=Var1, y=Var2, fill=value))+ 
#   geom_tile() +
#   scale_fill_manual(values = c("white", "black")) +
#   scale_y_reverse() +
#   labs(x = "", y = "", title = "NS-SGL") +
#   theme(legend.position = "none")
# 
# 
# res.NS.GL = Sparse_Group_Graphical_Lasso_CV(GeneProtein_NS, group_id, alpha = 0)
# res.NS.GL.SS <- SGL_StableSelection(GeneProtein_NS, group_id, res.NS.GL$lambda.fin, res.NS.GL$alpha.fin, res.NS.GL$edge)
# res.NS.GL.SS$prop
# sum(res.NS.GL.SS$edge >= 0.90)/2
# edge.NS.GL.SS <- res.NS.GL.SS$edge >= 0.90
# p.NS.GL <- ggplot(melt(edge.NS.GL.SS), aes(x=Var1, y=Var2, fill=value))+ 
#   geom_tile() +
#   scale_fill_manual(values = c("white", "black")) +
#   scale_y_reverse() +
#   labs(x = "", y = "", title = "NS-GL") +
#   theme(legend.position = "none")
# p.NS.GL
# 
# 
# res.Raw.SGL = Sparse_Group_Graphical_Lasso_CV(GeneProtein, group_id)
# res.Raw.SGL.SS <- SGL_StableSelection(GeneProtein, group_id, res.Raw.SGL$lambda.fin, res.Raw.SGL$alpha.fin, res.Raw.SGL$edge)
# res.Raw.SGL.SS$prop
# sum(res.Raw.SGL.SS$edge >= 0.90)/2
# edge.Raw.SGL.SS <- res.Raw.SGL.SS$edge >= 0.90
# p.Raw.SGL <- ggplot(melt(edge.Raw.SGL.SS), aes(x=Var1, y=Var2, fill=value))+ 
#   geom_tile() +
#   scale_fill_manual(values = c("white", "black")) +
#   scale_y_reverse() +
#   labs(x = "", y = "", title = "Raw-SGL") +
#   theme(legend.position = "none")
# p.Raw.SGL
# 
# plot_grid(p.Raw.SGL, p.NS.GL, p.NS.SGL, nrow = 1)
# 
# sum(edge.NS.SGL.SS + edge.Raw.SGL.SS == 2) / 2
# sum(edge.NS.GL.SS + edge.Raw.SGL.SS == 2) / 2

