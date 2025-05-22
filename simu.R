source('fn.R')    # load functions
library(fda)
library(fields)
library(clusteval)
library(TSclust)
library(scatterplot3d)
library(splines)
library(MASS)
library(matrixcalc) 
library(Matrix) 
library(gurobi)    # you need a valid gurobi license
library(numDeriv)
library(adagio)
library(mvtnorm)
library(pdist)
library(alphahull)
library(fdaoutlier)
library(sp)
library(matrixStats)
library(tclust)
library(mclust)
#######################################################################
########Simulation Experiment 1 and 2####################
#######################################################################

#######
n_iter = 100     #NUMBER OF ITERATION
nch = 20
ricrout = 0;rimsout = 0;riwdout = 0;ritvout = 0;ritrout = 0
ricrout2 = 0;rimsout2 = 0;riwdout2 = 0;ritvout2 = 0;ritrout2 = 0
set.seed(0)
tau = 0.5                     #RATIO OF CENTRAL CURVES
nc = 4                        #NUMBER OF CLUSTER
rate = 0.1                    #OUTLIER RATE, can be 0.1, 0.15, and 0.2
for(k in 1:n_iter){ 
  sim = generatesim2(type=1, T=200, rate=rate)    # reproduce experiment 2, comment this line
  #sim=generatesim(type=1, T=100, rate=rate)      # reproduce experiment 1, comment this line
  resultcr = robustcluster(sim,nc,'cr',tau=tau)$res
  resultms = robustcluster(sim,nc,'ms',tau=tau)$res
  resultwd = robustcluster(sim,nc,'ward',tau=tau)$res
  resulttr = tclust1(sim, nc, 1-tau)
  
  tru = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4)
  trutr = c(rep(1,150),rep(2,150),rep(3,150),rep(4,150))
  ind = which(resulttr == 0)
  trutr = trutr[-ind]
  resulttr = resulttr[-ind] 
  
  ricrout = ricrout + cluster.evaluation(resultcr,tru)/n_iter
  rimsout = rimsout + cluster.evaluation(resultms,tru)/n_iter
  riwdout = riwdout + cluster.evaluation(resultwd,tru)/n_iter
  ritrout = ritrout + cluster.evaluation(resulttr,trutr)/n_iter
  
  ricrout2 = ricrout2 + cluster_similarity(resultcr,tru)/n_iter
  rimsout2 = rimsout2 + cluster_similarity(resultms,tru)/n_iter
  riwdout2 = riwdout2 + cluster_similarity(resultwd,tru)/n_iter
  ritrout2 = ritrout2 + cluster_similarity(resulttr,trutr)/n_iter
  print(k) 
}
rimsout;ricrout;riwdout;ritrout;ritvout
rimsout2;ricrout2;riwdout2;ritrout2;ritvout2


