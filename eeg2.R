source('GCVsmoothing.R')
library(fda)
library(fields)
library(clusteval)
library(TSclust)
library(scatterplot3d)
library(splines)
library(MASS)
library(matrixcalc) 
library(Matrix) 
library(gurobi)
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
library(factoextra)


load('pgshort2.rdata')
rescr_ari1 = 0;resms_ari1 = 0;rescr_sim1 = 0;resms_sim1 = 0
rescr_ari2 = 0;resms_ari2 = 0;rescr_sim2 = 0;resms_sim2 = 0
rescr_ari3 = 0;resms_ari3 = 0;rescr_sim3 = 0;resms_sim3 = 0
rescr_ari4 = 0;resms_ari4 = 0;rescr_sim4 = 0;resms_sim4 = 0
#DCBEA
pgshort2 = pgshort
k=10

for(s in 1:k){
  set.seed(s)
  pgshort2[1:100,] = pgshort2[sample(1:100,100), ]
  pgshort2[101:200,] = pgshort2[sample(101:200,100), ]
  pgshort2[201:300,] = pgshort2[sample(201:300,100), ]
  pgshort2[301:400,] = pgshort2[sample(301:400,100), ]
  pgshort2[401:500,] = pgshort2[sample(401:500,100), ]
  
  data = array(0,c(20,250,25))
  for(i in 1:25){data[,,i] = pgshort2[((i-1)*20+1):(i*20),1:250]}
  data=data[1:15,,]
  
  
  task1 = data[,,c(21:25,16:20)]  #AE
  cr = robustcluster(task1,2,'cr',0.7)
  ms = robustcluster(task1,2,'ms',0.7)
  tru = c(rep(1,5),rep(2,5))
  rescr_ari1 = rescr_ari1 + cluster_similarity(cr$res,tru)/k
  resms_ari1 = resms_ari1 + cluster_similarity(ms$res,tru)/k
  rescr_sim1 = rescr_sim1 + cluster.evaluation(cr$res,tru)/k
  resms_sim1 = resms_sim1 + cluster.evaluation(ms$res,tru)/k #AE
  
  task2 = data[,,c(21:25,1:5)]  #AD
  cr = robustcluster(task2,2,'cr',0.7)
  ms = robustcluster(task2,2,'ms',0.7)
  tru = c(rep(1,5),rep(2,5))
  cluster_similarity(cr$res,tru);cluster_similarity(ms$res,tru)
  cluster.evaluation(cr$res,tru);cluster.evaluation(ms$res,tru)
  rescr_ari2 = rescr_ari2 + cluster_similarity(cr$res,tru)/k
  resms_ari2 = resms_ari2 + cluster_similarity(ms$res,tru)/k
  rescr_sim2 = rescr_sim2 + cluster.evaluation(cr$res,tru)/k
  resms_sim2 = resms_sim2 + cluster.evaluation(ms$res,tru)/k 
  
  task3 = data[,,c(1:5,16:20)]  #DE
  cr = robustcluster(task3,2,'cr',0.7)
  ms = robustcluster(task3,2,'ms',0.7)
  tru = c(rep(1,5),rep(2,5))
  cluster_similarity(cr$res,tru);cluster_similarity(ms$res,tru)
  cluster.evaluation(cr$res,tru);cluster.evaluation(ms$res,tru)
  rescr_ari3 = rescr_ari3 + cluster_similarity(cr$res,tru)/k
  resms_ari3 = resms_ari3 + cluster_similarity(ms$res,tru)/k
  rescr_sim3 = rescr_sim3 + cluster.evaluation(cr$res,tru)/k
  resms_sim3 = resms_sim3 + cluster.evaluation(ms$res,tru)/k 
  
  task4 = data[,,c(16:20,1:5,21:25)] #ADE
  cr = robustcluster(task4,3,'cr',0.7)
  ms = robustcluster(task4,3,'ms',0.7)
  tru = c(rep(1,5),rep(2,5),rep(3,5))
  cluster_similarity(cr$res,tru);cluster_similarity(ms$res,tru)
  cluster.evaluation(cr$res,tru);cluster.evaluation(ms$res,tru)
  rescr_ari4 = rescr_ari4 + cluster_similarity(cr$res,tru)/k
  resms_ari4 = resms_ari4 + cluster_similarity(ms$res,tru)/k
  rescr_sim4 = rescr_sim4 + cluster.evaluation(cr$res,tru)/k
  resms_sim4 = resms_sim4 + cluster.evaluation(ms$res,tru)/k 
  print('============================')
  print(s)}

c(resms_ari1,rescr_ari1)
c(resms_sim1,rescr_sim1) 

c(resms_ari3,rescr_ari3)
c(resms_sim3,rescr_sim3) 

c(resms_ari2,rescr_ari2) 
c(resms_sim2,rescr_sim2) 

c(resms_ari4,rescr_ari4) 
c(resms_sim4,rescr_sim4) 