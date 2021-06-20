rm(list=ls())
#1.Load and clean data------------------------------------------------------------------------------------------
source('FunGraph.R')
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)

geno_df <- read.csv('geno_df.csv',check.names = F)
ck_df <- read.csv('CK_df.csv')
salt_df <- read.csv('salt_df.csv')

#using selected pheno(N_XX)
pheno_id <- c(1,c(3:16))

#select sum_root and pheno_id
ck_df <- ck_df[which(ck_df$roottype=='sum_root'),pheno_id]
salt_df <- salt_df[which(salt_df$roottype=='sum_root'),pheno_id]

#get_mean value for duplicated id
get_data <- function(df){
  tmp <- do.call(cbind,sapply(1:(ncol(df)-1),function(c)aggregate( df[,c+1] ~ id, df, mean)))
  tmp <- tmp[,c(1,seq(2,ncol(tmp),2))]
  colnames(tmp) <- colnames(df)
  return(tmp)
}
pheno_df <- merge(get_data(ck_df),get_data(salt_df),by='id')

#geno_with_pheno
geno_df <- geno_df[,c(1:3,sapply(1:nrow(pheno_df),
                                 function(c)which(as.character(pheno_df$id)[c]==colnames(geno_df))))]

#2.run biFunMap-------------------------------------------------------------------------------------------------
FunMap_results <- get_FunClu_rewsult(geno_df[,-3:-1])
#save.image(file = 'FunMap_results.Rdata')

#3.run biFunCluster---------------------------------------------------------------------------------------------
#load('FunMap_results.Rdata')
#cluster=15
set.seed(2021)
input_pars <- get_init_par(data=FunMap_results$generic_effect,k=15,legendre_order=4)
cluster_result <- get_cluster(data=FunMap_results$generic_effect,k=15,input=input_pars,legendre_order=4)
#save.image('FunClu_results.Rdata')

#4.run ODE solving----------------------------------------------------------------------------------------------
#load('FunClu_results.Rdata')
legendre_order=4
#ODE of modules
module_result <- get_module_result(k=15,times=seq(1,14,length=30),order=3)
#output csv file for Cytoscape
write.csv(module_result[[1]][[1]][[1]],file = 'ck1.csv',row.names = FALSE)
write.csv(module_result[[1]][[1]][[2]],file = 'ck2.csv',row.names = FALSE)
write.csv(module_result[[1]][[2]][[1]],file = 'salt1.csv',row.names = FALSE)
write.csv(module_result[[1]][[2]][[2]],file = 'salt2.csv',row.names = FALSE)
#ODE of module 7
submodule_result <- get_submodule_result(cluster=7,time=seq(1,14),order=3)
#save.image("ODE_results.Rdata")
