rm(list=ls())
#load('FunGraph.Rdata')
#1.Load and clean data------------------------------------------------------------------------------------------
source('FunGraph.R')
source('FunGraph_plot.R')
library(mvtnorm)
library(pbapply)
library(parallel)
library(orthopolynom)
library(glmnet)
library(ggplot2)
library(reshape2)

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
biFunMap_results <- get_biFunMap_result(geno_df[,-3:-1])
#Out put Figure1
Figure1 <- get_figure1(get_p1_v2(times=1:14,data=pheno_df,geno_df2=geno_df[,-3:-1]),
                       get_plot1b(cbind(geno_df[,1:3],biFunMap_results[[1]])),
                       get_plot1c(biFunMap_results[[2]]))
#3.run biFunCluster---------------------------------------------------------------------------------------------
#cluster(k)=15,LOP=4
set.seed(2021)
biFunClu_intial_pars <- get_init_par(data=biFunMap_results$generic_effect,k=15,legendre_order=4)
biFunClu_results <- get_cluster(data=biFunMap_results$generic_effect,k=15,
                                input=biFunClu_intial_pars,legendre_order=4)
#Out put Figure2_B,use example_bic files
Figure2 <- get_figure2(BIC_df=read.csv(file='example_bic.csv'),cluster_result=biFunClu_results)
table(biFunClu_results$clustered_ck$cluster)
#4.LASSO-based variable selection-------------------------------------------------------------------------------
#demonstration for modules,k=15
module_effect <- t(sapply(1:15, function(c)legendre_fit(biFunClu_results$curve_par[c,1:5])))
rownames(module_effect) <- paste0('Module',1:15)
module_relationship <- lapply(1:15,function(c) get_interaction(t(module_effect),c))

#5.run ODE solving----------------------------------------------------------------------------------------------
legendre_order=4
#ODE of modules
module_result <- get_module_result(k=15,times=seq(1,14,length=30),order=3,cluster_result = biFunClu_results)
#plot Figure3_A
Figure3A <- get_fig3a_v2(k=15,times=seq(1,14,length=30),all_net=module_result[[2]],
                         cluster_result = biFunClu_results)
#output xlsx files for Cytoscape
library(writexl)
write_xlsx(list(module_result[[1]][[1]][[1]]),'ODE_ck1.xlsx')
write_xlsx(list(module_result[[1]][[1]][[2]]),'ODE_ck2.xlsx')
write_xlsx(list(module_result[[1]][[2]][[1]]),'ODE_salt1.xlsx')
write_xlsx(list(module_result[[1]][[2]][[2]]),'ODE_salt2.xlsx')
#ODE of module 7
SNP_result <- get_SNP_result(cluster=7,time=seq(1,14),order=3,cluster_result = biFunClu_results)
#plot Figure3_C
Figure3C <- get_figure3c(k=15,times=seq(1,14),all_net=SNP_result[[2]])
#5.if necessary,select a module to run sub_biFunclusr and ODE solving-------------------------------------------
#cluster(k)=8,LOP=4,further cluster module 13,key SNP is 'nn_np_2890'
sub_cluster = 13

sub_cluster_data = cbind(biFunClu_results$clustered_ck[which(biFunClu_results$clustered_ck$cluster==sub_cluster),-15],
                         biFunClu_results$clustered_salt[which(biFunClu_results$clustered_salt$cluster==sub_cluster),-15])
biFunClu_sub_intial_pars <- get_init_par(data=sub_cluster_data,k=8,legendre_order=4)
biFunClu_sub_results <- get_cluster(data=sub_cluster_data,k=8,
                                    input=biFunClu_sub_intial_pars,legendre_order=4)
sub_module_result <- get_module_result(k=8,times=seq(1,14,length=30),order=3,
                                       cluster_result = biFunClu_sub_results)

sub_SNP_result <- get_SNP_result(cluster=7,time=seq(1,14),order=3,
                                 cluster_result = biFunClu_sub_results)
#save.image('FunGraph.Rdata')


