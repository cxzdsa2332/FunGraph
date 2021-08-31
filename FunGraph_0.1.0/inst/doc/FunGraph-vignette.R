## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, fig.height = 3.5
)

## ----setup--------------------------------------------------------------------
library(FunGraph)

## ---- eval=FALSE--------------------------------------------------------------
#  View(geno)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(geno[,1:10]))

## ---- eval=FALSE--------------------------------------------------------------
#  View(pheno)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(pheno[,1:10]))

## -----------------------------------------------------------------------------
get_mean_curve_plot(pheno_df = pheno, times = 1:14)

#Notice wrong given inital parameters results in incorrect mean curve(blue line)
get_mean_curve_plot(pheno_df = pheno, times = 1:14, init_sd_par = c(2,0.1,2,0.1))

## -----------------------------------------------------------------------------
#geno[,-1:-2] remove redundant information for calculation
get_LR_effect(marker = 1, 
              geno_df = geno[,-1:-2], 
              pheno_df = pheno, 
              times = 1:14)

## ----eval = FALSE-------------------------------------------------------------
#  get_biFunMap_result(geno_df = geno[,-1:-2],
#                      pheno_df = pheno,
#                      times = 1:14)

## -----------------------------------------------------------------------------
get_manh_plot(geno_df = geno, LR = LR)

## -----------------------------------------------------------------------------
get_generic_effect_plot(generic_effect = generic_effect, number = 10)

## ----eval = FALSE-------------------------------------------------------------
#  #example use permutation 5 times, each permutation use first 3 SNPs
#  get_permutation(n = 5,
#                  geno_df = geno[1:3,-1:-2],
#                  pheno_df = pheno,
#                  times = 1:14)

## -----------------------------------------------------------------------------
set.seed(2021) #use same initial value

input <- get_init_par(data = generic_effect[1:100,],
                         k = 5,
                         legendre_order = 3,
                         times = 1:14)

c1 <- get_cluster(input = input, itermax = 10)

#plot the result
get_cluster_base_plot(c1$clustered_data[[1]])

## -----------------------------------------------------------------------------
module_data <- get_module_data(data_par = c1$curve_par, times = 1:14)

#select ck data for variable selection
get_interaction(data = module_data[[1]], col = 1)

#alternatively, use a cluster data(ck) for variable selection
get_interaction(data = c1$clustered_data[[1]][,-ncol(c1$clustered_data[[1]])], 
                col = 1, 
                reduction = TRUE)

## -----------------------------------------------------------------------------
module_data <- get_module_data(data_par = c1$curve_par, times = 1:14)
#for module, just input module data
module_ode1 <- get_ode_par(data = module_data[[1]],
                          times = 1:14,
                          order = 3,
                          reduction = FALSE,
                          parallel = FALSE)
#the result should be further processed
module_net1 <- get_all_net(module_ode1)

## -----------------------------------------------------------------------------
get_decomposition_plot(module_ode1,module_net1,1)
get_decomposition_plot(module_ode1,module_net1,3)

## -----------------------------------------------------------------------------
get_net_output(module_ode1,module_net1,write_data = FALSE)

## -----------------------------------------------------------------------------
#acquire max_effect to control color
max_effect <- cbind(get_max_effect(module_net1),get_max_effect(module_net1))
network_plot(k = module_net1,
             title = 'CK',
             max_effect = max_effect,
             save_plot = FALSE)

