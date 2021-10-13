## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  tidy = TRUE,
  cache = TRUE,
  collapse = TRUE,
  dev = "png",
  fig.width = 7, 
  fig.height = 3.5
)

## ----load-packages, include=FALSE---------------------------------------------
library(FunGraph)

## ---- eval=FALSE--------------------------------------------------------------
#  View(pheno)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(pheno[,1:10]))

## ---- eval=FALSE--------------------------------------------------------------
#  View(geno)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(head(geno[,1:10]))

## ---- eval=TRUE---------------------------------------------------------------
plot( x = seq(1,5,length=100), 
      y = get_mu(c(10,5,3,-2,1),seq(1,5,length=100)), 
      type = 'l',
      main="Modified Logistic Growth Curve", 
      xlab="Times", 
      ylab="Phenotypic Values")

## ---- eval=TRUE---------------------------------------------------------------
get_SAD1_covmatrix(c(2,0.5), n = 5)


## ---- eval=TRUE---------------------------------------------------------------
get_biSAD1(c(2,0.5,0.5,2), n = 5)


## ---- fig.height=3, fig.width=7-----------------------------------------------
get_mean_curve_plot(pheno_df = pheno, times = 1:14)

#Notice wrong given inital parameters results in incorrect mean curve(blue line)
get_mean_curve_plot(pheno_df = pheno, times = 1:14, init_sd_par = c(2,0.1,2,0.1))

## ---- eval=TRUE---------------------------------------------------------------
#geno[,-1:-2] remove redundant information for calculation
get_LR_effect(marker = 1, 
              geno_df = geno[,-1:-2], 
              pheno_df = pheno, 
              times = 1:14)

## ---- eval=FALSE--------------------------------------------------------------
#  result1 <- get_biFunMap_result(geno_df = geno[,-1:-2],
#                      pheno_df = pheno,
#                      times = 1:14)

## ---- eval=TRUE---------------------------------------------------------------
get_manh_plot(geno_df = geno, LR = LR, threshold = 20)

## ----eval = TRUE--------------------------------------------------------------
#example use permutation 5 times, each permutation use first 3 SNPs
perm <- get_permutation(n = 5, 
                        geno_df = geno[1:3,-1:-2],
                        pheno_df = pheno,
                        times = 1:14)

threthold <- sort(do.call(c, lapply(perm, max)))

## ---- eval=TRUE, fig.width = 10-----------------------------------------------
get_genetic_effect_plot(genetic_effect = genetic_effect, number = 10)

## -----------------------------------------------------------------------------
set.seed(2021) #use same initial value

input <- get_init_par(data = genetic_effect[1:1000,],
                         k = 5,
                         legendre_order = 4,
                         times = 1:14)

## -----------------------------------------------------------------------------
input$init_SAD_par

## -----------------------------------------------------------------------------
plot( x = seq(1,14), 
      y = legendre_fit(input$init_curve_par[1,],seq(1,14,)),
      type = 'l',
      main="Initial parameters fit for cluster 1", 
      xlab="Times", 
      ylab="Genetic Effect")


## -----------------------------------------------------------------------------
input$init_pro

## -----------------------------------------------------------------------------
c1 <- get_cluster(input = input, itermax = 100)

## ---- warning=FALSE-----------------------------------------------------------
get_cluster_base_plot(c1$clustered_data[[1]])

## ---- warning=FALSE-----------------------------------------------------------
set.seed(2021) 

BIC <- get_BIC(data = genetic_effect[1:1000,], 
             order = 4,
             times = 1:14,
             rep = 3,
             min_cluster = 2,
             max_cluster = 6,
             itermax = 10)

plot(x = BIC$k, y = BIC$BIC,
     type = 'l',
     main="BIC", 
     xlab="k", 
     ylab="BIC")
abline( v = BIC$k[which.min(BIC$BIC)],col="red", lwd=3, lty=2)

## ---- warning=FALSE-----------------------------------------------------------
module_data <- get_module_data(data_par = c1$curve_par, times = 1:14)

#select ck data for variable selection
get_interaction(data = module_data[[1]], col = 1)


## ---- warning=FALSE-----------------------------------------------------------
#alternatively, use a cluster data(ck) for variable selection
get_interaction(data = c1$clustered_data[[1]][,-ncol(c1$clustered_data[[1]])], 
                col = 1, 
                reduction = TRUE)

## ---- warning=FALSE-----------------------------------------------------------
#Prepare genetic effect dataset for modules
module_data <- get_module_data(data_par = c1$curve_par, times = 1:14)
#Solve ODE between modules
module_ode1 <- get_ode_par(data = module_data[[1]],
                           times = 1:14,
                           order = 3,
                           reduction = FALSE,
                           parallel = FALSE)
#the result should be further processed
module_net1 <- get_all_net(module_ode1)

## ---- warning=FALSE-----------------------------------------------------------
get_decomposition_plot(module_ode1,module_net1,1)

## ---- fig.height=6, fig.width=7, warning=FALSE--------------------------------
#max_effect control the overall plot colour (eg. ck and stress network)
max_effect <- cbind(get_max_effect(module_net1),get_max_effect(module_net1))

network_plot(k = module_net1, title = 'CK', max_effect = max_effect, save_plot = FALSE)

## ---- message=FALSE, warning=FALSE--------------------------------------------
m1_ck <- get_subset_data(data = c1$clustered_data[[1]], cluster = 1 )
m1_ck_ode <- get_ode_par(data = m1_ck, times = 1:14, order = 3, reduction = TRUE, parallel = TRUE)
m1_ck_net <- get_all_net(m1_ck_ode)
get_decomposition_plot(m1_ck_ode,m1_ck_net,2)

## ---- fig.height=6, fig.width=7, warning=FALSE--------------------------------
max_effect <- cbind(get_max_effect(m1_ck_net),get_max_effect(m1_ck_net))

network_plot(k = m1_ck_net, title = 'Module1_CK', max_effect = max_effect, save_plot = FALSE)

## -----------------------------------------------------------------------------
table(c1$clustered_data[[1]]$cluster)

## ---- warning=FALSE-----------------------------------------------------------
m1 <- cbind(get_subset_data(data = c1$clustered_data[[1]], cluster = 1 ),
             get_subset_data(data = c1$clustered_data[[2]], cluster = 1 ))

set.seed(2021) #use same initial value
input2 <- get_init_par(data = m1, k = 4, legendre_order = 4, times = 1:14)
c2 <- get_cluster(input = input2)
#view stress data
get_cluster_base_plot(c2$clustered_data[[2]])

table(c2$clustered_data[[2]]$cluster)

## ---- message=FALSE, warning=FALSE--------------------------------------------
submodule1_data <- get_module_data(data_par = c2$curve_par, times = 1:14)
submodule1_ode1 <- get_ode_par(data = submodule1_data[[1]],
                               times = 1:14,
                               order = 3,
                               reduction = FALSE,
                               parallel = TRUE)
submodule1_net1 <- get_all_net(submodule1_ode1)
max_effect1 <- cbind(get_max_effect(submodule1_net1),get_max_effect(submodule1_net1))

## ---- fig.height=6, fig.width=7-----------------------------------------------
network_plot(k = submodule1_net1, title = "Submodule1_CK", max_effect = max_effect1, save_plot = FALSE)

## ---- message=FALSE, warning=FALSE--------------------------------------------

m1_1_ck <- get_subset_data(data = c2$clustered_data[[1]], cluster = 2)
m1_1_ck_ode <- get_ode_par(data = m1_1_ck, times = 1:14, order = 3, reduction = TRUE, parallel = TRUE)
m1_1_ck_net <- get_all_net(m1_1_ck_ode)
max_effect2 <- cbind(get_max_effect(m1_1_ck_net),get_max_effect(m1_1_ck_net))

## ---- fig.height=6, fig.width=7-----------------------------------------------
network_plot(k = m1_1_ck_net, title = "Submodule1_1_CK", 
             max_effect = max_effect2, save_plot = FALSE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

