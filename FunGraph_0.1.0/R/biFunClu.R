#' @title generate inital parameters for biFunClu
#' @importFrom stats kmeans
#' @param data dataframe contain numeric values for further analysis
#' @param k scalar indicate the number of cluster
#' @param legendre_order scalar of legendre polynomials
#' @param times vector of time point
#' @param init_SAD_par vector of length 4 for biSAD1 covariance matrix
#' @return list contain initial parameters for biFunClu
#' @export
get_init_par <- function(data, k, legendre_order, times, init_SAD_par = c(1.06,0.25,1.15,0.18)){
  #get inital pars based on k-means
  n = ncol(data)/2
  init_cluster <- kmeans(data,centers = k,iter.max = 1e3)
  cuM <- init_cluster$centers

  init_curve_par <- cbind(t(sapply(1:k,function(c)get_legendre_par(cuM[c,1:n],legendre_order,times))),
                          t(sapply(1:k,function(c)get_legendre_par(cuM[c,(n+1):(2*n)],legendre_order,times))))
  init_pro <- table(init_cluster$cluster)/nrow(data)
  return_object <- list(data = data,
                        init_SAD_par = init_SAD_par,
                        init_curve_par = init_curve_par,
                        init_pro = init_pro,
                        times = times)
  return(return_object)
}

#' @title biFunClu function
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats optim
#' @importFrom utils write.csv
#' @param input list object from get_init_par function contain initial parameters and data
#' @param Delta scalar prevent Log-likelihood drastically change
#' @param itermax scalar constrain the maximum number of iteracion
#' @param write_data write result table for further analysis
#' @return Log-likelihood value
#' @export
get_cluster <- function(input, itermax = 100, Delta = 100,  write_data = TRUE){
  data = input$data
  times = input$times
  n = length(times)
  d = ncol(data)/2
  k = length(input$init_pro)
  legendre_order = 1/2*ncol(input$init_curve_par) - 1
  iter = 1
  mle <- function(par,data,prob){
    par1 <- par[1:4]
    par2 <- matrix(par[-c(1:4)],nrow = k,ncol = (legendre_order+1)*2)
    mu <- t( sapply(1:k, function(c) c(legendre_fit(par2[c,1:(legendre_order+1)],times),
                                       legendre_fit(par2[c,(ncol(par2)/2+1):ncol(par2)],times)) ) )
    temp_S <- sapply(1:k,function(c) dmvnorm(data,
                                             mu[c,],
                                             get_biSAD1(par1,n))*prob[c] )
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  cat(paste0("Start biFunClu Calculation ","\n","Cluster_number=",k," Legendre_order=", legendre_order,'\n'))
  while ( Delta > 1 && iter <= itermax ) {
    # initiation
    if(iter == 1){
      init_SAD_par <- input$init_SAD_par
      init_curve_par <- input$init_curve_par
      pro <- input$init_pro
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(old_par,data,pro)
    mu <- t( sapply(1:k, function(c)c(legendre_fit(init_curve_par[c,1:(legendre_order+1)],times),
                                      legendre_fit(init_curve_par[c,(ncol(init_curve_par)/2+1):ncol(init_curve_par)],times))) )
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             mu[c,],
                                             get_biSAD1(init_SAD_par,n))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_SAD_par <- new_par$par[1:4]
    init_curve_par <- matrix(new_par$par[-c(1:4)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    if (Delta > 20000)
      break
    cat("iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  }

  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)

  cluster <- apply(omega,1,which.max)
  clustered_data1 <- cbind(data[,1:d],cluster)
  clustered_data2 <- cbind(data[,(d+1):(2*d)],cluster)
  return_object <- list(SAD_par = init_SAD_par,
                        curve_par = init_curve_par,
                        pro = pro,
                        Loglikelihood = LL_mem,
                        BIC = BIC,
                        clustered_data = list(clustered_data1, clustered_data2))
  cat("Finish biFunClu Calculation")
  if (write_data==TRUE) {
    write.csv(clustered_data1,file = 'cluster_data1.csv')
    write.csv(clustered_data2,file = 'cluster_data2.csv')
    write.csv(init_curve_par,file = 'cluster_par.csv')
    cat('cluster result data written',sep="\n")
  }
  return(return_object)
}

#' @title plot cluster result base plot
#' @import ggplot2
#' @importFrom reshape2 melt
#' @param clustered_data clustered data with rownames in first column
#' @return a plot
#' @export
get_cluster_base_plot <- function(clustered_data){
  clustered_data <- cbind(rownames(clustered_data),clustered_data)
  colnames(clustered_data) <- c("marker", 1:(ncol(clustered_data)-2), "cluster")
  long_df <- melt(clustered_data,c("marker","cluster"))
  colnames(long_df) <- c("marker","cluster","time","effect")
  p <-  ggplot()+
    geom_line(long_df,mapping=
                aes(as.numeric(as.character(long_df$time)),long_df$effect,group=long_df$marker,
                                               colour = as.character(long_df$cluster)),alpha=1)+
    facet_wrap(long_df$cluster,scales = "fixed")+
    theme(legend.position="none") + xlab("Time")+ylab("genetic_effect")
  return(p)
}
