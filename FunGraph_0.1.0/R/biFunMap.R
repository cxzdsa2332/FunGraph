#' @title calculate Loglikelihood under h0
#' @importFrom mvtnorm dmvnorm
#' @param y_all dataframe contain all phenotypic data except missing marker(drop marker=9)
#' @param par vector for modified logistic growth curve
#' @param times vector for the time point
#' @return the Loglikelihood value
L0 <- function(y_all, par, times){
  n = length(times)
  mu <- c(get_mu(par[1:5],times),get_mu(par[6:10],times))
  SAD1<- get_biSAD1(par[11:14],n)
  L0 <- -sum(dmvnorm(y_all,mu,SAD1,log = T))
  return(L0)
}

#' @title calculate Loglikelihood under h1
#' @importFrom mvtnorm dmvnorm
#' @param pheno_0 dataframe contain all phenotypic data with marker=0
#' @param pheno_1 dataframe contain all phenotypic data with marker=1
#' @param par vector for modified logistic growth curve
#' @param times vector for the time point
#' @return the Loglikelihood value
L1 <- function(pheno_0, pheno_1, par, times){
  n = length(times)
  SAD1 <- get_biSAD1(par[21:24],n)
  L_0 <- -sum(dmvnorm(pheno_0,c(get_mu(par[1:5],times),get_mu(par[6:10],times)),SAD1,log = T))
  L_1 <- -sum(dmvnorm(pheno_1,c(get_mu(par[11:15],times),get_mu(par[16:20],times)),SAD1,log = T))
  LL <- L_0 + L_1
  return(LL)
}

#' @title calculate Loglikelihood under h1
#' @importFrom mvtnorm dmvnorm
#' @param pheno_0 dataframe contain all phenotypic data with marker=0
#' @param pheno_1 dataframe contain all phenotypic data with marker=1
#' @param pheno_2 dataframe contain all phenotypic data with marker=2
#' @param par vector for modified logistic growth curve
#' @param times vector for the time point
#' @return the Loglikelihood value
L2 <- function(pheno_0, pheno_1, pheno_2, par, times){
  n = length(times)
  SAD1 <- get_biSAD1(par[31:34],n)
  L_0 <- -sum(dmvnorm(pheno_0,c(get_mu(par[1:5],times),get_mu(par[6:10],times)),SAD1,log = T))
  L_1 <- -sum(dmvnorm(pheno_1,c(get_mu(par[11:15],times),get_mu(par[16:20],times)),SAD1,log = T))
  L_2 <- -sum(dmvnorm(pheno_2,c(get_mu(par[21:25],times),get_mu(par[26:30],times)),SAD1,log = T))
  LL <- L_0 + L_1 + L_2
  return(LL)
}

#' @title perform biFunMap for a SNP
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats optim
#' @param marker scalar indicate which SNP to use(row number)
#' @param geno_df dataframe of Genotypic data
#' @param pheno_df dataframe of Phenotypic data
#' @param times vector for the time point
#' @param init_sd_par vector of length four for biSAD1 covariance matrix
#' @param simplify default is FALSE, only used in permutation to return less result
#' @return the LR value, generic effect and calculated parameters
#' @export
get_LR_effect <- function(marker, geno_df, pheno_df, times,
                          init_sd_par = c(0.95,12,1.02,8),simplify = FALSE){
  n = length(times)

  pheno_0 <- pheno_df[which(geno_df[marker,]==0),]
  pheno_1 <- pheno_df[which(geno_df[marker,]==1),]
  pheno_2 <- pheno_df[which(geno_df[marker,]==2),]
  pheno_9 <- pheno_df[which(geno_df[marker,]==9),]
  y_all <- rbind(pheno_0,pheno_1,pheno_2)

  get_init_pars <- function(par){
    y <- as.numeric(colMeans(y_all))
    y1 <- c(get_mu(par[1:5],times),get_mu(par[6:10],times))
    ssr <- sum((y1-y)^2)
  }
  init_curve_par <- optim(c(mean(y_all[,14]),10,0.5,-2,2,mean(y_all[,28]),10,2,-2,1),get_init_pars)$par

  init_par <- c(init_curve_par,init_sd_par)

  NH_0 <- optim(init_par, L0, y_all = y_all, times = times)

  if (nrow(pheno_2)==0) {
    h1_pars <- c(NH_0$par[1:10],NH_0$par[1:10],NH_0$par[11:14])
    NH_1 <- optim(h1_pars, L1, pheno_0 = pheno_0, pheno_1 = pheno_1, times = times)
  }else{
    h2_pars <- c(NH_0$par[1:10],NH_0$par[1:10],NH_0$par[1:10],NH_0$par[11:14])
    NH_1 <- optim(h2_pars, L2, pheno_0 = pheno_0, pheno_1 = pheno_1,  pheno_2 = pheno_2,
                  times = times)
  }
  LR <- 2*(NH_0$value - NH_1$value)

  get_genetic_effect <- function(times){
    effect1 <- nrow(pheno_0)*(get_mu2(NH_1$par[1:10],times)^2) +
      nrow(pheno_1)*(get_mu2(NH_1$par[11:20],times)^2)
    effect2 <- nrow(pheno_0)*get_mu2(NH_1$par[1:10],times) +
      nrow(pheno_1)*get_mu2(NH_1$par[11:20],times)
    if (nrow(pheno_2)!=0) {
      effect1 <- effect1+nrow(pheno_2)*(get_mu2(NH_1$par[21:30],times)^2)
      effect2 <- effect2+nrow(pheno_2)*get_mu2(NH_1$par[21:30],times)
    }
    effect <- sqrt( effect1/nrow(y_all) - (effect2/nrow(y_all))^2 )
    return(effect)
  }
  generic_effect <- get_genetic_effect(times)

  if (simplify != FALSE) {
    return_object <- LR
  }
  else{
    return_object <- list(LR = LR,
                          generic_effect = generic_effect,
                          NH0_pars = NH_0$par,
                          NH1_pars = NH_1$par)
  }
  return(return_object)
}

#' @title permutation test in parallel mode
#' @import parallel
#' @import pbapply
#' @param n scalar the number of times permutation performed, default = 1000
#' @param geno_df dataframe of Genotypic data
#' @param pheno_df dataframe of Phenotypic data
#' @param times vector for the time point
#' @return a list of LR values, each list contain LR results equal to SNP numbers
#' @export
get_permutation <- function(n = 1000, geno_df, pheno_df, times){
  cat('Start permutation test',sep="\n")
  permutation_results <- list()
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(mvtnorm)})
  clusterExport(cl, c("geno_df","pheno_df","permutation_results","n"),envir=environment())
  for (i in 1:n) {
    new_pheno <- data.frame(pheno_df[sample(nrow(pheno_df[,])),],row.names = NULL)
    #use scrambled phenotypic data to run biFunMap
    permutation_results[[i]] <- pbsapply(cl=cl,1:nrow(geno_df),function(c)
      get_LR_effect(marker=c,geno_df,new_pheno,times=times,simplify = TRUE))
    cat(paste0('Finish permutation ',i,' times'),sep="\n")
  }
  stopCluster(cl)
  return(permutation_results)
}

#' @title perform biFunMap for all SNP in parallel mode
#' @import parallel
#' @import pbapply
#' @importFrom utils write.csv
#' @param geno_df dataframe of Genotypic data
#' @param pheno_df dataframe of Phenotypic data
#' @param times vector for the time point
#' @param write_data write result table for further analysis
#' @return a list of LR values
#' @export
get_biFunMap_result <- function(geno_df, pheno_df, times, write_data = TRUE){
  cat('Start biFunMap Calculation')
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {require(mvtnorm)})
  clusterExport(cl, c("geno_df","pheno_df","times"),envir = environment())
  result1 <- pblapply(cl=cl,1:nrow(geno_df),function(c)
    get_LR_effect(marker=c,geno_df,pheno_df,times=times))
  stopCluster(cl)

  LR_result <- do.call(c,lapply(1:length(result1),function(c) result1[[c]][[1]]))
  generic_effect <- do.call(rbind,lapply(1:length(result1),function(c) result1[[c]][[2]]))
  rownames(generic_effect) <- rownames(geno_df)
  colnames(generic_effect) <- colnames(pheno_df)

  return_object <- list(LR_result = LR_result,
                        generic_effect = generic_effect,
                        FunMap_pars = result1)

  if (write_data==TRUE) {
    write.csv(generic_effect,file = 'generic_effect.csv')
    cat('generic_effect data written',sep="\n")
  }
  return(return_object)
}

#' @title helper function to convert plot data
#' @param data dataframe of Phenotypic data
#' @return a long dataframe for plot
get_plot_data <- function(data){
  d = ncol(data)/2
  tmp1 <- cbind(data[,1:d],'ck')
  colnames(tmp1) <- c(1:d,'condition')
  tmp2 <- cbind(data[,(d+1):(2*d)],'stress')
  colnames(tmp2) <- c(1:d,'condition')
  tmp <- rbind(tmp1,tmp2)
  return(tmp)
}

#' @title plot mean curve for biFunMap results
#' @import ggplot2
#' @import cowplot
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats optim
#' @importFrom reshape2 melt
#' @param pheno_df dataframe of Phenotypic data
#' @param times vector for the time point
#' @param init_sd_par vector of length four for biSAD1 covariance matrix
#' @return a plot with mean curve for biFunMap
#' @export
get_mean_curve_plot <- function(pheno_df, times, init_sd_par = c(0.95,12,1.02,8)){
  d = ncol(pheno_df)/2
  y_all = pheno_df
  get_init_pars <- function(par){
    y <- as.numeric(colMeans(y_all))
    y1 <- c(get_mu(par[1:5],times),get_mu(par[6:10],times))
    ssr <- sum((y1-y)^2)
  }
  init_curve_par <- optim(c(mean(y_all[,14]),10,0.5,-2,2,mean(y_all[,28]),10,2,-2,1),get_init_pars)$par

  init_par <- c(init_curve_par,init_sd_par)

  NH_0 <- optim(init_par, L0, y_all = y_all, times = times)

  col=c("#F8766D","#619CFF")

  CK <- cbind(rownames(y_all),get_plot_data(y_all)[1:nrow(y_all),])
  colnames(CK) <- c('id',1:d,'condition')
  CK <- melt(CK,id.vars = c('id','condition'))
  CK_mean <- data.frame(cbind(1:d,get_mu(NH_0$par[1:5],times)))

  p1 <- ggplot()+geom_line(CK,mapping = aes(variable,value,group=id),color=col[2],alpha=0.15)+
    geom_line(CK_mean,mapping = aes(X1,X2),color=col[2],size=2)+
    theme_bw()+theme(panel.grid =element_blank())+
    xlab(NULL)+ylab('Phenotype')+
    annotate('text',x=7,y=190,label='Control',size=4,col=col[2])+
    theme(plot.margin = unit(c(0.5,0,0,0.5),"lines"))+
    theme(
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 10))

  salt <- cbind(rownames(y_all),get_plot_data(y_all)[(nrow(y_all)+1):(2*nrow(y_all)),])
  colnames(salt) <- c('id',1:d,'condition')
  salt <- melt(salt,id.vars = c('id','condition'))
  salt_mean <- data.frame(cbind(1:d,get_mu(NH_0$par[6:10],times)))

  p2 <- ggplot()+geom_line(salt,mapping = aes(variable,value,group=id),color=col[1],alpha=0.15)+
    geom_line(salt_mean,mapping = aes(X1,X2),color=col[1],size=2)+
    theme_bw()+theme(panel.grid =element_blank())+
    xlab(NULL)+ylab(NULL)+
    annotate('text',x=7,y=195,label='Stress',size=4,col=col[1])+
    scale_y_continuous(position = "right")+
    theme(plot.margin = unit(c(0.5,0.5,0,0),"lines"))+
    theme(
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 10))

  pp <- plot_grid(p1,p2,ncol=2)
  pp <- add_sub(pp, "Time (day)", hjust = 0.3,vjust=0)
  ggdraw(pp)
}


#' @title manhattan plot
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @param geno_df dataframe of Genotypic data(contain addioitonal information)
#' @param LR_result vector of length = nrow(geno_df), contain all LR value
#' @param threshold scalar for LR threthold
#' @return manhattan plot
#' @export
get_manh_plot <- function(geno_df, LR_result, threshold = 20) {
  geno_df$LR_result = LR_result
  geno_df$Marker_ID = rownames(geno_df)
  get_max <- function(i){
    chr <- paste0('lg',i)
    return(max(geno_df[geno_df$Linkage==chr,]$`Genetic_Distances(cM)`))
  }
  k = length(table(geno_df$Linkage))

  chr_max <- sapply(1:k, function(c)get_max(c))
  chr_max2 <- c(0,cumsum(chr_max))
  geno_df$x <- NA
  for (i in 1:k) {
    chr <- paste0('lg',i)
    geno_df[geno_df$Linkage==chr,]$x <- geno_df[geno_df$Linkage==chr,]$`Genetic_Distances(cM)`+chr_max2[i]
  }

  get_axisdf <- function(i){
    chr <- paste0('lg',i)
    return(1/2*(min(geno_df[geno_df$Linkage==chr,]$x)+max(geno_df[geno_df$Linkage==chr,]$x)))
  }
  axisdf <- sapply(1:k, function(c)get_axisdf(c))

  filtered_df <- geno_df[geno_df$LR_result>threshold,]

  p <- ggplot(geno_df, aes(x=x, y=LR_result)) +
    geom_point( aes(color=Linkage), alpha=0.8, size=2) +
    scale_color_manual(values = rep(c("#E2709A", "#CB4577",
                                      "#BD215B", "#970F42",
                                      "#75002B"), 22)) +
    scale_x_continuous(labels = 1:k, breaks= axisdf) +
    scale_y_continuous(expand = c(0, 0),limits = c(0,threshold+5)) +
    theme_bw() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),panel.grid =element_blank(),
          axis.line = element_line(color = "black")) +
    xlab("Chromosome") + ylab('LR')+
    geom_hline(yintercept = threshold, linetype="dashed")+
    geom_point(data=filtered_df, color="orange", size=2)+
    geom_label_repel( data=filtered_df,
                      aes(label=Marker_ID), size=3)
  return(p)
}

#' @title plot generic effect of random choose number = number
#' @import ggplot2
#' @importFrom reshape2 melt
#' @param generic_effect dataframe of calculated generic effect data
#' @param number scalar of number of SNPs' effect curve want to show(must <= nrow(generic_effect))
#' @return line plot of generic effect
#' @export
get_generic_effect_plot <- function(generic_effect, number){
  d = ncol(generic_effect)/2
  df <- generic_effect[sample(nrow(generic_effect), number),]
  df_ck <- df[,1:d]
  colnames(df_ck) <- 1:d
  df_ck <- cbind(rownames(df_ck),df_ck)
  df_ck <- melt(df_ck)
  df_ck <- data.frame(cbind(df_ck,'ck'))
  colnames(df_ck) <- c('id','time','effect','condition')

  df_salt <- df[,(d+1):(2*d)]
  colnames(df_salt) <- 1:d
  df_salt <- cbind(rownames(df_salt),df_salt)
  df_salt <- melt(df_salt)
  df_salt <- data.frame(cbind(df_salt,'salt'))
  colnames(df_salt) <- c('id','time','effect','condition')

  df2 <- rbind(df_ck,df_salt)

  mytheme <- theme_minimal() + theme(
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(colour = '#8E9775', fill=NA, size=1)
  )

  p <- ggplot(df2,mapping = aes(time,effect,group=condition,color=condition))+
    geom_line()+facet_wrap(~id,nrow=2)+scale_color_manual(values=c("#619CFF","#F8766D"))+
    mytheme + theme(legend.title = element_blank())+xlab('Time (day)')+ylab('Generic Effect')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p
}
