get_plot1a <- function(marker,times,data,geno_df2){
  requiredPackages = c("ggplot2","reshape2","cowplot",'mvtnorm')
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  pheno_df = data; t = times
  pheno_0 <- pheno_df[which(geno_df2[marker,]==0),-1]
  pheno_1 <- pheno_df[which(geno_df2[marker,]==1),-1]
  pheno_2 <- pheno_df[which(geno_df2[marker,]==2),-1]
  pheno_9 <- pheno_df[which(geno_df2[marker,]==9),-1]
  y_all <- rbind(pheno_0,pheno_1,pheno_2)
  #miu
  get_miu <- function(miu_par){
    miu <- miu_par[1]/(1 + miu_par[2] * exp(-miu_par[3] * t)) - (miu_par[4] * exp(-miu_par[5] * t))
    miu
  }
  #bi(SAD1)
  
  get_biSAD1 <- function(par){
    n=length(t)
    get_SAD1_covmatrix <- function(par){
      phi <- par[1]; gamma <- par[2]; 
      sigma <- array(dim=c(n,n))
      #formula 1, diag element
      diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
      #formula 2, non-diag element
      sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
      sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
      return(gamma^2*sigma)
    }
    sig1 <- get_SAD1_covmatrix(par[1:2])
    sig2 <- get_SAD1_covmatrix(par[3:4])
    sig12 <- array(0, dim=c(n,n))
    sigma1 <- cbind(sig1,sig12)
    sigma2 <- cbind(sig12,sig2)
    sigma <- rbind(sigma1,sigma2)
    return(sigma)
  }
  
  L0 <- function(par){
    miu = c(get_miu(par[1:5]),get_miu(par[6:10]))
    SAD1 = get_biSAD1(par[11:14])
    L0 = -sum(dmvnorm(y_all,miu,SAD1,log = T))
    L0
  }
  L1 <- function(par){
    SAD1 <- get_biSAD1(par[21:24])
    L_0 <- -sum(dmvnorm(pheno_0,c(get_miu(par[1:5]),get_miu(par[6:10])),SAD1,log = T))
    L_1 <- -sum(dmvnorm(pheno_1,c(get_miu(par[11:15]),get_miu(par[16:20])),SAD1,log = T))
    LL <- L_0+L_1
    return(LL)
  }
  L2 <- function(par){
    SAD1 <- get_biSAD1(par[31:34])
    L_0 <- -sum(dmvnorm(pheno_0,c(get_miu(par[1:5]),get_miu(par[6:10])),SAD1,log = T))
    L_1 <- -sum(dmvnorm(pheno_1,c(get_miu(par[11:15]),get_miu(par[16:20])),SAD1,log = T))
    L_2 <- -sum(dmvnorm(pheno_2,c(get_miu(par[21:25]),get_miu(par[26:30])),SAD1,log = T))
    LL <- L_0+L_1+L_2
    return(LL)
  }
  get_init_pars <- function(par){
    y <- as.numeric(colMeans(y_all))
    y1 <- c(get_miu(par[1:5]),get_miu(par[6:10]))
    ssr <- sum((y1-y)^2)
  }
  
  init_curve_par <- optim(c(mean(y_all[,14]),10,0.5,-2,2,mean(y_all[,28]),10,2,-2,1),get_init_pars)$par
  init_par <- c(init_curve_par,0.95,12,1.02,8)
  
  NH_0 <- optim(init_par,L0,method="Nelder-Mead")
  
  if (nrow(pheno_2)==0) {
    h1_pars <- c(NH_0$par[1:10],NH_0$par[1:10],NH_0$par[11:14])
    NH_1 <- optim(h1_pars,L1,method="Nelder-Mead")
  }else{
    h2_pars <- c(NH_0$par[1:10],NH_0$par[1:10],NH_0$par[1:10],NH_0$par[11:14])
    NH_1 <- optim(h2_pars,L2,method="Nelder-Mead")
  }
  LR <- 2*(NH_0$value - NH_1$value)
  
  col=c("#F8766D","#619CFF")
  darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
  }
  
  get_plot_data <- function(data){
    tmp1 <- cbind(data[,1:14],'ck')
    colnames(tmp1) <- c(1:14,'condition')
    tmp2 <- cbind(data[,15:28],'salt')
    colnames(tmp2) <- c(1:14,'condition')
    tmp <- rbind(tmp1,tmp2)
    return(tmp)
  }
  AA_ck <- cbind(rownames(pheno_0),get_plot_data(pheno_0)[1:nrow(pheno_0),])
  colnames(AA_ck) <- c('id',1:14,'condition')
  AA_ck <- melt(AA_ck,id.vars = c('id','condition'))
  AA_mean <- data.frame(cbind(1:14,get_miu(NH_1$par[1:5])))
  
  aa_ck <- cbind(rownames(pheno_1),get_plot_data(pheno_1)[1:nrow(pheno_1),])
  colnames(aa_ck) <- c('id',1:14,'condition')
  aa_ck <- melt(aa_ck,id.vars = c('id','condition'))
  aa_mean <- data.frame(cbind(1:14,get_miu(NH_1$par[11:15])))
  
  p1 <- ggplot()+geom_line(AA_ck,mapping = aes(variable,value,group=id),color=col[2],alpha=0.25)+
    geom_line(AA_mean,mapping = aes(X1,X2),color=darken(col[2]),size=2,linetype = "dashed")+
    geom_line(aa_ck,mapping = aes(variable,value,group=id),color=col[2],alpha=0.25)+
    geom_line(aa_mean,mapping = aes(X1,X2),color=darken(col[2]),size=2)+
    theme_bw()+theme(panel.grid =element_blank())+
    xlab(NULL)+ylab('Phenotype')+
    scale_y_continuous(limits=c(0,200))+
    annotate('text',x=7,y=190,label='Control',size=4,col=darken(col[2]))+
    theme(plot.margin = unit(c(0.5,0,0,0.5),"lines"))+
    theme(
      axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 10))
  
  AA_salt <- cbind(rownames(pheno_0),get_plot_data(pheno_0)[(nrow(pheno_0)+1):(2*nrow(pheno_0)),])
  colnames(AA_salt) <- c('id',1:14,'condition')
  AA_salt <- melt(AA_salt,id.vars = c('id','condition'))
  AA_mean2 <- data.frame(cbind(1:14,get_miu(NH_1$par[6:10])))
  
  aa_salt <- cbind(rownames(pheno_1),get_plot_data(pheno_1)[(nrow(pheno_1)+1):(2*nrow(pheno_1)),])
  colnames(aa_salt) <- c('id',1:14,'condition')
  aa_salt <- melt(aa_salt,id.vars = c('id','condition'))
  aa_mean2 <- data.frame(cbind(1:14,get_miu(NH_1$par[16:20])))
  
  p2 <- ggplot()+geom_line(aa_salt,mapping = aes(variable,value,group=id),color=col[1],alpha=0.25)+
    geom_line(AA_mean2,mapping = aes(X1,X2),color=darken(col[1]),size=2,linetype = "dashed")+
    geom_line(aa_salt,mapping = aes(variable,value,group=id),color=col[1],alpha=0.25)+
    geom_line(aa_mean2,mapping = aes(X1,X2),color=darken(col[1]),size=2)+
    theme_bw()+theme(panel.grid =element_blank())+
    xlab(NULL)+ylab(NULL)+
    annotate('text',x=7,y=195,label='Stress',size=4,col=darken(col[1]))+
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

get_plot1b <- function(df){
  requiredPackages = c("RColorBrewer","ggrepel","ggplot2")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  colnames(df)[4] <- 'LR_result'
  threshold=20
  manh_plot <- function(df, threshold) {
    get_max <- function(i){
      chr <- paste0('lg',i)
      return(max(df[df$Linkage==chr,]$`Genetic_Distances(cM)`))
    }
    k= length(table(df$Linkage))
    
    chr_max <- sapply(1:k, function(c)get_max(c))
    chr_max2 <- c(0,cumsum(chr_max))
    df$x <- NA
    for (i in 1:k) {
      chr <- paste0('lg',i)
      df[df$Linkage==chr,]$x <- df[df$Linkage==chr,]$`Genetic_Distances(cM)`+chr_max2[i]
    }
    
    get_axisdf <- function(i){
      chr <- paste0('lg',i)
      return(1/2*(min(df[df$Linkage==chr,]$x)+max(df[df$Linkage==chr,]$x)))
    }
    axisdf <- sapply(1:k, function(c)get_axisdf(c))
    
    filtered_df <- df[df$LR_result>threshold,]
    
    
    p <- ggplot(df, aes(x=x, y=LR_result)) + 
      geom_point( aes(color=Linkage), alpha=0.8, size=2) + 
      scale_color_manual(values = rep(c("#E2709A", "#CB4577", 
                                        "#BD215B", "#970F42", 
                                        "#75002B"), 22)) + 
      scale_x_continuous(label = 1:k, breaks= axisdf) + 
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
    return(p) # return the final plot
  }
  p <- manh_plot(df,20)
  p
}

get_plot1c <- function(generic_effect){
  requiredPackages = c("ggplot2","reshape2")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  set.seed(1)
  df <- generic_effect[sample(nrow(generic_effect),12),]
  df_ck <- df[,1:14]
  colnames(df_ck) <- 1:14
  df_ck <- melt(df_ck)
  df_ck <- data.frame(cbind(df_ck,'ck'))
  colnames(df_ck) <- c('id','time','effect','condition')
  
  df_salt <- df[,15:28]
  colnames(df_salt) <- 1:14
  df_salt <- melt(df_salt)
  df_salt <- data.frame(cbind(df_salt,'salt'))
  colnames(df_salt) <- c('id','time','effect','condition')
  
  df2 <- rbind(df_ck,df_salt)
  
  mytheme <- theme_minimal() + theme(
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks = element_blank(),
    #axis.title = element_blank(),
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(colour = '#8E9775', fill=NA, size=1)
  )
  
  p <- ggplot(df2,mapping = aes(time,effect,group=condition,color=condition))+
    geom_line()+facet_wrap(~id,nrow=2)+scale_color_manual(values=c("#619CFF","#F8766D"))+
    mytheme + theme(legend.title = element_blank())+xlab('Time (day)')+ylab('Generic Effect')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    theme(legend.position = "none")
  p
  
}

get_figure1 <- function(p1,p2,p3){
  
  requiredPackages = c("cowplot"，"reshape2")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  cat('Plotting Figure1','\n')
  p <- plot_grid(p1,p2,p3,labels = c('A', 'B', 'C'),ncol=1,nrow=3)
  ggsave("Figure1.pdf",p,width = 10, height = 12,units = "in")
  cat('Done')
  return(p)
}

get_figure2 <- function(BIC_df,cluster_result){
  requiredPackages = c("ggplot2","reshape2","RColorBrewer",'patchwork','orthopolynom','cowplot')
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  cat('Plotting Figure2','\n')
  clustered_df1 <- cbind(rownames(cluster_result$clustered_ck),cluster_result$clustered_ck)
  colnames(clustered_df1) <- c('ID',1:14,'cluster')
  clustered_df2 <- cbind(rownames(cluster_result$clustered_salt),cluster_result$clustered_salt)
  colnames(clustered_df2) <- c('ID',1:14,'cluster')
  
  
  long_df1 <- melt(clustered_df1,id.vars=c("ID","cluster"))
  colnames(long_df1) <- c("SNP","cluster","time","effect")
  long_df2 <- melt(clustered_df2,id.vars=c("ID","cluster"))
  colnames(long_df2) <- c("SNP","cluster","time","effect")
  
  par1 <- cluster_result$curve_par[,1:5]
  par2 <- cluster_result$curve_par[,6:10]
  
  get_mean_df <- function(data){
    legendre_fit <- function(par){
      x <- seq(1,14,length=30);legendre_order=4
      fit <- sapply(1:length(par),function(c)
        par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
      legendre_fit <- as.matrix(as.data.frame(polynomial.values(
        polynomials=fit,x=scaleX(x, u=-1, v=1))))
      x_interpolation <- rowSums(legendre_fit)
      return(x_interpolation)
    }
    mean_df <- sapply(1:nrow(data),function(c)legendre_fit(as.numeric(data[c,])))
    colnames(mean_df) <- c(1:nrow(data))
    mean_df <- melt(mean_df)
    colnames(mean_df) <- c("time","cluster","effect")
    mean_df$time <- seq(1,14,length=30)
    return(mean_df)
  }
  
  mean_df1 <- get_mean_df(par1)
  mean_df2 <- get_mean_df(par2)
  
  #cols <- brewer.pal(8,"Set3")
  cols1 <- colorRampPalette(c('#7dba91', '#59a590', '#40908e', '#287a8c', '#1c6488', '#254b7f' ))(15)
  cols2 <- colorRampPalette(c('#e98d6b', '#e3685c', '#d14a61', '#b13c6c', '#8f3371', '#6c2b6d' ))(15)
  cols1 <- c("#619CFF")
  cols2 <- c("#F8766D")
  
  darken <- function(color, factor=1.3){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
  }
  
  normalization <- function(x){((x-min(x))/(max(x)-min(x)) *0.075)+0.02}
  
  alpha <- 1/table(clustered_df1$cluster)
  alpha <- normalization(alpha)
  
  get_p <- function(i){
    cluster1 <- long_df1[which(long_df1$cluster==i),]
    cluster1$time <- as.numeric(as.character(cluster1$time))
    cluster_mean1 <- mean_df1[which(mean_df1$cluster==i),]
    cluster2 <- long_df2[which(long_df2$cluster==i),]
    cluster_mean2 <- mean_df2[which(mean_df2$cluster==i),]
    cluster2$time <- as.numeric(as.character(cluster2$time))
    
    p <- ggplot() + geom_line(cluster1,mapping=aes(time,effect,group=SNP),color=cols1,alpha=alpha[i])+
      geom_line(cluster_mean1,mapping=aes(time,effect),color=darken(cols1),size=1.5)+ 
      geom_line(cluster2,mapping=aes(time,effect,group=SNP),color=cols2,alpha=alpha[i])+
      geom_line(cluster_mean2,mapping=aes(time,effect),color=darken(cols2),size=1.5)+
      theme_bw()+theme(panel.grid =element_blank()) + 
      scale_y_continuous(limits = c(0,12),breaks=seq(0,10,4))+
      scale_x_continuous(limits = c(0,14),breaks=seq(0,14,3))+
      xlab('')+ylab('')+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))
    return(p)
  }
  
  p16 <- ggplot()+geom_line(BIC_df,mapping=aes(k,BIC),size=1.1)+ 
    theme_bw()+theme(panel.grid =element_blank()) +scale_x_continuous(limit=c(2,20),breaks=seq(0,20,5))+
    theme(legend.position = "none")+xlab("No. of Modules")+ylab("BIC")+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    geom_vline(aes(xintercept=15), colour="#F8766D", linetype="dashed",size=1)
  
  p_up <- function(i){
    p <- get_p(i)+theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.length.x = unit(-0.1,"cm"))
    return(p)
  }
  
  p_middle <- function(i){
    p <- get_p(i)+theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.length.x = unit(-0.1,"cm"),
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.length.y = unit(-0.1,"cm"))
    return(p)
  }
  
  p_down <- function(i){
    p <- get_p(i)+theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.length.y = unit(-0.1,"cm"))
    return(p)
  }
  p <- list()
  for (i in c(1,5,9)) {
    p[[i]] <- p_up(i)
  }
  for (i in c(2,3,4,6,7,8,10,11,12)) {
    p[[i]] <- p_middle(i)
  }
  for (i in c(14:15)) {
    p[[i]] <- p_down(i)
  }
  p[[13]] <- get_p(13)
  p0 <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
                  p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p16,ncol = 4,nrow = 4)+ 
    theme(plot.margin = margin(5, 5, 5, 30))
  p0 <- add_sub(p0, "Time (day)", hjust = 0.3,vjust=-0.5)
  p0 <- add_sub(p0, "Generic Effect",x=-0.05, y=3.5, angle = 90,vjust=3)
  ggsave("Figure2_B.pdf",ggdraw(p0),width = 10, height = 8,units = "in")
  cat('Done')
}

get_figure3a <- function(k,times,all_net,cluster_result){
  requiredPackages = c("ggplot2","reshape2","orthopolynom",'cowplot','patchwork')
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  cat('Plotting Figure3A','\n')
  legendre_fit <- function(par){
    x <- times
    fit <- sapply(1:length(par),function(c)
      par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
    legendre_fit <- as.matrix(as.data.frame(polynomial.values(
      polynomials=fit,x=scaleX(x, u=-1, v=1))))
    x_interpolation <- rowSums(legendre_fit)
    return(x_interpolation)
  }
  get_curve_data <- function(i){
    times <- times
    cluster_mean <- t(sapply(1:k, function(c)legendre_fit(cluster_result$curve_par[c,1:5])))
    rownames(cluster_mean) <- 1:k
    d1 <- data.frame(all_net[[i]][[6]],check.names = F)
    colnames(d1)[-1] <- all_net[[i]][[2]]
    d1$sum <- rowSums(d1)
    d1$origin <- as.numeric(cluster_mean[i,])
    d1$x <- times
    return(d1)
  }
  microbe <- lapply(1:k,function(c)get_curve_data(c))
  col <- c('#ff7171','#0dceda','#9ede73')
  darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
  }
  
  myplot <- function(i){
    d <- microbe[[i]]
    p1 <- ggplot(d,aes(x=x))+geom_line(d,mapping=aes(x=x,y=origin),color=col[2],size=1.5)+
      geom_line(d,mapping=aes(x=x,y=ind_effect),color=col[1],size=1.4)
    for (j in 1:(ncol(d)-4)) {
      p1 <- p1+geom_line(aes_string(y=d[,1+j]),color=col[3],size=1.2,alpha=1)
    }
    p1 <- p1+ theme(panel.grid = element_blank(), 
                    panel.background = element_rect(color = 'black', fill = 'transparent'), 
                    legend.title = element_blank()) +
      xlab(NULL) + ylab(NULL)+
      theme(plot.margin = unit(c(0,0,0,0),"lines"))+
      annotate('text',x=7,y=round(max(d[,-ncol(d)])),label=paste0('M',1:k)[i],size=5)+
      geom_hline(yintercept = 0, size = 0.6)
    return(p1)
  }
  plot0 <- function(i){
    tmp1 <- myplot(i)+theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.length.x = unit(-0.1,"cm"))
    p0 <- tmp1
    return(p0)
  }
  p <- list()
  for (i in 1:10) {
    p[[i]] <- plot0(i)
  }
  p[[11]] <- myplot(11)
  p[[12]] <- myplot(12)
  p[[13]] <- myplot(13)
  p[[14]] <- myplot(14)
  p[[15]] <- myplot(15)
  
  pp <- p[[1]]+p[[2]]+p[[3]]+p[[4]]+p[[5]]+p[[6]]+p[[7]]+p[[8]]+p[[9]]+p[[10]]+p[[11]]+
    p[[12]]+p[[13]]+p[[14]]+p[[15]]+ plot_layout(ncol = 5)
  pp <- add_sub(pp, "Time (day)", hjust = 0.3,vjust=-0.5)
  ggsave("Figure3_A.pdf",pp,width =15, height = 8,units = "in")
  cat('Done')
}

get_figure3c <- function(k,times,all_net){
  requiredPackages = c("igraph","pbapply")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  cat('Plotting Figure3C','\n')
  get_after <- function(i){
    temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
    temp[,1] <- i[[2]]
    temp[,2] <- i[[1]]
    temp[,3] <- i[[5]][2:(length(i[[2]])+1)]
    
    colnames(temp) <- c('from','to','dep_effect')
    temp <- data.frame(temp)
    temp[,3] <- as.numeric(as.character(temp[,3]))
    return(temp)
  }
  
  get_max_effect <- function(k){
    after <- do.call(rbind,lapply(k, get_after))
    max_dep_effect <- max(abs(after$dep_effect))
    
    temp <- aggregate(dep_effect ~ to, data = after, sum)
    all_dep_effect <- max(abs(temp$dep_effect))
    return(c(max_dep_effect,all_dep_effect))
  }
  
  max_effect <- t(sapply(1:length(all_net),function(c)get_max_effect(all_net[[c]])))
  
  network_plot <- function(k,title){
    #extra <- as.numeric(table(df$cluster))
    get_extra <- function(i){
      temp <- i[[5]][1]
      return(temp)
    }
    extra <- sapply(k,get_extra)
    
    after <- do.call(rbind,lapply(k, get_after))
    
    colfunc <- colorRampPalette(c("#619CFF", #ggplot blue
                                  "#ffdead", 
                                  "#F8766D"))#ggplot red
    #unchange-weight-data
    #links
    #colour_edge
    edge_number <- round(max(max_effect[,1]))
    edge_col <- data.frame(colfunc(2*edge_number+1),seq(-edge_number,edge_number,1))
    
    get_edge_colour <- function(i){
      temp <-  round(i)
      temp2 <- adjustcolor(edge_col[which(edge_col[,2]==temp),1], alpha=0.65)
      return(temp2)
    }
    links <- after
    colnames(links) <- c("from","to","weight")
    links$edge.colour <- pbsapply(links$weight,get_edge_colour) #add colour for links
    
    #nodes
    node_number <- round(max(max_effect[,2]))
    node_col <- data.frame(colfunc(2*node_number+1),seq(-node_number,node_number,1))
    get_vertex_colour <- function(i){
      temp <-  round(i)
      temp2 <- adjustcolor(node_col[which(node_col[,2]==temp),1], alpha=1)
      return(temp2)
    }
    
    nodes <- data.frame(unique(links[,2]),paste0("SNP",1:135),extra)
    colnames(nodes) <- c("id","name","ind_effect")
    nodes <- nodes[order(nodes[,1]),]
    nodes$influence <- aggregate(weight ~ to, data = links, sum)[,2]
    nodes$colour <- pbsapply(nodes$influence,get_vertex_colour) #add colour for links
    
    #normalization
    normalization <- function(x){(x-min(x))/(max(x)-min(x))*1.5+0.3} 
    
    #final plot
    links[,3] <- normalization(abs(links[,3]))
    nodes[,3:4] <- normalization(abs(nodes[,3:4]))
    net <- graph_from_data_frame( d=links,vertices = nodes,directed = T ) 
    
    #layout
    set.seed(1)
    l <- layout_randomly(net)
    
    plot.igraph(net,
                 vertex.label=V(net)$name,
                 vertex.label.color="black",
                 vertex.shape="circle", 
                 vertex.label.cex=V(net)$ind_effect*0.5,
                 vertex.size=V(net)$ind_effect*10,
                 edge.curved=0.05,
                 edge.color=E(net)$edge.colour,
                 edge.frame.color=E(net)$edge.colour,
                 edge.width=E(net)$weight*4,
                 vertex.color=V(net)$colour,
                 layout=l,
                 main=title,
                 margin=c(-.05,-.05,-.05,-.05)
    )
  }
  pdf("Fig3_C_net.pdf",width=15,height=8)
  layout(matrix(1:2, 1, 2, byrow = TRUE))
  network_plot(all_net[[1]],"Control")
  network_plot(all_net[[2]],"Stress")
  dev.off()
  cat('Done')
}