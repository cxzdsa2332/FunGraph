get_LR_effect <- function(marker,times,data){
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
  
  #plot(c(get_miu(NH_0$par[1:5]),get_miu(NH_0$par[6:10])))
  #equation_4
  get_genetic_effect <- function(times){
    get_miu <- function(miu_par,times){
      t = times
      miu <- miu_par[1]/(1 + miu_par[2] * exp(-miu_par[3] * t)) - (miu_par[4] * exp(-miu_par[5] * t))
      miu
    }
    get_mu2 <- function(par){ mu <- c(get_miu(par[1:5],times),get_miu(par[6:10],times))}
    
    effect1 <- nrow(pheno_0)*(get_mu2(NH_1$par[1:10])^2) +
      nrow(pheno_1)*(get_mu2(NH_1$par[11:20])^2)
    effect2 <- nrow(pheno_0)*get_mu2(NH_1$par[1:10]) +
      nrow(pheno_1)*get_mu2(NH_1$par[11:20])
    if (nrow(pheno_2)!=0) {
      effect1 <- effect1+nrow(pheno_2)*(get_mu2(NH_1$par[21:30])^2)
      effect2 <- effect2+nrow(pheno_2)*get_mu2(NH_1$par[21:30])
    }
    effect <- sqrt( effect1/nrow(y_all) - (effect2/nrow(y_all))^2 )
    return(effect)
  }
  #plot(get_genetic_effect(seq(1,14,0.1)))
  generic_effect <- get_genetic_effect(times)
  #result
  return_object <- list(LR,generic_effect,NH_0$par,NH_1$par)
  return(return_object)
}

get_FunClu_rewsult <- function(data){
  geno_df2 <- data
  core.number <- detectCores()
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(mvtnorm)})
  clusterExport(cl, c("get_LR_effect","geno_df2","pheno_df"),envir=environment())
  result1 <- pblapply(cl=cl,1:nrow(geno_df2),function(c)get_LR_effect(c,1:14,pheno_df))
  stopCluster(cl)
  
  LR_result <- do.call(c,lapply(1:length(result1),function(c) result1[[c]][[1]]))
  generic_effect <- do.call(rbind,lapply(1:length(result1),function(c) result1[[c]][[2]]))
  rownames(generic_effect) <- geno_df[,3]
  colnames(generic_effect) <- colnames(pheno_df)[-1]
  return_object <- list(LR_result,generic_effect,result1)
  names(return_object)<-c("LR_results", "generic_effect","FunMap_pars")
  return(return_object)
}

get_init_par <- function(data,k,legendre_order){
  get_legendre_par <- function(y,legendre_order,x){
    #lm_method
    get_legendre_matrix <- function(x,legendre_order){
      legendre_coef <- legendre.polynomials(n=legendre_order, normalized=F)
      legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
        polynomials=legendre_coef,x=scaleX(x, u=-1, v=1))))
      colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
      return(legendre_matrix[,2:(legendre_order+1)])
    }
    legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
    return(legendre_par)
  }
  #get inital pars based on k-means
  init_cluster <- kmeans(data,centers = k,iter.max = 1000)
  cuM <- init_cluster$centers
  
  init_curve_par <- cbind(t(sapply(1:k,function(c)get_legendre_par(cuM[c,1:14],legendre_order,1:14))),
                          t(sapply(1:k,function(c)get_legendre_par(cuM[c,15:28],legendre_order,1:14))))
  init_SAD_par <- c(1.06,0.25,1.15,0.18)
  init_pro <- table(init_cluster$cluster)/nrow(data)
  return_object <- list(init_SAD_par,init_curve_par,init_pro)
  names(return_object)<-c("init_SAD_par","init_curve_par","init_pro")
  return(return_object)
}

get_cluster <- function(data,k,input,legendre_order){
  Delta <- 100; iter <- 0; itermax <- 100;
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
  legendre_fit <- function(par){
    x <- 1:(ncol(data)/2)
    fit <- sapply(1:length(par),function(c)
      par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
    legendre_fit <- as.matrix(as.data.frame(polynomial.values(
      polynomials=fit,x=scaleX(x, u=-1, v=1))))
    x_interpolation <- rowSums(legendre_fit)
    return(x_interpolation)
  }
  get_biSAD1 <- function(par){
    n=ncol(data)/2
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
  mle <- function(par,data,prob){
    par1 <- par[1:4]
    par2 <- matrix(par[-c(1:4)],nrow = k,ncol = (legendre_order+1)*2)
    miu <- t( sapply(1:k, function(c)c(legendre_fit(par2[c,1:(legendre_order+1)]),
                                       legendre_fit(par2[c,(ncol(par2)/2+1):ncol(par2)]))) )
    temp_S <- sapply(1:k,function(c)dmvnorm(data,
                                            miu[c,],
                                            get_biSAD1(par1))*prob[c] )
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  
  while ( Delta > 1 && iter <= itermax ) {
    # initiation
    if(iter == 0){
      init_SAD_par <- input[[1]]
      init_curve_par <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_SAD_par,init_curve_par)
    LL_mem <- mle(old_par,data,pro)
    miu <- t( sapply(1:k, function(c)c(legendre_fit(init_curve_par[c,1:(legendre_order+1)]),
                                       legendre_fit(init_curve_par[c,(ncol(init_curve_par)/2+1):ncol(init_curve_par)]))) )
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             miu[c,],
                                             get_biSAD1(init_SAD_par))*pro[c] )
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
    cat('\n',"iter=",iter,"LL=",L_Value,'\n')
    iter <- iter+1; LL_mem <- L_Value
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  #plot-----------
  cluster <- apply(omega,1,which.max)
  clustered_ck <- data.frame(row.names(data),data[,1:14],cluster)
  clustered_salt <- data.frame(row.names(data),data[,15:28],cluster)
  clustered_data <- clustered_ck
  get_plot <- function(clustered_data){
    colnames(clustered_data) <- c("marker",1:14,"cluster")
    long_df <- melt(clustered_data,c("marker","cluster"))
    colnames(long_df) <- c("marker","cluster","time","effect")
    p <-  ggplot()+geom_line(long_df,mapping=aes(as.numeric(as.character(time)),effect,group=marker,
                                                 colour= as.character(cluster)),alpha=1)+
      facet_wrap(long_df$cluster,scales = "fixed")+ 
      theme(legend.position="none") + xlab("Time")+ylab("generic_effect")
    return(p)
  }
  p1 <- get_plot(clustered_ck); p2 <- get_plot(clustered_salt)
  clustered_ck <- clustered_ck[,-1];clustered_salt <- clustered_salt[,-1]
  return_object <- list(init_SAD_par,init_curve_par,pro,LL_mem,BIC,clustered_ck,clustered_salt,p1,p2)
  names(return_object)<-c("SAD_par", "curve_par", "pro", "LL", 
                          "BIC", "clustered_ck","clustered_salt","plot1","plot2")
  return(return_object)
}

get_module_result <- function(k,times,order){
  legendre_fit <- function(par){
    x <- times
    fit <- sapply(1:length(par),function(c)
      par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
    legendre_fit <- as.matrix(as.data.frame(polynomial.values(
      polynomials=fit,x=scaleX(x, u=-1, v=1))))
    x_interpolation <- rowSums(legendre_fit)
    return(x_interpolation)
  }
  get_legendre_par <- function(times,order,par) {
    get_interaction <- function(data,col){
      n <- nrow(data)
      clean_data <- data
      gene_list <- list()
      m <- clean_data[,col]
      M <- clean_data[,-col]
      x_matrix <- M
      x_matrix <- as.matrix(x_matrix)
      #vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
      #x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
      #x_matrix <- as.matrix(x_matrix)
      name <- colnames(clean_data)
      ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", 
                             family="gaussian",nfold = 10,alpha = 0)
      best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
      
      fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", family="gaussian",
                           nfold = 10,alpha = 1,
                           penalty.factor = 1/best_ridge_coef,
                           keep = TRUE)
      best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
      
      gene_list_one <- list()
      gene_list_one[[1]] <- name[col]
      gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
      gene_list_one[[3]] <- best_alasso_coef1@x[-1]
      gene_list[[col]] <- gene_list_one
      
      return(gene_list_one)
    }
    cluster_mean <- t(sapply(1:k, function(c)legendre_fit(par[c,])))
    rownames(cluster_mean) <- 1:k
    module_relationship <- pblapply(1:k,function(c) get_interaction(t(cluster_mean),c))
    #----------------------
    get_effect <- function(pars,effect,times,order,y0){
      if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
      LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
      d_LOP_fit <-  sapply(1:length(pars),function(c)
        pars[c]*polynomial.derivatives(LOP)[[c+1]])
      h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
      #rk4 for legendre with step=h
      LOP_rk4 <- function(x0,y0){
        f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
        k1 <- f(x0,y0) 
        k2 <- f(x0+h/2,y0+h/2*k1)
        k3 <- f(x0+h/2,y0+h/2*k2)
        k4 <- f(x0+h,y0+h*k3)
        y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
        return(y)
      }
      #dy_LOP, the increasment of each step
      dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],y0))
      #dy_LOP*y= main or sub effect
      dy_fit <- effect*c(0,dy[1:(length(times)-1)])
      return(cumsum(dy_fit))
    }
    
    ode_optimize <- function(pars,ind,dep,times,data,order,effect){
      ind_pars <- matrix(pars,ncol=order)[1,]
      dep_pars <- matrix(pars,ncol=order)[-1,]
      inital_value <- data[,ind][1]
      ind_effect <- get_effect(ind_pars,data[,ind],times,order,inital_value)+inital_value 
      if ( is.null(nrow(dep_pars)) ) {
        dep_effect <- get_effect(dep_pars,data[,dep],times,order,0)
        y <- ind_effect+dep_effect
      }else{
        dep_effect <- sapply(1:length(dep), function(c)
          get_effect(dep_pars[c,],data[,dep[c]],times,order,0))
        y <- ind_effect+rowSums(dep_effect)
      }
      ssr <- sum((data[,ind]-y)^2)
      #add penalty
      alpha=5e-5
      ridge <- sum((data[,ind]-y)^2+alpha*(sum(ind_pars^2)+sum(dep_pars^2)))
      return(ridge)
    }
    
    get_value <- function(effect,data,times,order){
      #input
      ind <- data[[1]]
      dep <- data[[2]]
      ind_no <- as.numeric(which(colnames(effect)==ind))
      dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
      init_pars <- rep(0.001,(length(ind_no)+length(dep_no))*order)
      result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,
                      times=times,data=effect,order=order,
                      method = "BFGS", control=list(maxit=500,trace=T))
      par_after <- matrix(result$par,length(ind)+length(dep),order)
      return(par_after)
    }
    
    core.number <- detectCores()
    cl <- makeCluster(getOption("cl.cores", core.number))
    clusterEvalQ(cl, {library(orthopolynom)})
    clusterExport(cl, c("get_value","ode_optimize","get_effect","times","get_interaction",
                        "cluster_mean","module_relationship","order"),envir=environment())
    lop_par <- pblapply(1:nrow(cluster_mean),function(c)get_value(t(cluster_mean),
                                                                  module_relationship[[c]],times,order),cl=cl)
    stopCluster(cl)
    return(list(lop_par,module_relationship))
  }
  
  all_lop_par <- lapply(1:2,function(c)get_legendre_par(times=times,order=order,
                                                        par=cluster_result$curve_par[,list(1:5,6:10)[[c]]]))
  #output for result-
  get_output <- function(relationship,par,effect,times,order){
    get_effect <- function(pars,effect,times,order,y0){
      if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
      LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
      d_LOP_fit <-  sapply(1:length(pars),function(c)
        pars[c]*polynomial.derivatives(LOP)[[c+1]])
      h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
      #rk4 for legendre with step=h
      LOP_rk4 <- function(x0,y0){
        f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
        k1 <- f(x0,y0) 
        k2 <- f(x0+h/2,y0+h/2*k1)
        k3 <- f(x0+h/2,y0+h/2*k2)
        k4 <- f(x0+h,y0+h*k3)
        y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
        return(y)
      }
      #dy_LOP, the increasment of each step
      dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],y0))
      #dy_LOP*y= main or sub effect
      dy_fit <- effect*c(0,dy[1:(length(times)-1)])
      return(cumsum(dy_fit))
    }
    output <- list()
    output[[1]] <- relationship[[1]]  
    output[[2]] <- relationship[[2]]  
    output[[3]] <- par[1,]
    output[[4]] <- par[2:nrow(par),]
    ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
    dep_no <- as.numeric(sapply(1:length(output[[2]]), 
                                function(c) which(colnames(effect)==output[[2]][c])))
    inital_value <- effect[,ind_no][1]
    ind_effect <- get_effect(as.numeric(output[[3]]),effect[,ind_no],times,order,inital_value)+inital_value
    if (length(dep_no)==1) {
      dep_effect <- get_effect(as.numeric(output[[4]]),effect[,dep_no],times,order,0)
    }else{
      dep_effect <- sapply(1:length(dep_no), function(c)
        get_effect(as.numeric(output[[4]][c,]),effect[,dep_no[c]],times,order,0))
      colnames(dep_effect) <- dep_no
    }
    #------------
    all_effect <- cbind(ind_effect,dep_effect)
    effect_mean <- apply(all_effect,2,mean)
    output[[5]] <- effect_mean
    output[[6]] <- all_effect
    return(output)
  }
  
  get_net <- function(i){
    cluster_mean <- t(sapply(1:k, function(c)legendre_fit(cluster_result$curve_par[c,list(1:5,6:10)[[i]]])))
    rownames(cluster_mean) <- 1:k
    module_relationship <- all_lop_par[[i]][[2]]
    net <- pblapply(1:k,function(c)
      get_output(module_relationship[[c]],all_lop_par[[i]][[1]][[c]],
                 t(cluster_mean),times=times,order=order))
    return(net)
  }
  
  all_net <- pblapply(1:2,function(c)get_net(c))
  j=all_net[[1]]
  get_net_output <- function(j){
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
    links <- do.call(rbind,lapply(j, get_after))
    links$from <- paste0('M',links$from)
    links$to <- paste0('M',links$to)
    get_link_color <- function(i){
      tmp <- links$dep_effect[i]
      if (tmp >= 0 ) {
        tmp2 <- '+'
      } else {
        tmp2 <- '-'
      }
      return(tmp2)
    }
    links$effect_type <- sapply(1:nrow(links),function(c)get_link_color(c))
    links$dep_effect <- abs(links$dep_effect)
    get_ind <- function(i){
      temp <- i[[5]][1]
      return(temp)
    }
    nodes <- data.frame(unique(links[,2]),paste0('M',1:k),sapply(j,get_ind))
    colnames(nodes) <- c("id","name","received_effect")
    nodes$size <- table(cluster_result$clustered_ck$cluster)
    return(list(links,nodes[,c(2,4)]))
  }
  ODE_result <- lapply(1:2,function(c)get_net_output(all_net[[c]]))
  return(list(ODE_result,all_net))
}

get_submodule_result <- function(cluster,times,order){
  df <- list(cluster_result$clustered_ck[cluster_result$clustered_ck$cluster==cluster,-15],
             cluster_result$clustered_salt[cluster_result$clustered_salt$cluster==cluster,-15])
  get_legendre_par <- function(times,order,i) {
    get_interaction <- function(data,col){
      n <- nrow(data)
      clean_data <- data
      gene_list <- list()
      m <- clean_data[,col]
      M <- clean_data[,-col]
      x_matrix <- M
      x_matrix <- as.matrix(x_matrix)
      #vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
      #x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
      #x_matrix <- as.matrix(x_matrix)
      name <- colnames(clean_data)
      ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", 
                             family="gaussian",nfold = 10,alpha = 0)
      best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
      
      fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse", family="gaussian",
                           nfold = 10,alpha = 1,
                           penalty.factor = 1/best_ridge_coef,
                           keep = TRUE)
      best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
      
      gene_list_one <- list()
      gene_list_one[[1]] <- name[col]
      gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
      gene_list_one[[3]] <- best_alasso_coef1@x[-1]
      gene_list[[col]] <- gene_list_one
      
      return(gene_list_one)
    }
    cluster_mean <- df[[i]]
    k=nrow(cluster_mean)
    module_relationship <- pblapply(1:k,function(c) get_interaction(t(cluster_mean),c))
    #----------------------
    get_effect <- function(pars,effect,times,order,y0){
      if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
      LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
      d_LOP_fit <-  sapply(1:length(pars),function(c)
        pars[c]*polynomial.derivatives(LOP)[[c+1]])
      h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
      #rk4 for legendre with step=h
      LOP_rk4 <- function(x0,y0){
        f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
        k1 <- f(x0,y0) 
        k2 <- f(x0+h/2,y0+h/2*k1)
        k3 <- f(x0+h/2,y0+h/2*k2)
        k4 <- f(x0+h,y0+h*k3)
        y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
        return(y)
      }
      #dy_LOP, the increasment of each step
      dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],y0))
      #dy_LOP*y= main or sub effect
      dy_fit <- effect*c(0,dy[1:(length(times)-1)])
      return(cumsum(dy_fit))
    }
    
    ode_optimize <- function(pars,ind,dep,times,data,order,effect){
      ind_pars <- matrix(pars,ncol=order)[1,]
      dep_pars <- matrix(pars,ncol=order)[-1,]
      inital_value <- data[,ind][1]
      ind_effect <- get_effect(ind_pars,data[,ind],times,order,inital_value)+inital_value 
      if ( is.null(nrow(dep_pars)) ) {
        dep_effect <- get_effect(dep_pars,data[,dep],times,order,0)
        y <- ind_effect+dep_effect
      }else{
        dep_effect <- sapply(1:length(dep), function(c)
          get_effect(dep_pars[c,],data[,dep[c]],times,order,0))
        y <- ind_effect+rowSums(dep_effect)
      }
      ssr <- sum((data[,ind]-y)^2)
      #add penalty
      alpha=5e-5
      ridge <- sum((data[,ind]-y)^2+alpha*(sum(ind_pars^2)+sum(dep_pars^2)))
      return(ridge)
    }
    
    get_value <- function(effect,data,times,order){
      #input
      ind <- data[[1]]
      dep <- data[[2]]
      ind_no <- as.numeric(which(colnames(effect)==ind))
      dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
      init_pars <- rep(0.001,(length(ind_no)+length(dep_no))*order)
      result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,
                      times=times,data=effect,order=order,
                      method = "BFGS", control=list(maxit=500,trace=T))
      par_after <- matrix(result$par,length(ind)+length(dep),order)
      return(par_after)
    }
    
    core.number <- detectCores()
    cl <- makeCluster(getOption("cl.cores", core.number))
    clusterEvalQ(cl, {library(orthopolynom)})
    clusterExport(cl, c("get_value","ode_optimize","get_effect","times","get_interaction",
                        "cluster_mean","module_relationship","order"),envir=environment())
    lop_par <- pblapply(1:nrow(cluster_mean),function(c)get_value(t(cluster_mean),
                                                                  module_relationship[[c]],times,order),cl=cl)
    stopCluster(cl)
    return(list(lop_par,module_relationship))
  }
  
  all_lop_par <- lapply(1:2,function(c)get_legendre_par(times=times,order=order,c))
  
  #output for result-
  get_output <- function(relationship,par,effect,times,order){
    get_effect <- function(pars,effect,times,order,y0){
      if ( length(pars) != order ) {warning("legendre_pars != legendre_order")}
      LOP <-  legendre.polynomials(order, normalized=F) #Legendre polynomials
      d_LOP_fit <-  sapply(1:length(pars),function(c)
        pars[c]*polynomial.derivatives(LOP)[[c+1]])
      h <- scaleX(times,u=-1,v=1)[2]-scaleX(times,u=-1,v=1)[1] #per step h
      #rk4 for legendre with step=h
      LOP_rk4 <- function(x0,y0){
        f <- function(x,y){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
        k1 <- f(x0,y0) 
        k2 <- f(x0+h/2,y0+h/2*k1)
        k3 <- f(x0+h/2,y0+h/2*k2)
        k4 <- f(x0+h,y0+h*k3)
        y <- y0+h/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
        return(y)
      }
      #dy_LOP, the increasment of each step
      dy <- sapply(1:length(times),function(c)LOP_rk4(scaleX(times,u=-1,v=1)[c],y0))
      #dy_LOP*y= main or sub effect
      dy_fit <- effect*c(0,dy[1:(length(times)-1)])
      return(cumsum(dy_fit))
    }
    output <- list()
    output[[1]] <- relationship[[1]]  
    output[[2]] <- relationship[[2]]  
    output[[3]] <- par[1,]
    output[[4]] <- par[2:nrow(par),]
    ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
    dep_no <- as.numeric(sapply(1:length(output[[2]]), 
                                function(c) which(colnames(effect)==output[[2]][c])))
    inital_value <- effect[,ind_no][1]
    ind_effect <- get_effect(as.numeric(output[[3]]),effect[,ind_no],times,order,inital_value)+inital_value
    if (length(dep_no)==1) {
      dep_effect <- get_effect(as.numeric(output[[4]]),effect[,dep_no],times,order,0)
    }else{
      dep_effect <- sapply(1:length(dep_no), function(c)
        get_effect(as.numeric(output[[4]][c,]),effect[,dep_no[c]],times,order,0))
      colnames(dep_effect) <- dep_no
    }
    #------------
    all_effect <- cbind(ind_effect,dep_effect)
    #effect_mean <- all_effect[5,]
    effect_mean <- apply(all_effect,2,mean)
    output[[5]] <- effect_mean
    output[[6]] <- all_effect
    return(output)
  }
  
  get_net <- function(i){
    cluster_mean <- df[[i]]
    module_relationship <- all_lop_par[[i]][[2]]
    net <- pblapply(1:nrow(cluster_mean),function(c)
      get_output(module_relationship[[c]],all_lop_par[[i]][[1]][[c]],
                 t(cluster_mean),times=times,order=order))
    return(net)
  }
  
  all_net <- pblapply(1:2,function(c)get_net(c))
  get_net_output <- function(j){
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
    links <- do.call(rbind,lapply(j, get_after))
    get_link_color <- function(i){
      tmp <- links$dep_effect[i]
      if (tmp >= 0 ) {
        tmp2 <- '+'
      } else {
        tmp2 <- '-'
      }
      return(tmp2)
    }
    links$effect_type <- sapply(1:nrow(links),function(c)get_link_color(c))
    get_ind <- function(i){
      temp <- i[[5]][1]
      return(temp)
    }
    nodes <- data.frame(unique(links[,2]),paste0('SNP',1:length(unique(links[,2]))),sapply(j,get_ind))
    colnames(nodes) <- c("id","name","received_effect")
    nodes$influence <- aggregate(dep_effect ~ to, data = links, sum)[,2]
    return(list(links,nodes))
  }
  ODE_result <- lapply(1:2,function(c)get_net_output(all_net[[c]]))
  return(list(ODE_result,all_net))
}
