# -----------------------------------rand_index--------------------------------------
randfun <- function(lvel, dv) {
  lem <- abs(sapply(lvel, function(x) x-lvel))
  lem[lem > 1] <- 1
  olem <- abs(sapply(dv, function(x) x-dv))
  olem[olem > 1] <- 1
  ada <- sum(abs(lem-olem))/2
  tnp <- choose(dim(lem)[1],2)
  ada/tnp
}
# ------------------------------------------------------------------------------------

#----------------------------pmf of G-hyg distribution--------------------------------
pmf <- function(M){
  cols <- colSums(M)
  nume <- prod(as.numeric(unlist(apply(M, 1, function(v) dmultinom(v, size = NULL, rep(1/ncol(M),ncol(M)), log = FALSE)))))
  deno <- dmultinom(cols, size = NULL, rep(1/ncol(M),ncol(M)), log = FALSE)
  nume/deno  
}
#-------------------------------------------------------------------------------------

#------------------- Holm (1979) test------------------------------------------------

Holm <- function(pvalues, alpha){
  n_hypo <- length(pvalues)
  criticalValues <- sapply(1:n_hypo, function(i) alpha/(n_hypo-i+1))
  sorted <- sort(pvalues)
  rejected <- (sorted < criticalValues)
  orp <- as.vector(order(pvalues))
  
  if(sum(rejected)==n_hypo){
    output <- as.numeric(rejected)
    return(output)
  }
  else{
    FF <- which(rejected==FALSE)
    rejected[min(FF):n_hypo] <- FALSE 
    output <- as.numeric(rejected[order(orp)])
    return(output)
  }
}
#-----------------------------------------------------------------------------------

#-----------------Benjamini-Hochberg (1995) test------------------------------------

BenHoch <- function(pvalues, alpha){ 
  n_hypo <- length(pvalues)
  criticalValues <- sapply(1:n_hypo, function(i) (i/n_hypo)*alpha)
  sorted <- sort(pvalues)
  rejected <-  (sorted <= criticalValues)
  orp <- as.vector(order(pvalues))
  
  if(sum(rejected)==0){
    output <- as.numeric(rejected)
    return(output)
  }
  else{
    FT <- which(rejected==TRUE)
    rejected[1:max(FT)] <- TRUE 
    output <- as.numeric(rejected[order(orp)])
    return(output)
  }
}
#-----------------------------------------------------------------------------------

# ------------------------------------rand_index_test-----------------------------------------------------------
RItest <- function(M, labels, sizes, n_clust, randomization = TRUE, clust_alg = "knwClustNo", kmax = 2*n_clust, s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(labels) || !is.numeric(sizes) || anyNA(M) || anyNA(labels) || anyNA(sizes))
    stop("all entries of obsevations, cluster membership index of observations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  if(nrow(M)!=length(labels)||length(labels)!=sum(sizes))
    stop("Number of observations and number of membership of observations must be same")
  N <- nrow(M)
  if(clust_alg=="knwClustNo"){
    dvec <- as.numeric(gMADD(s_psi, s_h, n_clust, lb, M))
  }
  else if(clust_alg=="estClustNo"){
    # maxClNo <- 2*n_clust
    dvec_di_mat <- as.matrix(gMADD_DI(s_psi, s_h, kmax, lb, M)) 
    est_no_cl <- which.max(dvec_di_mat[ ,(N+1)])
    dvec <- dvec_di_mat[est_no_cl,1:N]
  }
  obs_ct <- unname(table(labels, dvec))
  rand_id <- randfun(labels, dvec)

  sim_ri <- rep(0, n_sts)
  for(j in 1:n_sts){
    sc_tab <- as.matrix(rctab(obs_ct))
    s_dvec <- as.numeric(c(unlist(apply(sc_tab, 1, function(v) rep(0:(ncol(sc_tab)-1), v)))))
    sim_ri[j] <- randfun(labels, s_dvec)
  }
  
  pvalue_RI <-  mean(sim_ri<=rand_id)
  # R_alpha <- as.numeric(quantile(sim_ri, alpha))
  R_alpha <- sort(sim_ri)[floor(alpha*n_sts)]
  dec_r_x <- 0
  # if(rand_id <= R_alpha ){dec_r_x <- 1}
  
  ri_gamma <- 0
 if(randomization==TRUE){
   
   if(rand_id < R_alpha){
     dec_r_x <- 1
   } 
   
   if(rand_id==R_alpha){
     u <-runif(1) 
     s0 <- mean(sim_ri==R_alpha)
     s1 <- mean(sim_ri<R_alpha)
     
     ri_gamma <- (alpha - s1)/s0     
     if(u < ri_gamma){dec_r_x <- 1}
   }
 }
 else if(randomization==FALSE){
   if(rand_id < R_alpha){
     dec_r_x <- 1
   } 
 }
  #output
  output<-list()
  if(clust_alg=="knwClustNo"){
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedRI <- rand_id
    output$RICutoff <- R_alpha
    output$randomGamma <- ri_gamma
    output$estPvalue <- pvalue_RI
    output$decisionRI <- dec_r_x
  }
  else if(clust_alg=="estClustNo"){
    output$estClustNo <- est_no_cl
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedRI <- rand_id
    output$RICutoff <- R_alpha
    output$randomGamma <- ri_gamma
    output$estPvalue <- pvalue_RI
    output$decisionRI <- dec_r_x
  }
  return(output)
}
# ------------------------------------------------------------------------------------------------

# ------------------------------Fisher_Exact_Independence_test------------------------------------------------------
FStest <- function(M, labels, sizes, n_clust, randomization = TRUE, clust_alg = "knwClustNo", kmax = 2*n_clust, s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(labels) || !is.numeric(sizes) || anyNA(M) || anyNA(labels) || anyNA(sizes))
    stop("all entries of obsevations, cluster membership index of observations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  if(nrow(M)!=length(labels)||length(labels)!=sum(sizes))
    stop("Number of observations and number of membership of observations must be same")
  N <- nrow(M)
  if(clust_alg=="knwClustNo"){
    dvec <- as.numeric(gMADD(s_psi, s_h, n_clust, lb, M))
  }
  else if(clust_alg=="estClustNo"){
    # maxClNo <- 2*n_clust
    dvec_di_mat <- as.matrix(gMADD_DI(s_psi, s_h, kmax, lb, M)) 
    est_no_cl <- which.max(dvec_di_mat[ ,(N+1)])
    dvec <- dvec_di_mat[est_no_cl,1:N]
  }
  obs_ct <- unname(table(labels, dvec))
  fpmf <- pmf(obs_ct)

  fsim_pmf <- rep(0, n_sts)
  for(j in 1:n_sts){
    fsc_tab <- as.matrix(rctab(obs_ct))
    fsim_pmf[j] <- pmf(fsc_tab)
  }
  
  pvalue_FS <- mean(fsim_pmf<=fpmf)
  # f_alpha <- as.numeric(quantile(fsim_pmf, alpha))
  f_alpha <- sort(fsim_pmf)[floor(alpha*n_sts)]
  decfisher <- 0
  # if(fpmf <= f_alpha){decfisher <- 1}
  
  fgamma <- 0
  if(randomization==TRUE){
    
    if(fpmf < f_alpha){
      decfisher <- 1
    } 
    
    if(fpmf==f_alpha){
      u1 <-runif(1) 
      fs0 <- mean(fsim_pmf==f_alpha)
      fs1 <- mean(fsim_pmf<f_alpha)
      
      fgamma <- (alpha - fs1)/fs0     
      if(u1 < fgamma){decfisher <- 1}
    }
  }
  else if(randomization==FALSE){
    if(fpmf < f_alpha){
      decfisher <- 1
    } 
  }
  #output
  output<-list()
  if(clust_alg=="knwClustNo"){
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedProb <- fpmf
    output$FCutoff <- f_alpha
    output$randomGamma <- fgamma
    output$estPvalue <- pvalue_FS
    output$decisionF <- decfisher
  }
  else if(clust_alg=="estClustNo"){
    output$estClustNo <- est_no_cl
    output$estClustLabel <- dvec
    output$obsCtyTab <- obs_ct
    output$ObservedProb <- fpmf
    output$FCutoff <- f_alpha
    output$randomGamma <- fgamma
    output$estPvalue <- pvalue_FS
    output$decisionF <- decfisher
  }
  return(output)
}

# ------------------------------------------------------------------------------------------

# ------------------------------Multiple_test_rand_index_test------------------------------------
MTRItest <- function(M, labels, sizes, k_max, multTest = "Holm", s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  n_hypo <- (k_max-1)
  ri_vec <- rep(0, n_hypo)
  p_values_vec <- rep(0, n_hypo)
  contingenyTabs <- list()
  for (i in 1:n_hypo) {
    testFun <- RItest(M, labels, sizes, n_clust = (i+1), randomization = TRUE, clust_alg = "knwClustNo", kmax = 2*(i+1), s_psi, s_h, lb, n_sts, alpha)
    ri_vec[i] <- testFun$ObservedRI
    p_values_vec[i] <- testFun$estPvalue
    contingenyTabs[[i]] <- testFun$obsCtyTab 
  }
  
  dec_mtri <- 0
  if(multTest=="Holm"){
    multest_dec <- Holm(p_values_vec,alpha)
    if(sum(multest_dec)!=0){
      dec_mtri <- 1
    }
  }
  else if(multTest=="BenHoch"){
    multest_dec <- BenHoch(p_values_vec, alpha)
    if(sum(multest_dec) != 0){
      dec_mtri <- 1
    }
  }
  
    output<-list()
    
    output$RIvec <- ri_vec
    output$Pvalues <- p_values_vec
    output$decisionMTRI <- dec_mtri
    output$contTabs <- contingenyTabs
    output$mulTestdec <- multest_dec
    return(output)
}

# ------------------------------------------------------------------------------------------

# ------------------------------Multiple_test__Fisher_exact_independent_test------------------------------------
MTFStest <- function(M, labels, sizes, k_max, multTest = "Holm", s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  n_hypo <- (k_max-1)
  fpmf_vec <- rep(0, n_hypo)
  p_values_vec <- rep(0, n_hypo)
  contingenyTabs <- list()
  for (i in 1:n_hypo) {
    testFun <- FStest(M, labels, sizes, n_clust = (i+1), randomization = TRUE, clust_alg = "knwClustNo", kmax = 2*(i+1), s_psi, s_h, lb, n_sts, alpha)
    fpmf_vec[i] <- testFun$ObservedProb
    p_values_vec[i] <- testFun$estPvalue
    contingenyTabs[[i]] <- testFun$obsCtyTab 
  }
  
  dec_mtFS <- 0
  if(multTest=="Holm"){
    multest_dec <- Holm(p_values_vec,alpha)
    if(sum(multest_dec)!=0){
      dec_mtFS <- 1
    }
  }
  else if(multTest=="BenHoch"){
    multest_dec <- BenHoch(p_values_vec, alpha)
    if(sum(multest_dec) != 0){
      dec_mtFS <- 1
    }
  }
  
  output<-list()
  
  output$fpmfvec <- fpmf_vec
  output$Pvalues <- p_values_vec
  output$decisionMTFS <- dec_mtFS
  output$contTabs <- contingenyTabs
  output$mulTestdec <- multest_dec
  return(output)
}

# ------------------------Aggregate_rand_index_test----------------------------------------
ARItest <- function(M, sizes, randomization = TRUE, clust_alg = "knwClustNo", kmax = 4, multTest = "Holm", s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(sizes) || anyNA(M) || anyNA(sizes))
    stop("all entries of obsevations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  no_clust <- length(sizes)
  cumS <- rep(0,no_clust+1)
  cumS[1] <- 0
  for(kk in 2:(no_clust+1)){
    cumS[kk] <- cumS[kk-1] + sizes[kk-1]
  }
  
  aggRI <- 1
  epvalues <- 0
  minrow <- rep(1,n_sts)
  for(i in 1:(no_clust-1)){
    mm <- sizes[i]
    for(j in (i+1):no_clust){
      nn <- sizes[j]
      level <- c(rep(0,mm), rep(1,nn))
      obsM <- M[c((cumS[i]+1):cumS[i+1],(cumS[j]+1):cumS[j+1]), ]
      testFun <- RItest(obsM, labels=level, sizes = c(mm,nn), n_clust = 2, randomization, clust_alg, kmax, s_psi, s_h, lb, n_sts, alpha) 
      aggRI <- min(aggRI,testFun$ObservedRI)
      epvalues <- c(epvalues,testFun$estPvalue)
      
      sim_ri <- rep(0, n_sts)
      for(jj in 1:n_sts){
        sc_tab <- as.matrix(rctab(testFun$obsCtyTab))
        s_dvec <- c(unlist(apply(sc_tab, 1, function(v) rep(0:(ncol(sc_tab)-1), v))))
        sim_ri[jj] = randfun(level, s_dvec)
      }
      minrow <- pmin(minrow,sim_ri)
    }
  }
  
  # AR_alpha <- as.numeric(quantile(minrow, alpha))
  AR_alpha <- sort(minrow)[floor(alpha*n_sts)]
  dec_ARI <- 0
  # if(aggRI <= AR_alpha ){dec_ARI <- 1}
  
  randGamma <- 0
  if(randomization==TRUE){
    
    if(aggRI< AR_alpha){
      dec_ARI <- 1
    } 
    
    if(aggRI==AR_alpha){
      u <-runif(1)
      s0 <- mean(minrow==AR_alpha)
      s1 <- mean(minrow<AR_alpha)
      
      randGamma <- (alpha - s1)/s0     
      if(u < randGamma){dec_ARI <- 1}
    }
  }
  else if(randomization==FALSE){
    if(aggRI < AR_alpha){
      dec_ARI <- 1
    } 
  }
  lpvs <- no_clust*(no_clust-1)/2
  pvalues <- epvalues[2:(lpvs+1)]
  hyposet <- t(combn(1:no_clust,2))
  
  if(multTest == "Holm"){
    multest_dec <- Holm(pvalues,alpha)
    multipletest <- data.frame(Population = hyposet, rejected = as.logical(multest_dec), pvalues=pvalues)
  }
  
  else if(multTest == "BenHoch"){
    multest_dec <- BenHoch(pvalues,alpha)
    multipletest <- data.frame(Population = hyposet, rejected = as.logical(multest_dec), pvalues=pvalues)
  }
  
  #output
  output<-list()
  output$ARIStat <- aggRI
  output$ARICutoff <- AR_alpha
  output$randomGamma <- randGamma
  output$decisionARI <- dec_ARI
  output$multipleTest <- multipletest
 
  # if(clust_alg=="knwClustNo"){
  #   output$ARIStat <- aggRI
  #   output$ARICutoff <- AR_alpha
  #   output$randomGamma <- randGamma
  #   output$decisionARI <- dec_ARI
  #   output$multipleTest <- multipletest
  # }
  # else if(clust_alg=="estClustNo"){
  #   output$ARIStat <- aggRI
  #   output$ARICutoff <- AR_alpha
  #   output$randomGamma <- randGamma
  #   output$decisionARI <- dec_ARI
  #   output$multipleTest <- multipletest
  #   }
  return(output)
}
  
# ------------------------------------------------------------------------------------------

# ------------------------Aggregate_Fisher_Exact_Independence_test------------------------------
AFStest <- function(M, sizes, randomization = TRUE, clust_alg = "knwClustNo", kmax = 4, multTest = "Holm", s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05){
  if(is.data.frame(M))
    M <- as.matrix(M)
  if(!is.numeric(M) || !is.numeric(sizes) || anyNA(M) || anyNA(sizes))
    stop("all entries of obsevations and number of observations from each of the population must be finite")
  if(any(sizes < 2L))
    stop("Each of population must have at least 2 obsevations")
  no_clust <- length(sizes)
  cumS <- rep(0,no_clust+1)
  cumS[1] <- 0
  for(kk in 2:(no_clust+1)){
    cumS[kk] <- cumS[kk-1] + sizes[kk-1]
  }
  
  aggAFS <- 1
  epvalues <- 0
  minrow <- rep(1,n_sts)
  
  for(i in 1:(no_clust-1)){
    mm <- sizes[i]
    for(j in (i+1):no_clust){
      nn <- sizes[j]
      level <- c(rep(0,mm), rep(1,nn))
      obsM <- M[c((cumS[i]+1):cumS[i+1],(cumS[j]+1):cumS[j+1]), ]
      testFun <- FStest(obsM, labels=level, sizes = c(mm,nn), n_clust = 2, randomization, clust_alg, kmax, s_psi, s_h, lb, n_sts, alpha) 
      aggAFS <- min(aggAFS,testFun$ObservedProb)
      epvalues <- c(epvalues,testFun$estPvalue)
      
      fsim_pmf <- rep(0, n_sts)
      for(jj in 1:n_sts){
        fsc_tab <- as.matrix(rctab(testFun$obsCtyTab))
        fsim_pmf[jj] = pmf(fsc_tab)
      }
      
      minrow <- pmin(minrow,fsim_pmf)
    }
  }
  
  # AF_alpha <- as.numeric(quantile(minrow, alpha))
  AF_alpha <- sort(minrow)[floor(alpha*n_sts)]
  dec_AFS <- 0
  # if(aggAFS <= AF_alpha ){dec_AFS <- 1}
  
  randGamma <- 0
  if(randomization==TRUE){
    
    if(aggAFS< AF_alpha){
      dec_AFS <- 1
    } 
    
    if(aggAFS==AF_alpha){
      u1 <-runif(1) 
      fs0 <- mean(minrow==AF_alpha)
      fs1 <- mean(minrow<AF_alpha)
      
      randGamma <- (alpha - fs1)/fs0     
      if(u1 < randGamma){dec_AFS <- 1}
    }
  }
  else if(randomization==FALSE){
    if(aggAFS< AF_alpha){
      dec_AFS <- 1
    } 
  }
  lpvs <- no_clust*(no_clust-1)/2
  pvalues <- epvalues[2:(lpvs+1)]
  hyposet <- t(combn(1:no_clust,2))
  
  if(multTest == "Holm"){
    multest_dec <- Holm(pvalues,alpha)
    multipletest <- data.frame(Population = hyposet, rejected = as.logical(multest_dec), pvalues=pvalues)
  }
  
  else if(multTest == "BenHoch"){
    multest_dec <- BenHoch(pvalues,alpha)
    multipletest <- data.frame(Population = hyposet, rejected = as.logical(multest_dec), pvalues=pvalues)
  }
  
  #output
  output<-list()
  output$AFSStat <- aggAFS
  output$AFCutoff <- AF_alpha
  output$randomGamma <- randGamma
  output$decisionAFS <- dec_AFS
  output$multipleTest <- multipletest
  
  # if(clust_alg=="knwClustNo"){
  #   output$AFSStat <- aggAFS
  #   output$AFCutoff <- AF_alpha
  #   output$randomGamma<- randGamma
  #   output$decisionAFS <- dec_AFS
  #   output$multipleTest <- multipletest
  # }
  # else if(clust_alg=="estClustNo"){
  #   output$AFSStat <- aggAFS
  #   output$AFCutoff <- AF_alpha
  #   output$randomGamma <- randGamma
  #   output$decisionAFS <- dec_AFS
  #   output$multipleTest <- multipletest
  #   }
  return(output)
}
# ------------------------------------------------------------------------------------------

