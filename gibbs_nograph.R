Rcpp::sourceCpp("MCMC/prob_ind_new.cpp")
Rcpp::sourceCpp("MCMC/norm_weights.cpp")
Rcpp::sourceCpp("MCMC/marg_cpp.cpp")
source("MCMC/move_decomposable.R")


library(abind)
library(plyr)
library(prodlim)

Gibbs_nograph = function(X, S, burn, a_alpha, b_alpha, a){
  
  ###########
  ## INPUT ##
  ###########
  
  # X           : (n, q) data matrix (nb this is a numeric matrix to use the rcpp function)
  # S           : number of MCMC iterations
  # burn        : the burn-in period
  # a           : the common hyperparameter of the Hyper Dirichlet prior
  # a_alpha     : hyperparameter of the concentration parameter alpha0 of the DP 
  # b_alpha     : hyperparameter of the concentration parameter alpha0 of the DP
 
  
  ############
  ## OUTPUT ##
  ############
  
  # Graphs_chain : a list of length S, where each element is a (q,q,K) array containing the K graphs
  # Xi_chain     : an (n,S) matrix containing n cluster indicators for each of the S MCMC iterations
  # alpha0_chain : an (S,1) vecor collecting posterior draws of the concentration parameter alpha_0
  # simil_mat    : an (n,n) posterior similarity matrix
  # graph_probs  : a (q,q,n) array with n (q,q) matrices each collecting subject-specific PPIs
  
  
  n = nrow(X)
  q = ncol(X)
  colnames(X) = as.character(paste0("X",1:q))
  
  I.cal  = sapply(1:ncol(X), function(j) length(unique(X[,j])))
  I.size = prod(I.cal)
  
  Xi_chain = matrix(NA, n, S)
  Graphs_chain = vector(mode = "list", length = S)
  
  alpha0_chain = rep(NA,S)
  
  S_base    = S
  burn_base = burn
  
  ## Sample from the baseline over graphs
  
  #out_baseline = sample_baseline_dec(S = S_base, burn = burn_base,
  #                                   q = q, a_pi = a_pi, b_pi = b_pi)$G_chain
  
  
  ########################
  ## Set initial values ##
  ########################
  
  ## Number of clusters and concentration parameter
  
  K_0 = 2
  alpha0 = rgamma(1, a_alpha, b_alpha); alpha0_chain[1] = alpha0
  
  ## Graphs
  
  G_0 = array(0, c(q, q, K_0))
  
  
  Graphs_chain[[1]] = G_0
  
  ## Cluster indicators
  
  xi = sample(K_0, n, replace = TRUE)
  
  while(length(table(xi)) < K_0){
    xi = sample(K_0, n, replace = TRUE)
  }
  
  Xi_chain[,1] = xi
  
  Graphs = G_0
  r = table(xi)
  
  
  #####################
  ## MCMC iterations ##
  #####################
  
  pb = txtProgressBar(min = 2, max = S, style = 3)
  
  for(t in 2:S){
    
    
    ###############################
    ## Update cluster indicators ##
    ###############################
    
    
    G_star = matrix(0, q, q)
    Graphs = abind(Graphs, G_star)
    K_star = dim(Graphs)[3]
    
    logProbs = matrix(nrow = n, ncol = K_star) 
    
    for(k in 1:(K_star-1)){
      
      cliques_separators_k = mpd(Graphs[,,k])
      
      cl_k = lapply(1:length(cliques_separators_k$cliques), function(i) as.integer(cliques_separators_k$cliques[[i]]))
      se_k = lapply(1:length(cliques_separators_k$separators), function(i) as.integer(cliques_separators_k$separators[[i]]))
      

      logProbs[,k] = sapply(1:n, function(i) prob_ik_nonempty_test(Xk = X[xi == k,, drop = FALSE], cl_list = cl_k,
                                                                   se_list = se_k, 
                                                                   member = (xi[i] == k), x = X[i,], a = a, I_cal = I.cal))
      
    }
    
    cliques_separators_kstar = mpd(Graphs[,,K_star])
    
    cl_kstar = lapply(1:length(cliques_separators_kstar$cliques), function(i) as.integer(cliques_separators_kstar$cliques[[i]]))
    se_kstar = lapply(1:length(cliques_separators_kstar$separators), function(i) as.integer(cliques_separators_kstar$separators[[i]]))
    
    #logProbs[,K_star] = prob_ik_empty(cl.list = cl_kstar, se.list = se_kstar, a = a, alpha0 = alpha0, I.cal = I.cal)
    logProbs[,K_star] = prob_ik_empty_test(cl_list = cl_kstar, se_list = se_kstar, a = a, alpha0 = alpha0, I_cal = I.cal)
    probs = t(sapply(1:nrow(logProbs), function(i) normalize_weights_test(logProbs[i,])))
    
    xi_star = sapply(1:n, function(i) sample(1:(K_star), size = 1, prob = probs[i,]))
    
    labs = as.integer(names(table(xi_star)))
    K_star = length(labs) # effective new number of clusters
    
    Graphs  = array(Graphs[,,labs], c(q,q,K_star))
    xi_star = as.factor(xi_star); levels(xi_star) = 1:K_star # update labels
    
    r = table(xi_star)
    xi = c(xi_star)
    K = dim(Graphs)[3] # number of non-empty clusters
    
    
    ####################################
    ## Update concentration parameter ##
    ####################################
    
    eta = rbeta(1, alpha0 + 1, n)
    
    alpha0 = c(rgamma(1, shape = a_alpha + K, 
                      rate = b_alpha - log(eta)), 
               rgamma(1, shape = a_alpha + K - 1,
                      rate = b_alpha - log(eta)))[sample(c(1,2), 1, prob = c(a_alpha + K - 1, n*(b_alpha - log(eta))))]
    
    alpha0_chain[t] = alpha0
    
    Xi_chain[,t]      = xi
    
    setTxtProgressBar(pb, t)
    close(pb)
    
  }
  
  
  #########################
  ## Posterior summaries ##
  #########################
  
  ## Posterior similarity matrix
  
  simil_mat = matrix(0, nrow = n, ncol = n)
  
  for(t in (burn + 1):S){
    
    simil_mat = simil_mat + (matrix(Xi_chain[,t], nrow = n, ncol = n) == t(matrix(Xi_chain[,t], nrow = n, ncol = n)))*1
    
  }
  
  simil_probs = simil_mat/(S - burn)
  
  
  return(list(Xi_chain = Xi_chain, alpha0_chain = alpha0_chain,
              simil_probs = simil_probs, graph_probs = NULL))
  
  
}



