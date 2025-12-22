source("MCMC/sample_from_baseline.R")
Rcpp::sourceCpp("MCMC/prob_ind_new.cpp")
Rcpp::sourceCpp("MCMC/norm_weights.cpp")
Rcpp::sourceCpp("MCMC/marg_cpp.cpp")
#source("MCMC/marg_dec.R")
source("MCMC/move_decomposable.R")
#source("MCMC/prob_indicator_ik_rcpp.R")
#source("MCMC/normalize_weights.R")

library(abind)
library(plyr)
library(prodlim)

Gibbs_collapsed = function(X, S, burn, a_pi, b_pi, a_alpha, b_alpha, a){
  
  ###########
  ## INPUT ##
  ###########
  
  # X           : (n, q) data matrix (nb this is a numeric matrix to use the rcpp function)
  # S           : number of MCMC iterations
  # burn        : the burn-in period
  # a           : the common hyperparameter of the Hyper Dirichlet prior
  # a_alpha     : hyperparameter of the concentration parameter alpha0 of the DP 
  # b_alpha     : hyperparameter of the concentration parameter alpha0 of the DP
  # a_pi        : hyperparameter of the Beta prior on probability of edge inclusion
  # b_pi        : hyperparameter of the Beta prior on probability of edge inclusion
  
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
  
  out_baseline = sample_baseline_dec(S = S_base, burn = burn_base,
                                     q = q, a_pi = a_pi, b_pi = b_pi)$G_chain
  
  
  ########################
  ## Set initial values ##
  ########################
  
  ## Number of clusters and concentration parameter
  
  K_0 = 2
  alpha0 = rgamma(1, a_alpha, b_alpha); alpha0_chain[1] = alpha0
  
  ## Graphs
  
  G_0 = out_baseline[,,sample(1:(S_base - burn_base), 2)]
  
  duplicated_graphs = c(FALSE, TRUE)
  
  while(!sum(duplicated_graphs) == 0){
    
    G_0 = out_baseline[,,sample(S_base - burn_base,2)]
    
    Graphs_tmp = sapply(1:dim(G_0)[3], function(i) c(G_0[,,i]))
    duplicated_graphs = duplicated(t(Graphs_tmp))
    
  }
  
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
    
    #duplicated_graphs = c(rep(FALSE, dim(Graphs)[3]), TRUE)
    
    #while(duplicated_graphs[length(duplicated_graphs)] == TRUE){
      
    #  G_star = out_baseline[,,sample(S_base - burn_base, 1)]
    #  Graphs = abind(Graphs, G_star)
      
    #  Graphs_tmp = sapply(1:dim(Graphs)[3], function(i) c(Graphs[,,i]))
    #  duplicated_graphs = duplicated(t(Graphs_tmp))
      
    #}
    G_star = out_baseline[,,sample(S_base - burn_base, 1)]
    Graphs = abind(Graphs, G_star)
    K_star = dim(Graphs)[3]
    
    logProbs = matrix(nrow = n, ncol = K_star) 
    
    for(k in 1:(K_star-1)){
      
      cliques_separators_k = mpd(Graphs[,,k])
      
      cl_k = lapply(1:length(cliques_separators_k$cliques), function(i) as.integer(cliques_separators_k$cliques[[i]]))
      se_k = lapply(1:length(cliques_separators_k$separators), function(i) as.integer(cliques_separators_k$separators[[i]]))
      
      #logProbs[,k] = sapply(1:n, function(i) prob_ik_nonempty(Xk = X[xi == k,, drop = FALSE], cl.list = cl_k, se.list = se_k, 
      #                                                        member = (xi[i] == k), x = X[i,], a = a, I.cal = I.cal))
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
    
    
    ###################
    ## Update graphs ##
    ###################
    
    set = 1:K
    
    for(k in set){
      
      ## Propose new graph
      
      G = Graphs[,,k]
      G_move = move_dec(A = G)
      
      G_star        = G_move$A_new  # adjacency matrix of the proposed graph
      nodes_star    = G_move$nodes  # nodes (u,v) involved in the local move
      type.operator = G_move$O_new  # the type of the applied operator in the local move (1: Insert; 2: Delete)
      
      Xk =(X[xi == k, , drop = FALSE])
      Nk = plyr::count(Xk)
      
      ## Compute logprior ratio --> in cpp
      
      #logprior.star = lgamma(sum(G_star)/2 + a_pi) + 
      #  lgamma(q*(q-1)/2 - sum(G_star)/2 + b_pi)
      
      #logprior = lgamma(sum(G)/2 + a_pi) + 
      #  lgamma(q*(q-1)/2 - sum(G)/2 + b_pi)
      
      #logprior.ratio = logprior.star - logprior
      logprior.ratio = logprior_ratio(as.vector(G_star), as.vector(G), a_pi, b_pi, q)
      
      ## Compute the marginal likelihood ratio
      
      if(type.operator == 1){
        
        # Type 1 : edge insertion
        # Find C.star : unique clicque in G_star containing both u and v
        
        cliques.star = mpd(G_star)$cliques
        
        C.star = cliques.star[sapply(X = lapply(X = cliques.star, FUN = intersect, nodes_star), FUN = length) == 2]
        C.star = as.numeric(unlist(C.star))
        
        C.u = setdiff(C.star, nodes_star[1])
        C.v = setdiff(C.star, nodes_star[2])
        C.0 = setdiff(C.star, nodes_star)
        
        m.ratio = marg_S_test(S = C.star, N = Nk, a = a, I_cal = I.cal) + 
          marg_S_test(S = C.0, N = Nk, a = a, I_cal = I.cal) -
          marg_S_test(S = C.u, N = Nk, a = a, I_cal = I.cal) -
          marg_S_test(S = C.v, N = Nk, a = a, I_cal = I.cal)
        
      } else{
        
        # Type 2 : edge deletion
        # Find C.star : unique clicque in G containing both u and v
        
        cliques.star = mpd(G)$cliques
        
        C.star = cliques.star[sapply(X = lapply(X = cliques.star, FUN = intersect, nodes_star), FUN = length) == 2]
        C.star = as.numeric(unlist(C.star))
        
        C.u = setdiff(C.star, nodes_star[1])
        C.v = setdiff(C.star, nodes_star[2])
        C.0 = setdiff(C.star, nodes_star)
        
        #m.ratio = marg_S(S = C.u, N = Nk, a = a, I.cal = I.cal) + 
        #  marg_S(S = C.v, N = Nk, a = a, I.cal = I.cal) -
        #  marg_S(S = C.star, N = Nk, a = a, I.cal = I.cal) -
        #  marg_S(S = C.0, N = Nk, a = a, I.cal = I.cal)
        
        m.ratio = marg_S_test(S = C.u, N = Nk, a = a, I_cal = I.cal) + 
          marg_S_test(S = C.v, N = Nk, a = a, I_cal = I.cal) -
          marg_S_test(S = C.star, N = Nk, a = a, I_cal = I.cal) -
          marg_S_test(S = C.0, N = Nk, a = a, I_cal = I.cal)
        
        
      }
      
      ## Metropolis Hastings ratio and acceptance/rejection
      
      ratio_k = min(0, m.ratio + logprior.ratio)
      
      if(log(runif(1)) < ratio_k){
        G = G_star
      }
      
      Graphs[,,k] = G
      
    }
    
    Graphs_chain[[t]] = Graphs
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
  
  ## Subject-specific Posterior Probabilities of edge Inclusion (PPIs)
  
  graph_probs = array(NA, c(q, q, n))

  for(i in 1:n){
 
    probs_i = sapply((burn + 1): S, function(t) c(Graphs_chain[[t]][,,Xi_chain[i,t]]))
    graph_probs[,,i] = matrix(rowMeans(probs_i), q, q)

  }
  
  return(list(Graphs_chain = Graphs_chain, Xi_chain = Xi_chain, alpha0_chain = alpha0_chain,
              simil_probs = simil_probs, graph_probs = NULL))
  
  
}



