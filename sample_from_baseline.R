library(gRbase)
library(pcalg)
source("MCMC/move_decomposable.R")

sample_baseline_dec = function(S, burn, q, a_pi, b_pi){
  
  ## Function to sample from the baseline over the space of decomposable graphs
  
  ## INPUT
  
  ## S    : number of draws
  ## burn : burn in period
  ## q    : number of nodes in the graphs
  
  ## b_pi, b_pi : hyper-parameters of the Beta prior on probability of edge inclusion pi
  
  ## OUTPUT
  
  ## (S - burn) DAGs sampled from the space of DAGs with q nodes
  
  G_chain = array(NA, c(q, q, S))
  
  # set intial value
  
  G = matrix(0, q, q); G_chain[,,1] = G
  
  for(s in 2:S){
    
    move_star = move_dec(A = G)
    G_star    = move_star$A_new
    
    # multiplicity correction (log)prior
    
    logprior.new = lgamma(n.edge(G_star) + a_pi) +
      lgamma(q*(q-1)/2 - n.edge(G_star) + b_pi)
    
    logprior.old = lgamma(n.edge(G) + a_pi) +
      lgamma(q*(q-1)/2 - n.edge(G) + b_pi)
    
    logprior = logprior.new - logprior.old
    
    # acceptance ratio
    
    ratio = min(0, logprior)
    
    # accept graph
    
    if(log(runif(1)) < ratio){
      
      G = G_star
      
    }
    
    G_chain[,,s] = G
    
  }
  
  return(list(G_chain = G_chain[,,(burn + 1):S]))
  
}
