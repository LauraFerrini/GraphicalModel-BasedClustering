library(gRbase)

n.edge = function(A){
  sum(A)/2
}

check_decomp = function(A){
  
  check = gRbase::mcs(as(A,'graphNEL'))
  return(out = length(check) > 0)
  
}

# Two possible actions

actions = c("iu", "du")

# Insert an undirected edge

iu = function(A, nodes){
  A[nodes[1],nodes[2]] = A[nodes[2],nodes[1]] = 1
  return(A)
}

# Delete an undirected edge

du = function(A, nodes){
  A[nodes[1],nodes[2]] = A[nodes[2],nodes[1]] = 0
  return(A)
}

move_dec = function(A){
  
  ## INPUT
  
  ## A: adjacency matrix of decomposable graph G
  
  ## OUTPUT
  
  ## A_new: adjacency matrix of a new adjacent decomposable graph
  ## O_new: the type of applied operator (1 for insert, 2 for delete an edge in A)
  ## nodes: the two nodes involved in the move
  
  A_und = A
  A_und[upper.tri(A, diag = TRUE)] = NA
  
  # two types of moves: (1) insert an undirected edge, (2) delete an undirected edge
  
  iu_set = c()
  du_set = c()
    
  # set of nodes for iu
  
  nodes = which(A_und == 0,TRUE)
  if(length(nodes) != 0){
    iu_set = cbind(1, nodes)
  }
  
  # set of nodes for du
  
  nodes = which(A_und == 1,TRUE)
  if(length(nodes) != 0){
    du_set = cbind(2, nodes)
  }
  
  # set of possible moves
  
  O = rbind(iu_set, du_set)
  
  repeat {
    
    a   = sample(dim(O)[1],1)
    act = actions[O[a,1]] # type of the sampled operator (1 or 2)
    
    act_to_eval   = paste0(actions[O[a,1]],"(A=A,c(",as.vector(O[a,2]),",",as.vector(O[a,3]),"))")
    graph_to_eval = eval(parse(text = act_to_eval))
    
    if(check_decomp(graph_to_eval) == TRUE){
      break
    }
  }
  
  action = paste0(actions[O[a,1]],"(A=A,c(",as.vector(O[a,2]),",",as.vector(O[a,3]),"))")
  A_new = eval(parse(text = action))

  return(list(A_new = A_new, O_new = O[a,1], nodes = c(O[a,2], O[a,3])))
  
}