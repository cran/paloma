#
#    Internal procedures
#


TestPseudoMixNet <- function( nexp, N, Alpha=NULL, Pi=NULL, motif.size,
                             directed=FALSE ) {
  # Graphics
  # par(mar=c(2, 4, 0, 1))
  # par(mfrow=c(4, 2))
  
  # Generating Pi
  nclass <- length( Alpha )
  
  exp <- rep(0,nexp)
    
  means <- rep(0 ,nexp+1 )
  
  sum <- 0
  motifs <- list()
  motifs[[1]] <- list()
  motifs[[2]] <- list()
  # motif[[1]] : adjacency
  # motif[[2]] : count
  n.motif <- 0
  for( i in 1:nexp) {
    res <- countMotifInPseudoMixNet( N, Alpha, Pi, motif.size,
                                  directed=directed )
    if( ! is.null( res ) ) {
      for( k in 1:length( res$adjacency ) ) {
        found <- FALSE
        if ( length( motifs[[1]] ) > 0 ){
          for( j in 1:length( motifs[[1]] ) ) {
            if ( all( motifs[[1]][[j]] == res$adjacency[[ k ]] ) ) {
              motifs[[2]] [[ j ]] [i] <- res$count[[ k ]]
              found <- TRUE
            }
          }
        }
        if ( ! found ) {   
          n.motif <- n.motif + 1
          motifs[[1]] [[ n.motif ]] <- res$adjacency[[ k ]]
          motifs[[2]] [[ n.motif ]] <- rep(0, nexp)
          motifs[[2]] [[ n.motif ]][i] <- res$count[[ k ]]
        }
      }
    }
  }

  
  for( i in 1:length( motifs[[1]] ) ) {
    mat     <- matrix( motifs[[1]][[i]], motif.size, motif.size)
    edges   <- AdjMat2Edges( mat, directed=directed)
    ## mean.th <- computePseudoMixNetMean( N, Alpha, Pi,
    ##                                    edges, directed=directed )
    ##
    mean.th <- computeModelMean( Alpha, Pi,
                                edges, directed=directed )
    
    for( li in 1:motif.size) {
    str <- "  "
      for( lj in 1:motif.size) {
        str <- paste( str, mat[li, lj], sep=" " )
      }
      cat("\n", str)
    }

    sum.motif <- sum( motifs[[2]][[ i ]] )
    mean <- sum.motif / nexp
    std.dev <- sum( (motifs[[2]][[ i ]] - mean )^2 )
    std.dev <- sqrt( std.dev / nexp )
    fmt <- "             %d %f  %f %f \n"
    str <- sprintf( fmt,  sum.motif,
           mean, std.dev, mean.th ) 
    cat(str)
    
  }
    
}


# Generate Minet graph miture
# alpha : number of nodes per class. The number of classes  is len(alpha)
# Pi : Pi[i,j] (symmetric/unsymmetric) matrix probabilities to connect nodes
#     from class i to nodes belonging to class j 
# symmetric : specifies if the returned edges list is undirected
# diag : if TRUE removes loops fron edge list
#
# Return an edge list[2, n_edges]
#
PseudoMixNet <- function( N, Alpha, Pi, storage="matrix",
                          directed = FALSE, diag=FALSE){

  # Case where Pi is a scalar
  Pi <- as.matrix( Pi )

  Alpha.n <- Alpha
  
  # n nodes is used to take into account the rounding effects
  n <- sum( Alpha.n)
    
  start.index <- vector()
  start.index[1] <- 1
  
  for( i in 1:( length(Alpha.n) )) {
    start.index[i+1] <- start.index[i] +  Alpha.n[i] 
  }
  adj <- matrix( 0, n, n)   
  for( i in 1:length(Alpha.n) ) {
    for( j in 1:length(Alpha.n) ) {
      i.beg = start.index[i] 
      i.end = start.index[i+1]-1
      j.beg = start.index[j] 
      j.end = start.index[j+1]-1
      if( (Alpha.n[i] != 0) & (Alpha.n[j] != 0) )
      adj[ i.beg:i.end, j.beg:j.end] <- matrix(
                             rbinom( Alpha.n[i] * Alpha.n[j], 1, Pi[i,j]),
                                                  Alpha.n[i],  Alpha.n[j]) 
    }
  }

  # 
  # Building Edge matrix
  #

  # Diagonnal terms
  if(!diag){
    for (i in 1:n){
    adj[i , i] = 0   
    }
  }
  if (storage == "edges") {
    if ( !directed  ){
      m <- t( which( (adj==1) & (upper.tri(adj)), arr.ind=TRUE) )
    } else {
      m <- t( which( ( t(adj) == 1 ) , arr.ind=TRUE) )
    }
  } else if( storage == "matrix" ) {
    if ( !directed ) {
      m <- array( 0, dim( adj ) ) 
      adj[ lower.tri(adj) ] <- 0
      m <- adj + t( adj )
    } else {
      m <- adj
    }
  }
  return(m)
}


#
#                    TEST
#
countMotifInPseudoMixNet <- function( N, Alpha=NULL, Pi=NULL, k,
                                      directed=FALSE, display=FALSE ) {
# Generates graph according to a model (MixNet, EDD, ER)
# and counts the motif occurences.
# Used to test graphic models
# 
# Alpha, PI : used by MixNet and ER
# Alpha : contains degree distribution for EDD
  
  G <- PseudoMixNet(N, Alpha, Pi, directed = directed, diag=FALSE,
                      storage="edges")

  res <- countAllMotifs( G, k, directed )

  cat( " MOTIF LIST \n" )
  for( i in 1:length( res$count ) ) {
    cat( "m : ", res$adjacency[[i]] , " ", res$count[i], "\n")
  }

  return ( res )
}


computePValues <- function( G, k, NodeToClass, Pi, pvalue=1.0, directed=FALSE ) {
  # Transpose G, renumber nodes and add an end index at the end of the list
  G <- testGraph(G )

  n <- max(G) + 1 
  n.edges <- dim( G )[2]

  #  Compute pi erdos (used if NoteToClass is not specified)
  pi.erdos <- n.edges / (n * (n -1))

  if( ! directed )
      pi.erdos <- pi.erdos * 2

  res <- testPseudoMixNetParameters( NodeToClass, Pi, n, pi.erdos )
  NodeToClass <- res[[1]]
  Pi          <- res[[2]]
  
  nb.classes <- dim( Pi ) [1]
  
  directed.int <- directed*1
  res <- .Call( "get_M_Mp_PValues_C", G, n, k, NodeToClass, Pi, nb.classes,
               as.double(pvalue), directed.int, PACKAGE="paloma")
}
