
getCanonic <- function( m, directed=FALSE ) {
# motif : vector or matrix
  
  k = round( sqrt( length( m ) ) )
  
  m <- testVector( m, k, "adjacency matrix")

  directed.int <- directed*1

  return( .Call( "getCanonic_C", m, k, directed.int, PACKAGE="paloma") )
}

getAllExtensions <- function( G, k, directed=FALSE  ) {
# G : edge list
# motif : vector or matrix
  
  # Transpose G, renumber nodes and add end indexes at the end of the list
  G <- testGraph(G )

  # Number of nodes
  n <- max(G) + 1
  
  directed.int <- directed*1

  res <- .Call( "get_M_Mp_Occurrences_C", G, n, k, directed.int,
               PACKAGE="paloma")
  
  # Names the list fields 
  res <- list( list( k=c(k), motifs=c(res) ) )
  res <- setOccRecordNames( res, directed )

  return( res )
}

getExtensions <- function( G, motif=NULL, del.class=NULL, directed=FALSE ) {
# G : edge list
# motif : vector or matrix
# del.class : node index
  
  # Transpose G, renumber nodes and add an end index at the end of the list
  G <- testGraph(G )

  # Number of graph nodes and edges
  n <- max(G) + 1
  n.edges <- dim( G )[2]

  # Get motif node number
  k = round( sqrt( length( motif ) ) )
  
  # Test k value
  if ( k < 2) {
    cat("Error: bad values of k \n")
    stop
  }
  
  motif <- testVector( motif, k, arg.name="canonic" )

  del.class <- testDelClass ( del.class )
  
  directed.int <- directed*1
  
  res <- .Call( "get_particular_M_Mp_Occurrences_C", G, n, k,
               motif, del.class,
               directed.int, PACKAGE="paloma")

  # Names the list fields 
  res <- setOccRecordNames( res, directed )
  
  return( res )
}

getExceptional <- function( G,
                            kmin=2, kmax=2,
                            NodeToClass=NULL, Pi=NULL,
                            pvalue=0.001,
                            directed=FALSE,
                            h.filter=TRUE) {

  # Transpose G, renumber nodes and add an end index at the end of the list
  G <- testGraph(G )

  # Number of graph nodes and edges
  n <- max(G) + 1
  n.edges <- dim( G )[2]

  #  Compute pi erdos (used if NoteToClass is not specified)
  pi.erdos <- n.edges / (n * (n -1))

  if( ! directed )
      pi.erdos <- pi.erdos * 2

  res <- testPseudoMixNetParameters( NodeToClass, Pi, n, pi.erdos )
  NodeToClass <- res[[1]]
  Pi          <- res[[2]]
  
  testPValue( pvalue )

  # kmin, kmax
  if ( kmin > kmax ) {
    kmax = kmin
  }
  if ( kmin < 2) {
    stop("bad values of kmin, kmax \n")
  }

  nb.classes <- dim( Pi ) [1]
  
  directed.int <- directed*1
  
  res <- .Call( "get_Exceptional_M_Mp_C", G, n, kmin, kmax,
                NodeToClass, Pi, nb.classes,
                pvalue,
                directed.int, PACKAGE="paloma")
  
  # Names the list fields
  res <- setExceptRecordNames( res, directed )

  # Filters motifs induced by smaller motifs 
  if( h.filter ) {
    res <- filterPValue( res )
  }
  return( res )
}

countAllMotifs <- function( G, k, directed=FALSE ) {
# Count induced (exact) motifs
# G : edge list of the network
# k : motif node number

  # Transpose G, renumber nodes and add an end index at the end of the list
  G <- testGraph(G )
  n <- max( G) + 1
  
  if ( k < 2) {
    stop("bad values of the motisize (k)")
  }
  
  G.edge.nbr <- dim( G ) [2]
  
  directed.int <- directed*1
  
  res <- .Call( "countAllExactMotifs_C", G, n, k, directed.int,
                PACKAGE="paloma")

  names( res ) <- c( "adjacency", "count")

  # No motifs case
  if( res$count[1] == 0 ) res <- NULL
  return( res )
}

computeModelMean <- function( NodeToClass, Pi, motif, directed=FALSE ) {
# Compute expectation E of a motif under a MixNet model
# NodeToClass, Pi : MixNet model parameters
# motif : adjacency matrix
  
  # Number of nodes
  N <- length( NodeToClass )

  pi.erdos <- 0
  res <- testPseudoMixNetParameters( NodeToClass, Pi, N, pi.erdos )
  
  NodeToClass <- res[[1]]
  Pi          <- res[[2]]

  # Number of classes
  Q <- dim( Pi ) [1]

  # Alpha : Number of nodes per Class
  Alpha <- rep(0, Q)
  for (i in 1:Q) {
   Alpha [i] <- length( which( NodeToClass == (i-1) ) )
  }
  
  directed.int <-  directed*1

  motif <- testMotif( motif, directed )
  
  nb.edges <- dim(motif)[2]
  
  res <- .Call( "computePseudoMixNetMean_C", N, Q,
               Alpha, Pi,
               motif, nb.edges, directed.int, PACKAGE="paloma")
  return(res)
}




