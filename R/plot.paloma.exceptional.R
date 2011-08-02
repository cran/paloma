plot.paloma.exceptional <- function( x, ... ) {

  if( is.null( x ) ) {
    stop("NULL object")
  }
  if( ! is.paloma.exceptional( x ) ) {
    stop("not an exceptional.motif object")
  }

  
  # Frames to display
  nb.frame <- 4
  nb.col <- 2
  nb.lin <- 2
  if ( nb.frame == 1) {
    nb.col <- 1
    nb.lin <- 1
  } else if ( nb.frame == 2) {
    nb.col <- 2
    nb.lin <- 1
  }
  
  par(mar=c(1, 0, 2, 0))
  par(mfrow=c(nb.lin, nb.col))

  k.set <- 1:length( x )
  nb.plot <- 0
  
  for ( k in k.set ) {
    k.value  <- x[[k]]$k
    directed <- x[[k]]$directed
    Labels <- paste( 0:( k.value-1) )

    x.k <- x[[k]]$motifs
    m.set <-  1:length(x.k$m.mp )
    for( i in m.set ) {
      m <-  x.k$adjacency[[i]]
      m.mp.set <- 1:length(x.k$m.mp[[i]] )
      for( j in m.mp.set ) {
        del.class <-  x.k$m.mp[[i]][[j]]$Del
        pvalue <- x.k$m.mp[[i]][[j]]$PValue
        NodeToClass <- rep( 1, k.value)
        NodeToClass[ del.class + 1] <- 2
        title <- paste( "k = ",k.value, ", p-value = ",
                       sprintf( "%.2e", pvalue), sep="" )
        Gplot( t( matrix( m, k.value, k.value) ) , class=NodeToClass,
              label= Labels, main=title, directed=directed, ...)
        # Next plot
        nb.plot <- nb.plot + 1

        if ((nb.plot %% 4) == 0 )
          readline(prompt="\nPress return for next plot...")

      }
    }
  }
  
  par( mfrow=c(1, 1) )
}
