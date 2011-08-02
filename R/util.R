#
#  Test the validity of the edge list
# 
testGraph <- function( G ) {
# G : edge list (graph)  
  max.G <- max(G)
  min.G <- min(G)

  n <-  max.G + 1 - min.G

  if ( dim( G ) [2] != 2 ) {
    stop( "Not an edge list edges[2,n]")
  }

  G <- t(G) - min.G
  G <- cbind( G, c(-1,-1) )

  return (G)
  
}

#
#  Test the validity of the pvalue
# 
testPValue <- function( pvalue ) {
  if( pvalue < 0 ) {
    stop("bad p-value")
  }
}

#
#  Test the validity of the deletion class
# 
testDelClass <- function ( del.class ) {

  arg.name = "del.class"
  if(  is.null( del.class )  ) {
    stop( "bad ", arg.name, " argument")
  } else if( is.matrix( del.class ) ) {
    del.class <- as.vector( del.class ) 
  }

  return ( del.class[1] )
}

#
#  Test the validity of a motif
#  In : Adj. matrix (vector or matrix)
#  Out : edge list
# 
testMotif <- function( m, directed ) {
  
  if( is.null( m ) ) {
    stop("bad motif argument")
  } else if ( is.vector( m ) ) {
     k = round( sqrt( length( m ) ) )
     m <- matrix( m, k, k )
     edges <- AdjMat2Edges( m, directed=directed )
  } else if (is.matrix( m ) ) {
    if( dim(m)[1] != dim(m)[2] ) {
      stop("motif not a squared matrix")
    } else {
      edges <- AdjMat2Edges( m, directed=directed )
    }
  } else {
    stop("bad motif argumet")
  }
  edges <- edges - 1
  return( edges ) 
}

#
#  Type-cast matrix into a vector
#
testVector <- function ( vect, k, arg.name="vector" ) {

  if( is.null( vect ) ) {
    stop( "bad ",arg.name, " argument \n")
  } else if( is.matrix( vect ) ) {
    vect <- as.vector( vect ) 
  }
  if( length( vect) != k*k) {
    stop( "bad ",arg.name, " array size \n")
  }
  
  return( vect )
}

#
#  Tests the validity of a set
#  Returns its intersection with 1: default.length
# 
testSet <- function ( vect, default.length ) {

  if( is.null( vect ) ) {
    if ( default.length == 0 ) {
      return( NULL )
    } else {
      return(  1:default.length )
    }
  } else {
      return( intersect( vect, 1:default.length ) )
  } 
}

#
#  Test the validity of NodeToClass
#  and Pi
#  If NodeToClass, Pi are NULL there
#  are filled with an ER model
testPseudoMixNetParameters <- function ( NodeToClass, Pi, nnodes, pi.erdos ) {
  
  if(  is.null( NodeToClass ) & is.null( Pi ) ) {
    NodeToClass <- rep( 0, nnodes )
    if( pi.erdos == 0 ) {
      stop("NodeToClass is required")
    }
    Pi <- matrix ( c( pi.erdos ), 1, 1)
  } else if(  is.null( NodeToClass ) | is.null( Pi ) ) {
    stop( "for 'MixNet' model, NodeToClass and Pi are required ")
  } else if ( ! is.vector( NodeToClass ) | ! is.matrix( Pi ) ){
    stop( "NodeToClass must be a vector and Pi must be a matrix ")
  } else if( length( NodeToClass ) != nnodes ) {
    stop( "Error: NodeToClass and Pi have not the same dimensions \n")
  } else {
    min.alpha <- min( NodeToClass )
    max.alpha <- max( NodeToClass )
    if(  ( max.alpha -  min.alpha + 1 ) != dim( Pi ) [1] ) {
      cat( "Error: Bad number of classes in the NodeToClass vector \n")      
      cat( "     : Indexes not contiguous\n")
      stop()
    }
    NodeToClass <- NodeToClass - min.alpha
  }
  return( list( NodeToClass, Pi) )
}

#
#  Convert an Adjacency matrix to an edge list
#
AdjMat2Edges<-function(x, directed=FALSE, verbose=FALSE ) {

  # 
  if (dim(x)[1] == dim(x)[2]){
    
    # Adjacency matrix case
    nbrNodes<-dim(x)[1]
    NotConnected <- which( (rowSums( x ) + colSums(x)) == 0)
    if ( length(NotConnected) > 0 ) {
        if (verbose)
         cat("Some nodes are not connected to the network \n")
    } 
    if ( directed ) {
      m <- t( which( (t(x)==1) , arr.ind=TRUE) )
    } else {
      m <- t( which( (x==1) & (upper.tri(x)), arr.ind=TRUE) )
    }
  } else {
    stop( "Error: Bad matrix dimensions " )
  }

  return( m )
}

#
# Names the list fields
# and set a class attribute
#
setExceptRecordNames <- function( res, directed ) {

  # Remove NULL
  for( k in length( res ):1 ) {
    if ( is.null( res[[k]] ) ) {
      res[[k]] <- NULL
    }
  }

  if ( length( res ) != 0 ) {
    for( k in length( res ):1 ) {
      adj  <- res[[k]][[2]][[1]]
      m.mp <- res[[k]][[2]][[2]]
      for( i in (length( adj ) ):1 ) {
        # Remove adjacency
        if ( is.null( adj[[i]] ) ) {
          adj[[i]] <- NULL
        }         
        # Remove m list
        if ( is.null( m.mp[[i]] ) ) {
          m.mp[[i]] <- NULL
        } else {
          for( j in (length( m.mp[[i]] ) ):1 ) {
            if ( is.null( m.mp[[i]][[j]] ) ) {
              m.mp[[i]][[j]] <- NULL
            }
          }
        }
      }
      res[[k]][[2]][[1]] <- adj
      res[[k]][[2]][[2]] <- m.mp
    }

  
    for( k in 1:length( res ) ) {
      res[[k]]<-list ( res[[k]][[1]], directed, res[[k]][[2]] )
      names( res[[k]] ) <- c( "k", "directed", "motifs")
      names( res[[k]]$motifs ) <- c( "adjacency", "m.mp")
    
      # Loop over motif m 
      for( i in 1:length( res[[k]]$motifs$m.mp ) ) {
        for( j in 1:length( res[[k]]$motifs$m.mp[[i]] ) ) {
        
          names( res[[k]]$motifs$m.mp[[i]][[j]] ) <-
            c( "Del", "Filter", "PValue" )  
        }
      }
    }
    class(res) <- "paloma.exceptional"
  } else {
    
    res <- NULL
  }
  
  return( res )
}

#
# Names the list fields
# and set a class attribute
#
setOccRecordNames <- function( res, directed ) {
  
  for( k in 1:length( res ) ) {
    res[[k]]<-list (res[[k]][[1]], directed, res[[k]][[2]] )
    names( res[[k]] ) <- c( "k", "directed", "motifs")
    names( res[[k]]$motifs ) <- c( "adjacency", "m.mp")
    
    # Loop over motif m 
    for( i in 1:length( res[[k]]$motifs$m.mp ) ) {
      for( j in 1:length( res[[k]]$motifs$m.mp[[i]] ) ) {
      
        names( res[[k]]$motifs$m.mp[[i]][[j]] ) <-
                  c( "DelClass", "Occur", "Ext" )
        
      }
    }
  }
  
  class(res) <- "paloma.occurrences"
  
  return( res )
}

# Test class
is.paloma.occurrences <- function( x ) {
  
  if (class(x)=="paloma.occurrences") TRUE else FALSE
}

# Test class
is.paloma.exceptional <- function( x ) {
  
  if (class(x)=="paloma.exceptional") TRUE else FALSE
}

# Display object
print.paloma.occurrences <- function( x, k=NULL, m=NULL, mp=NULL, ... ) {
# x : paloma.occurrences object
# k, m, mp : k, motif, deletion classes indexes
  
  k.set <- testSet ( k, length( x ) )
  if( ! is.paloma.occurrences( x ) ) {
    stop("not an occurrence object")
  }

  for ( k in k.set ) {
    cat( " \n" )
    cat( " \n" )
    cat ( "motif size, k =",  x[[k]]$k , ", directed =", x[[k]]$directed, "\n" )
    x.k <- x[[k]]$motifs
    # Loop over motif m 
    m.set <-  testSet ( m, length(x.k$m.mp ) )
    for( i in m.set ) {
      cat( " \n" )
      cat( "  m : ", x.k$adjacency[[i]] , "\n")
      m.mp.set <-  testSet ( mp, length(x.k$m.mp[[i]] ) )
      for( j in m.mp.set ) {
        cat( "    m' : \n" )
        cat( "      Del class :", x.k$m.mp[[i]][[j]]$DelClass ,"\n" )
        cat( "      Occurences: \n" )
        for ( l in 1:length( x.k$m.mp[[i]][[j]]$Occur ) ) {
          cat( "        ", x.k$m.mp[[i]][[j]]$Occur[[l]] ," [" )
          
          for ( m in 1:length( x.k$m.mp[[i]][[j]]$Ext[[l]] ) ) {
            if (m != 1)
              cat(", ")
            cat( x.k$m.mp[[i]][[j]]$Ext[[l]][[m]] )
          }
          cat("] \n")
        }
      }
    }
  }
}

# Display object
print.paloma.exceptional <- function( x, k=NULL, m=NULL, mp=NULL, ... ) {
# x : paloma.exceptionnal object
# k, m, mp : k, motif, deletion class indexes

  k.set <- testSet ( k, length( x ) )
  if( ! is.paloma.exceptional( x ) ) {
    stop("not an paloma.exceptional object")
  }

  for ( k in k.set ) {
    cat( " \n" )
    cat( " \n" )
    cat ( "motif size, k =",  x[[k]]$k , ", directed =", x[[k]]$directed, "\n" )
    x.k <- x[[k]]$motifs
    # Loop over motif m 
    m.set <-  testSet ( m, length(x.k$m.mp ) )
    for( i in m.set ) {
      cat( " \n" )
      cat( "  m : ", x.k$adjacency[[i]] , "\n")
      m.mp.set <-  testSet ( mp, length(x.k$m.mp[[i]] ) )
      for( j in m.mp.set ) {
        cat( "    m' : \n" )
        cat( "      Del class :", x.k$m.mp[[i]][[j]]$Del ,"\n" )
        cat( "      H. Filter :", x.k$m.mp[[i]][[j]]$Filter[1] ,"\n" )
        cat( "      p-value   :", x.k$m.mp[[i]][[j]]$PValue ,"\n" )
      }
    }
  }
}

#
#  Remove over-represented solution (motif, del. class)
#  if a smaller motif m' is overrepresented wit a
#  compatible deletion class
#
filterPValue <- function( ext.k ) {

  if ( ! is.null( ext.k ) ) {
    for ( k in length( ext.k ):1 ) {
      res <- ext.k[[k]]$motifs
      
      remove.list <- vector()
      for( i in length( res$m.mp ):1 ) {
        for( j in length( res$m.mp[[i]] ):1 ) {
          if (  res$m.mp[[i]][[j]]$Filter[1] == 0) {
            res$m.mp[[i]][[j]] <- NULL
          }
        }
        if ( length( res$m.mp[[i]] ) == 0 ) {
          res$m.mp[[i]] <- NULL
          res$adjacency[[i]] <- NULL
        }
      }
      
      if ( length(res$m.mp) == 0 ) {
        res$m.mp      <- NULL
        res$adjacency <- NULL
        ext.k[[k]]$k      <- NULL
        ext.k[[k]]$motifs <- NULL
        ext.k[[k]] <- NULL
        
      } else {
        ext.k[[k]]$motifs <- res
      }
    }

    if( length (ext.k ) == 0 )
      ext.k <- NULL
  }
  return( ext.k )
}  

packageLoaded <- function(name) {
  return ( 0 != length( grep( paste("^package:", name, "$", sep=""), search())) )
}
