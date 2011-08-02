\name{getExtensions}
\alias{getExtensions}

\title{
 Search occurrences and extensions in a network for a given motif 
}

\description{
  Provides all the \emph{positions} (k-1 vertices) and all their
  possible one node extensions for a given motif (induced motif of k
  nodes) in the input graph.
  
  The output list is sorted according to the deletion class
  (removing one node form the motif). Then for each
  deletion class the
  motif occurrences are stored according to their \emph{positions} (k-1
  vertices, one
  node belonging to the given deletion class is removed) with all the
  possible one node extensions for a specific \emph{position} (k-1 vertices).
}

\usage{
getExtensions(G, motif, del.class, directed = FALSE)
}

\arguments{
  \item{G}{
    the graph is described as an edge list in a matrix
    \code{edges[ nbr_of_edges, 2]}. The first column contains the source
    nodes, whereas the second one concerns the destination nodes. The node
    indexes must be contiguous.
  }
  \item{motif}{
    adjacency matrix, it could be a vector or matrix type.
    The canonical form is required.
  }
  \item{del.class}{
    node index in the canonical form \code{motif} which selects a
    deletion class. If \code{-1} is set then all deletion classes
    are considered.
  }

  \item{directed}{
    specifies if the graph is directed (or not).
  }
}



\value{
  The result is returned as an \code{paloma.occurrences} class.
  \describe{
    \item{}{
      \code{rv[[1]]} list concerning a motif size (the \emph{kth} ).
    }
    \item{}{
      \code{rv[[1]]$k} motif node number. 
    }
    \item{}{
      \code{rv[[1]]$directed} the motis are directed or not (logical). 
    }
   \item{}{
      \code{rv[[1]]$motifs$adjacency[[1]]} adjacency matrix of
      the chosen motif  (canonical adjacency matrix).
    }    
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[1]]} list describing couples
      (chosen motif, deletion class).
    }
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[1]][[j]] } describes all
      \emph{positions} (k-1 vertices) in the graph with all the possible extensions (at one node) for
      the chosen motif and \emph{jth} deletion class.
    }
    \item{}{   
      \code{rv[[1]]$motifs$m.mp[[1]][[j]]$DelClass } nodes having the
      same deletion class (node indexes of the canonical motif) for the
      chosen motif and \emph{jth} deletion class.
    }
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[1]][[j]]$Occ[[l]]} \emph{lth} 
      \emph{position} (k-1 vertices) for the chosen motif with the \emph{jth}
      deletion class.
    }
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[1]][[j]]$Ext[[l]]} all possible
      extensions of the \emph{lth} \emph{position}
  (\ldots{}\code{$m.mp[[1]][[j]]$Occ[[l]]}, containing k-1 vertices). 
    }
 }
 Summary :
    \tabular{ll}{
      rv[[1]] \tab motifs having the same motif size \cr
      rv[[1]]$k \tab the motif size value \cr
      rv[[1]]$k \tab motifs are directed (or not) \cr
      rv[[1]]$motifs \tab motifs list for the given size \cr
      rv[[1]]$motifs$adjacency[[1]] \tab chosen motif adjacency matrix \cr
      rv[[1]]$motifs$m.mp[[1]] \tab occurrences and extensions for the
      chosen motif \cr
      rv[[1]]$motifs$m.mp[[1]][[j]] \tab occurrences and extensions for
  the chosen motif \cr
      \tab and a deletion class (\code{jth}) \cr
      
      rv[[1]]$motifs$m.mp[[1]][[j]]$DelClass \tab elements of the deletion class  \cr
      rv[[1]]$motifs$m.mp[[1]][[j]]$Occ[[l]] \tab \code{lth} \emph{position} (k-1 nodes) \cr  
      rv[[1]]$motifs$m.mp[[1]][[j]]$Ext[[l]] \tab possible extensions for the
      \emph{lth} \emph{position} (k-1 nodes)\cr
    }
  
    Remark: for performance reasons, the node indexes have been
    renumbered in the
    \emph{position} (k-1 nodes) field (\code{Occ}) and the extension field
    (\code{Ext}). 
    To get initial numbering add the original minimum index.
  }

\author{G. Grasseau, E. Birmele}

\section{Warning }{
Setting a great value to the motifs size could generate
huge occurrence lists and thus cause memory allocation problems.
}


\seealso{
\code{\link{paloma-package}}
}

\examples{

## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The connectivity matrix is stored in "example$pi"
# and the mapping node/class is stored in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# Print exceptionnal motifs
out


# Get all positions (with k-1 nodes) with the deletion class of node 0
# for the 11th exceptional motif of size k=5
u <- getExtensions( example$graph, out[[3]]$motifs$adjacency[[11]], 0 )

# Print occurrences
u

}

\keyword{ graphs }
