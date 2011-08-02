\name{getAllExtensions}
\alias{getAllExtensions}

\title{
 Searchs all occurrences and extensions in a network. 
}

\description{
  Search all induced motifs for a given node number \code{k} in a graph.
  
  The same motif occurrences  (having the same canonical from) are put
  together. For this occurrence list several sets are built according
  to the deletion class (removing one node form the motif). Then for
  each deletion class the
  motif occurrences are stored according to their \emph{position} (k-1
  vertices, one
  node belonging to the given deletion class is removed) with all the
  possible one node extensions for the considered \emph{position} (k-1
  vertices).
  
  The motif discovery algorithm implemented is the one used in the
  \code{FanMod} code and the canonical form is computed with the Nauty
  library (see \code{\link{paloma-package}})
}

\usage{
getAllExtensions(G, k, directed = FALSE)
}

\arguments{
  \item{G}{
    the graph is described as an edge list in a matrix
    \code{ edges[nbr_of_edges, 2] }. The first column contains the source
    nodes, whereas the second one concerns the destination nodes. The node
    indexes must be contiguous and the graph must not contain
    self-loops nor redundant edges (see \code{demo(paloma)} to have an
    example with \pkg{igraph} package).
  }
  \item{k}{
    is the node number of the searched motifs.   
  }
  \item{directed}{
    specifies if the graph is directed (or not).
  }
}



\value{
  The result is returned as a list.
  \describe{
    \item{}{
      \code{rv[[1]]} list concerning a motif size (the \emph{kth}).
    }
    \item{}{
      \code{rv[[1]]$k} motif node number. 
    }
    \item{}{
      \code{rv[[1]]$directed} directed motifs or not (logical). 
    }
   \item{}{
      \code{rv[[1]]$motifs$adjacency[[i]]} \emph{ith} adjacency matrix of
      an existing motif (canonical adjacency matrix).
    }    
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[i]]} list describing couples
      (\emph{ith} motif, deletion class).
    }
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[i]][[j]] } describes all
      \emph{positions} (k-1 vertices) in the graph with all the
      possible extensions (with one node) for
      the \emph{ith} motif and \emph{jth} deletion class.
    }
    \item{}{   
      \code{rv[[1]]$motifs$m.mp[[i]][[j]]$DelClass } nodes having the
      same deletion class (node indexes of the canonical motif) for the
      \emph{ith} motif and \emph{jth} deletion class.
    }
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[i]][[j]]$Occ[[l]]} \emph{lth} 
      \emph{position} (k-1 vertices) for the \emph{ith} motif with the \emph{jth}
      deletion class.
    }
    \item{}{
      \code{rv[[1]]$motifs$m.mp[[i]][[j]]$Ext[[l]]} all possible
      extensions of the \emph{lth} \emph{position}
  (\ldots{}\code{$m.mp[[i]][[j]]$Occ[[l]]} - containing k-1 vertices). 
    }
 }
 Summary :
    \tabular{ll}{
      rv[[k]] \tab motifs having the same motif size \cr
      rv[[k]]$k \tab the motif size value \cr
      rv[[k]]$directed \tab the motifs are directed (or not) \cr
      rv[[k]]$motifs \tab motifs list for the given size \cr
      rv[[k]]$motifs$adjacency[[i]] \tab \code{ith} motif adjacency matrix \cr
      rv[[k]]$motifs$m.mp[[i]] \tab occurrences and extensions for the
      \code{ith} motif \cr
      rv[[k]]$motifs$m.mp[[i]][[j]] \tab occurrences and extensions given
      a motif (\code{ith}) \cr
      \tab and a deletion class (\code{jth}) \cr
      
      rv[[k]]$motifs$m.mp[[i]][[j]]$DelClass \tab elements of the deletion class  \cr
      rv[[k]]$motifs$m.mp[[i]][[j]]$Occ[[l]] \tab \code{lth} \emph{position} (k-1 vertices) \cr  
      rv[[k]]$motifs$m.mp[[i]][[j]]$Ext[[l]] \tab possible extensions for the
      \emph{lth} \emph{position} \cr
    }
  
    Remark: for performance reasons, the node indexes have been
    renumbered in the
    \emph{position} (k-1 vertices) field (\code{Occ}) and the extension field
    (\code{Ext}). 
    To get initial numbering add the original minimum index.
  }

\author{G. Grasseau, E. Birmele}

\section{Warning }{
Setting a great value to the motifs size (\code{k > 7}) could generate
huge occurrence lists and thus cause memory allocation problems .
}


\seealso{
\code{\link{paloma-package}}, \code{\link{getExtensions}}
}

\examples{

## Load the data set
data(example)

# Get all occurrences and extensions of motifs
# containing 3 nodes in a graph.
out <- getAllExtensions( example$graph, 3, directed=TRUE )

# Number of nodes
out[[1]]$k

# Considering the second motif
out[[1]]$motifs$adjacency[[2]]
motif <- out[[1]]$motifs$m.mp[[2]]

# There is 3 deletion class for this motifs
# Select the third one
# D class (node list)
motif[[3]]$DelClass  # deletion class (node list)

# Second occurrence in this class (k-1 nodes)
motif[[3]]$Occur[[2]]

# Possible extensions of this occurrence
# in the same deletion class
motif[[3]]$Ext[[2]]
}

\keyword{ graphs }
