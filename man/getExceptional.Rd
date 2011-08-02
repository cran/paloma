\name{getExceptional}
\alias{getExceptional}
\title{
 Seeks for exceptional motifs under a statistical block model.
}  
\description{
  Extracts all motifs having \code{[kmin,kmax]} number of nodes which are
  over-represented according to the statistical block model.
  This statistical model assumes that the graph nodes belong
  to classes and the probability to connect two nodes belonging to two
  classes is described by the connectivity matrix \code{Pi[i,j]}
  (probability to connect a node belonging to the class \code{j} to a
  node belonging to the class \code{i}).
}
\usage{
getExceptional( G,
                kmin = 2, kmax = 2,
                NodeToClass = NULL, Pi = NULL,
                pvalue = 0.001,
                directed = FALSE,
                h.filter = TRUE
               )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{
    the graph is described as an edge list in a matrix
    \code{ edges[ nbr_of_edges, 2] }. The first column contains the source
    nodes, whereas the second one concerns the destination nodes. The node
    indexes must be contiguous and the graph must not contain
    self-loops nor redundant edges (see \code{demo(paloma)} to have an
    example with \pkg{igraph} package).
}
  \item{kmin}{
    minimum size (number of nodes) of the motifs.
}
  \item{kmax}{
    maximum size (number of nodes) of the motifs.
}
  \item{NodeToClass}{
    for each node gives its class index. These class indexes must be
    contiguous and the number of classes have to be equal to the \code{Pi}
    matrix dimensions.
    If  \code{ NodeToClass} and \code{Pi} are \code{NULL}, all the nodes are
    set to one class. The \code{Pi} value is computed from the graph.
  }
  \item{Pi}{
    is the connectivity matrix. \code{Pi[i,j]}
    is the probability to connect a node belonging to the class \code{j} to a
    node belonging to the class \code{i}.
    If  \code{ NodeToClass} and \code{Pi} are \code{NULL}, all the nodes are
    set to once class. The  \code{Pi} value is computed from the graph.
}
  \item{pvalue}{
    critical probability from which a motif with a deletion class is
     considered as an over-represented motif under the statistical
     model. Because \code{pvalue} is in fact an upper bound of the actual
     \emph{p-value}, it can be greater than 1.
     
  }
  \item{directed}{
    specifies if the graph is directed (or not).
}
  \item{h.filter}{
    specifies if the hierarchical filter is activated (or not)
  }
}
\details{
  The hierarchical filter \code{h.filter} is used to remove from the
  over-represented motifs the solutions 
  (motif of \code{k} nodes, deletion class) induced by an
  over-represented (motif of \code{k}' nodes, deletion class) solution
  with a smaller number of nodes (\code{k}'\code{<k}).
}
\value{
  The result is returned an \code{paloma.exceptional} object:
  \describe{
    \item{}{
      \code{rv[[k]]} list concerning a motif size (the \emph{kth} ).
    }
    \item{}{
      \code{rv[[k]]$k} motif node number. 
    }
    \item{}{
      \code{rv[[k]]$directed} the motifs are directed or not (logical). 
    }
    \item{}{
      \code{rv[[k]]$motifs$adjacency[[i]]} \emph{ith} adjacency matrix of
      an over-represented motif (canonical adjacency matrix).
      
    }
    \item{}{
      \code{rv[[k]]$motifs$m.mp[[i]]} list describing couples
      (\emph{ith} motif, deletion class).
    }
    \item{}{
      \code{rv[[k]]$motifs$m.mp[[i]][[j]]$Del} nodes indexes (of the canonical
      adjacency matrix) belonging to the same class. Select a couple
    (\emph{ith} motif, \emph{jth} deletion class)
    }
    \item{}{
      \code{rv[[k]]$motifs$m.mp[[i]][[j]]$Filter[1] = 0} specifies if the
      exceptionally of the couple (\emph{ith} motif, \emph{jth}
      deletion class) is induced by a smaller over-represented
      sub-motif. 
      Other values of Filter are for internal uses.
    }
    \item{}{
      \code{rv[[k]]$motifs$m.mp[[i]][[j]]$PValue} p-value of the couple
      (\emph{ith} motif, \emph{jth}
      deletion class) under the statistical block model. 
    }  
  
  }
  Summary :
  \tabular{ll}{
    rv[[k]] \tab motifs having the same motif size \cr
    rv[[k]]$k \tab the motif size value \cr
    rv[[k]]$directed \tab motifs are directed (or not) \cr
    rv[[k]]$motifs \tab motifs list for the given size \cr
    rv[[k]]$motifs$adjacency[[i]] \tab \code{ith} motif adjacency matrix \cr
    rv[[k]]$motifs$m.mp[[i]] \tab deletion class list for motif \code{i}\cr
    rv[[k]]$motifs$m.mp[[i]][[j]]$Del \tab \code{jth}  deletion class \cr
    rv[[k]]$motifs$m.mp[[i]][[j]]$Filter[1] \tab for \code{h.filter}   \cr
    rv[[k]]$motifs$m.mp[[i]][[j]]$PValue \tab p-value of the couple \cr
    \tab (\code{ith} motif, \code{jth} deletion class) \cr
    }
    
}
\author{G. Grasseau, E. Birmele}

\section{Warning }{
Setting a great value to the motifs size could generate
huge occurrence lists and thus cause memory allocation problems.
}

\seealso{
\code{\link{paloma-package}}, \code{\link{getExtensions}}
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

# The same without the hierarchical filter
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001, h.filter=FALSE )

# Print exceptional motifs
out

# Without classes
out <- getExceptional( example$graph, 3, 5, pvalue=0.001)

}
\keyword{ graphs }
\keyword{ models }
\keyword{ cluster }
