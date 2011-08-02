\name{countAllMotifs}
\alias{countAllMotifs}

\title{
Counts the number of all motifs having k nodes in a given graph
}

\description{
 counts all induced motifs for a given size in a graph. The algorithm
  used is the one used in the \code{FanMod} code (see references in
  \code{\link{paloma-package}}). The motif counters are updated
  according to the canonical
  form of the found motifs (the canonical is computed with 
  the Nauty algorithm - see the \code{\link{paloma-package}} references).
}

\usage{
countAllMotifs( G, k, directed = FALSE )
}

\arguments{
  \item{G}{
    the graph is described as an edge list in a matrix
    \code{ edges[ nbr_of_edges, 2] }. The first column contains the source
    nodes, whereas the second one maps the destination nodes. The node
    indexes must be contiguous and the graph must not contain
    self-loops nor redundant edges (see \code{demo(paloma)} to have an
    example with \pkg{igraph} package).
  }
  \item{k}{
    is the node number of the required motifs.
  }
  \item{directed}{
    specifies if the graph is directed (or not).
  }
}

\value{
  Returns a list containing an adjacency matrix list
  (\code{output$adjacency}) and their respective counts
  (\code{output$count}).
}

\author{G. Grasseau}


\seealso{
\code{\link{paloma-package}}
}

\examples{
## Load the data set
data(example)

# Count all motifs with 3 nodes in a graph
out <- countAllMotifs( example$graph, 3, directed=TRUE )

# Number of discovery motifs
motif_nbr <- length( out$adjacency)

# Getting the fifth motif
motif <- matrix( out$adjacency[[5]], 3, 3  )

# Number of occurrences of this motif
nb <- out$count[5]
}

\keyword{ graphs }

