\name{computModelMean}
\alias{computeModelMean}
\title{Computes the expected occurrence number of a motif}

\description{
  For a given statistical block model, computes the expected occurrence
  counts of a motif. This statistical model is described by the number
  of node classes, the node belonging to classes and the connectivity
  matrix \code{Pi[i,j]}
  (probability to connect a node belonging to the class \code{j} to a
  node belonging to the class \code{i}).
}
\usage{
computeModelMean( NodeToClass, Pi, motif, directed = FALSE )
}
\arguments{
  \item{NodeToClass}{
    for each node gives its class index. These class indexes must be
    contiguous and the number of classes have to be egal to the \code{Pi}
    matrix dimensions.
  }
  
  \item{Pi}{
    is the connectivity matrix. \code{Pi[i,j]}
    is the probability to connect a node belonging to the class \code{j} to a
    node belonging to the class \code{i}.
  }

  \item{motif}{
    is the adjacency matrix of the motif. \code{motif} can be vector or
    matrix type.
  }

  \item{directed}{
    specifies if the graph is directed (or not).
  }
}

\value{
  returns the expected number of motif occurrences for the given model.
}
\author{G. Grasseau, E. Birmele}

\seealso{
  \code{\link{paloma-package}} and \code{\link{getExceptional}}.
}
  
\examples{
## Load the data set
data(example)

# Get over-represented motifs under a statistical
# block model with two classes
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )


# Get the first over-represented motif
m  <- out[[1]]$motifs$adjacency[[1]]

# Compute the expected counts of the motif occurrences
computeModelMean( example$classes, example$pi, m)

}
\keyword{ stats }
\keyword{ models }
\keyword{ cluster }