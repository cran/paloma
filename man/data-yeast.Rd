\name{yeast}
\alias{yeast}

\docType{data}

\title{Graph example for package tests} 

\usage{data(yeast)}

\description{%
  This data set describes the Yeast regulatory network. The
  first field contains the edge list of the network (194 vertices with
  166 edges). The second one contains the
  class indexes for each node (there are 8 classes). The last one
  is the connectivity matrix between the eight classes. The vertex
  classification has been performed by the \pkg{mixer} package.
}

\examples{
## load the data set
data(yeast)

# Get exceptional motifs
out <- getExceptional( yeast$graph, 3, 4,
                       yeast$classes, yeast$pi, 0.001, directed=TRUE )

# Show the exceptional motifs
plot( out )

}

\keyword{datasets}
  