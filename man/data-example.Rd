\name{example}
\alias{example}

\docType{data}

\title{Graph example for package tests} 

\usage{data(example)}

\description{%
  Data set to use or test \code{paloma}. The first field contains the
  edge list of the network. The second one contains the
  class indexes for each node (there are 2 classes). And the last one
  is the connectivity matrix between the two classes.
}

\examples{
## load the data set
data(example)

out <- getExceptional( example$graph, 3, 4,
                       example$classes, example$pi, 0.001 )
}

\keyword{datasets}
  