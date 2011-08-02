\name{plot.paloma.exceptional}
\alias{plot.paloma.exceptional}

\title{Plots the over-represented motif list}

\description{Plots the over-represented motifs founds 
  by the \code{\link{getExceptional}} function.
}

\usage{%
\method{plot}{paloma.exceptional}( x, ... )
}

\arguments{
  \item{x}{
    "paloma.exceptional" class output of  \code{getExceptional}
  }

  \item{\dots}{
    Additional arguments for generic \code{plot}
    (such as \code{main ="title"}).
  }
}

\details{
  Exceptional motifs are displayed according to the motif size order.
  The nodes belonging to the deletion class are red.
}

\author{G. Grasseau}

\seealso{ \code{\link{getExceptional}}, \code{\link{paloma-package}}.}

\examples{

## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The  the connectivity matrix is stored in "example$pi"
# and the mapping node/class is stored in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# Get the all the over-represented motifs (motif, deletion class)
plot( out )

}

\keyword{plot}
