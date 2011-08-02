\name{print.paloma.exceptional}
\alias{print.paloma.exceptional}

\title{Prints the over-represented motifs list}

\description{Displays the outputs of functions
  returning an "paloma.exceptional" object \code{\link{getAllExtensions}}
}

\usage{%
\method{print}{paloma.exceptional}(x, k=NULL, m=NULL, mp=NULL, ...)
}

\arguments{
  \item{x}{
    "paloma.exceptional" class output of  \code{getExceptional}
  }

  \item{k}{%
    deals with the motif size, selects a subset of the initial list \code{x}.
  }
  
  \item{m}{%
    selects a subset in the initial list corresponding to a subset of motifs.
  }

  \item{mp}{%
    selects a subset of deletion class (or extension type).
  }

  \item{\dots}{
    Additional arguments for generic \code{print}
    (such as \code{main ="title"}).
  }
}

\details{
  Displays the different levels of the occurrences list :
  \describe{
    \item{\bold{motif size}}{
      in an "paloma.exceptional" list
      different motif sizes can be stored. Display the
      all the couple (motif, deletion class) of the \emph{kth} motif size.
    }
    \item{\bold{motif}}{
      Display the all the the couple (motif, deletion class) of the \emph{mth} motif.
    }
    \item{\bold{deletion class}}{
         Display the all the the couple (motif, deletion class) of the
  \emph{mpth} deletion class.
   }
  }
}
\author{G. Grasseau}

\seealso{ \code{\link{getExceptional}}, \code{\link{paloma-package}}.}

\examples{

## Load the data set
data(example)

# Get couple (motif, deletion class) over-represented
# under a statistical block model with two classes
# The  the connectivity matrix is store in "example$pi"
# and the mapping node/class is store in "example$classes"
out <- getExceptional( example$graph, 3, 5, example$classes,
          example$pi, pvalue=0.001 )

# Get the all the over-represented motifs (motif, deletion class)
print( out )
# or
out

# Display  the couple (motif, deletion class) for the first motif
print( out, k=1, m=1 )

# The same with only the second deletion class
print( out, k=1, m=1, mp=2 )

}

\keyword{print}
