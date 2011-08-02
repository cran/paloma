 \name{print.paloma.occurrences}
\alias{print.paloma.occurrences}

\title{Prints the motifs occurrences and their extensions}

\description{Displays the outputs of functions
  returning an "paloma.occurrences" object (\code{\link{getAllExtensions}}  and \code{\link{getExtensions}}).
}

\usage{%
\method{print}{paloma.occurrences}( x, k=NULL, m=NULL, mp=NULL, ...)
}

\arguments{
  \item{x}{
    "paloma.occurrences" class output of \code{getAllExtensions}
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
      in an "paloma.occurrences" list
      different motif sizes can be stored. Displays 
      all of the occurrences and extensions of the \emph{kth} motif size.
    }
    \item{\bold{motif}}{
      Display the all the the occurrences and extensions of the \emph{mth} motif.
    }
    \item{\bold{deletion class}}{
      For this deletion class all the occurrences and extensions are displayed.
   }
  }
}
\author{G. Grasseau}

\seealso{ \code{\link{getAllExtensions}},
  \code{\link{getExtensions}}, \code{\link{paloma-package}}.}

\examples{
## Load data set
data( example )

# Get the all the extensions for a motif size 
occ <- getAllExtensions( example$graph, 3, directed=TRUE )

# Display extensions (occurrences)
print( occ )
# or
occ

# Display the occurrences and extensions of the first motif
# Two deletion classes are displayed 
print( occ, k=1, m=1 )

# Display the occurrences and extensions of the second deletion class of
# the first motif  
print( occ, k=1, m=1, mp=2 )

}

\keyword{print}
