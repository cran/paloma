\name{getCanonic}
\alias{getCanonic}

\title{
  Compute the canonical adjacency matrix
}

\description{
  Returns a unique adjacency matrix called canonical form. This involves
  that two isomorphic adjacency matrices have exactly the same canonical
  form. This routine uses the \code{nauty} library.
}

\usage{
getCanonic( m, directed = FALSE)
}

\arguments{
  \item{m}{
    adjacency matrix, it could be a vector or matrix type.
  }
  \item{directed}{
    specifies if the adjacency matrix is directed or not.
  }
}

\value{

  Returns the canonical form as a vector. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Details for \code{nauty} library can be found at
  \url{http://cs.anu.edu.au/~bdm/nauty}
}

\author{G. Grasseau}

\seealso{
\code{\link{paloma-package}}.
}

\examples{
# Define a 3 node motif
motifV <-  matrix( c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3) 

# Compute canonical form
cf <- getCanonic( motifV )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphs }
