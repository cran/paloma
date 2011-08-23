\name{paloma-package}
\alias{paloma-package}
\alias{paloma}
\docType{package}
\title{
  Discovering motifs in networks which are locally over-represented. 
}
\description{
  The \pkg{paloma package} deals with  discovering motifs in networks which are
  locally over-represented. All subgraphs considered are induced subgraphs.

  The main idea is to consider a motif \emph{m} (small graph) and a
  sub-motif \emph{m'} of \emph{m} obtained by deleting a vertex of
  \emph{m}. 
  An occurrence of \emph{m'} in the network of
  interest is called a \emph{position}. An extension of a position is
  a vertex which addition yields an occurrence of \emph{m}. 
  For each position, the number of extensions is computed and compared
  to its mean in a random graph model.

  As shown in [Bir2010], a p-value on a statistic equivalent to the maximal ratio (observed
  extensions)/(expected extensions) over all positions
  can be computed using Poisson approximations.
  It allows to detect motifs that are locally over-represented with
  respect to one of their sub-motifs somewhere in the network.
  

  The random graph model used is the block model, defined by a set of
  vertex classes (blocks) and a connectivity matrix \code{Pi}. Every
  vertex belongs to one class and  the
  probability of existence of an edge between a
  vertex of class \code{i} and a vertex of class \code{j} is given by \code{Pi[j,i]}.
  This family of random graph models contains the classical Erdos-Renyi
  model where all edges are drawn with probability \code{p} by taking only
  one class of vertices. 

  Some vertices of \emph{m} may be isomorphic, and thus lead to the same
  sub-motif occurrences and extensions. Such vertices are grouped in
  classes called \emph{deletion classes}. Every motif \emph{m} is then
  tested to be over-represented or not with respect to each of its
  deletion classes.

  
 The main functions in this package are:
 \itemize{
   \item \code{getExceptional} computes the motifs which are locally
   over-represented in a network. Its main arguments are the network of
   interest, the motif sizes to be explored and the p-value threshold.

   The default random
 graph model is the Erdos-Renyi one but better classes and
 connectivities can be given.

 A hierarchical filtering procedure allows to
 avoid repetitions in the final list. A local over-representation of a
 motif of size 4 may for example be redundant with a local over-representation of a
 motif of size 3 which was already detected (see [Bir2010]).

 
 \item \code{getAllExtensions} and \code{getExtensions} list the
 positions and extensions for a deeper view on the local over-representation.

 }

 Some functions are added to simplify the use of the former functions:
 \itemize{
   \item \code{print.exceptional.motif} returns a simplified printing of
  the results of \code{getExceptional}.
   \item \code{plot.paloma.exceptional} plots the results of
 \code{getExceptional}. Every locally over-represented
 motif is shown with its p-value and the nodes of the considered deletion
 class are colored in red.
   \item \code{print.paloma.occurrences} returns  a simplified printing
   of the results of \code{getAllExtensions} or \code{getExtensions}.
 }
 
 


 The package can deal with directed and undirected networks and two
 example are provided:
 \itemize{
   \item \code{demo(read_graphs)}
   \item \code{demo(yeast)}
  }

  Due to large numbers of occurrences, the main functions
  (\code{getAllExtensions} and \code{getExtensions}) are time CPU and memory
  consuming. For a network with 1838 nodes and 5928 edges, requiring
  exceptional motifs from 2 nodes to 4 nodes takes about 600 seconds. 

}


\section{Acknowledgement}{
This work has been supported by the French Agence Nationale de la
Recherche under grant NeMo ANR-08-BLAN-0304-01.
}

\author{
\itemize{
    \item Etienne Birmele \email{Etienne.Birmele@genopole.cnrs.fr},
    \item Gilles Grasseau \email{Gilles.Grasseau@genopole.cnrs.fr},
  }
}

\references{
[Bir2010] Birmele E.,
\emph{Detection of local network motifs.}
arXiv preprint.
\url{http://arxiv.org/pdf/1007.1410}




  [FanMod] S. Wernicke and F. Rasche. FANMOD: a tool for fast network motif detection.
  Bioinformatics, 22(9):1152-1153, 2006.
  
  [FanMod] S. Wernicke.
  Efficient detection of network motifs.
  IEEE/ACM Transactions on Computational Biology and Bioinformatics,
  3(4):347-359, 2006.
  
  [Nauty] http://cs.anu.edu.au/~bdm/nauty/
  

}
\keyword{ package }

