# include <iostream>

# include "math.h"
# include "nr-permutation.h"
# include "sparse-adjacency.h"
# include "combinatorics.h"
# include "util.h"

using namespace std;

class RandomGraphModel {
 
 protected:
  int64_t N_;

 public:
  
  virtual double computeOccurenceProbability( const SparseAdjacency &motif) = 0;
  double computeMean( const SparseAdjacency   &motif, 
		      NonRedondantPermutation &NRPerm );
  // Two order momentum 
  double computeMoment2( const SparseAdjacency   &motif, 
			 NonRedondantPermutation &NRPerm );
};

class ErdosRenyi : public RandomGraphModel
{

 protected:
  bool directed_;
  double pi_;

 public:
  ErdosRenyi( int64_t N, double pi, bool directed ) {N_ = N, pi_ = pi; 
    directed_ = directed; };
  double computeOccurenceProbability( const SparseAdjacency &motif );
};

// Expected Degree Distribution model
class EDD: public RandomGraphModel  {

 private:
  bool directed_;
  double *DegreeDistribution_;
  int NbrDegree_;
  double *Moments_;
  int NbrMoments_;
  double gamma_;

 public:
  EDD( const int64_t N, const double *DegreeDistribution, const int NbrDegree, bool directed );
  void setMoments( int max_degree );
  double computeOccurenceProbability( const SparseAdjacency &motif ) ;
  ~EDD( );
};

class Mixnet: public RandomGraphModel  {

 private:
  bool directed_ ;
  double *pi_;
  double *alpha_;
  int Q_;

 public:
  Mixnet( const int64_t  N, const int Q, 
	  const double *alpha, const double *pi, 
	  bool directed );
  double computeOccurenceProbability( const SparseAdjacency &motif ) ;
  ~Mixnet( );
};

// Etienne B. Model
class PseudoMixnet: public RandomGraphModel  {

 private:
  bool directed_ ;
  double *pi_;
  // Number of nodes belonging to class i
  int *R_;
  int Q_;

 public:
  PseudoMixnet( const int64_t  N, const int Q, 
	  const int *R, const double *pi, 
	  bool directed );
  double computeOccurenceProbability( const SparseAdjacency &motif ) ;
  double computeLambdaU( const SparseAdjacency & motif, int s, 
			 int *NodeToClass );
  double computeLambdaU( const SparseAdjacency & G, 
			 int *NodeList, 
			 int NodeListSize, 
			 int NodeRemoved, 
			 int *NodeToClass );
  ~PseudoMixnet( );
};
