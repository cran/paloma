# include <iostream>
# include <iomanip>

# include "graphical-models.h"

using namespace std;

EDD::EDD( const int64_t N, const double *DegreeDistribution, const int NbrDegree, bool directed ) {

  N_ = N;
  NbrDegree_ =  NbrDegree;
  directed_ = directed;

  if ( directed_ ) {
    DegreeDistribution_ = new double[2*NbrDegree_];
  } else {
    DegreeDistribution_ = new double[NbrDegree_];
  }

  // Undirected case or Incoming arcs
  double sum = 0.0;
  for ( int i = 0; i < NbrDegree_; i++) {
    sum += DegreeDistribution[i];
  }
  if( sum == 0)  
    LOGMSG( 1, std::cout, "Bad degree distribution", "" );
  for ( int i = 0; i < NbrDegree_; i++) {
    DegreeDistribution_[i] = DegreeDistribution[i]/sum;
  }

  if (directed_) {
    //
    // Degree distribution for Outcoming arcs
    sum = 0;
    for ( int i = NbrDegree_ ; i < 2*NbrDegree_; i++) {
      sum += DegreeDistribution[i];
    }
    if( sum == 0)  
      LOGMSG( 1, std::cout, "Bad degree distribution", "" );
    for ( int i = NbrDegree_; i < 2*NbrDegree_; i++) {
      DegreeDistribution_[i] = DegreeDistribution[i]/sum;
    }
  }
  // Set the 10 first moments
  NbrMoments_ = 0;
  Moments_    = 0;
  setMoments( 10 );

  // Compute gamma
  // For directed case
  // GammaInt = GammaOut = Gamma
  gamma_ = 1.0 / (Moments_[1] * ( N_ - 1) ); 

 // Gamma normalisation 
 // constraint : d_max^2 <= (n-1) E(D)
 double var = 1.0 / ( (NbrDegree_ - 1.0)*(NbrDegree_ - 1.0) );
 gamma_ = MIN( gamma_ , var );
}

EDD::~EDD( ) {

  delete [] DegreeDistribution_ ;
  delete [] Moments_ ;
}

Mixnet::Mixnet( const int64_t N, const int Q, const double *alpha, const double *pi, bool directed ) {

  // Note
  // pi[i,j] index of line i is countiguous in memory( column storage mode)
  // pi[i,j] is the proba. to conect j -> i in direct case.
 
  Q_ =  Q;
  N_ =  N;
  directed_ = directed;
  alpha_ = new double[Q];

  pi_ = new double[Q*Q];

  for ( int i = 0; i < Q_; i++) {
    alpha_[i] = alpha[i];
  }

  for ( int i = 0; i < Q_*Q_; i++) {
    pi_[i] = pi[i];
  }
}

Mixnet::~Mixnet( ) {

  delete [] alpha_ ;
  delete [] pi_ ;
}

PseudoMixnet::PseudoMixnet( const int64_t N, const int Q, const int *R, const double *pi, bool directed ) {

  // Note
  // pi[i,j] index of line i is countiguous in memory( column storage mode)
  // pi[i,j] is the proba. to conect j -> i in direct case.
 
  Q_ =  Q;
  N_ =  N;
  directed_ = directed;
  R_ = new int[Q];

  pi_ = new double[Q*Q];

  for ( int i = 0; i < Q_; i++) {
    R_[i] = R[i];
  }

  for ( int i = 0; i < Q_*Q_; i++) {
    pi_[i] = pi[i];
  }
}

PseudoMixnet::~PseudoMixnet( ) {

  delete [] R_ ;
  delete [] pi_ ;
}

double
ErdosRenyi::computeOccurenceProbability( const SparseAdjacency & motif ) {

  // Compute m++/2 for symetric matrix or
  // m++ for directed matrix
  int sumMatrix = motif.getNbrOfValues();
  double res = pow( pi_, sumMatrix );
  return res;
}

void 
EDD::setMoments( int max_degree ) {

  if ( (max_degree+1) > (NbrMoments_) ) {
    if ( Moments_ != 0 ) delete [] Moments_;
    NbrMoments_ = max_degree + 1;

    // Compute distribution moments
    if( directed_ )
      Moments_ = new double[2*NbrMoments_];
    else
      Moments_ = new double[NbrMoments_];

    // Undirected case  or Incoming moments
    Moments_[0] = 1.0;
    for(int m=1; m < NbrMoments_ ; m++) {
      Moments_[m] = 0.0;
      for(int i=0; i < NbrDegree_; i++) {
	Moments_[m] += DegreeDistribution_[i] * pow ( static_cast<double>(i), 
						      static_cast<double>(m) );
      }
    }
    if (directed_) {

      // Outcoming moments
      Moments_[ NbrMoments_ ] = 1.0;
      for(int m = NbrMoments_+1; m < 2*NbrMoments_ ; m++) {
	Moments_[m] = 0.0;
	for(int i = NbrDegree_; i < 2*NbrDegree_; i++) {
	  Moments_[m] += DegreeDistribution_[i] * 
	    pow ( static_cast<double> (i - NbrDegree_), 
		  static_cast<double> (m - NbrMoments_) );
	}
      }
    }
  }
}

double
EDD::computeOccurenceProbability( const SparseAdjacency & motif ) {

  if ( directed_ == motif.getSymmetry()) {
    LOGMSG( 1, std::cout, "Model and motif cant be compared (symmetry)", "" );
    return 0.0;
  }
  // Compute m++/2 for undirected matrix
  // or m++ for directed matrix
  int sumMatrix = motif.getNbrOfValues();
  double res = pow( gamma_, sumMatrix );

  if ( directed_ ) {

    const int *InDegrees  = motif.getRowSums( );
    const int *OutDegrees = motif.getColSums( );

    // Get motif degree max
    int64_t max_degree =  0;
    for(int u=0; u < motif.getNbrRow(); u++) {
      max_degree = MAX( InDegrees[u], max_degree);
    } 
    for(int u=0; u < motif.getNbrRow(); u++) {
      max_degree = MAX( OutDegrees[u], max_degree);
    } 
    setMoments( max_degree );
    
    // Compute Occurence probability mu
    for(int u=0; u < motif.getNbrRow(); u++) {
      // Incoming
      res = res * Moments_[ InDegrees[u] ];
      // Outcoming
      res = res * Moments_[ OutDegrees[u] + NbrMoments_ ];
    } 
    // delete [] InDegrees;
    // delete [] OutDegrees;

    return( res );
  } else {
    const int *degrees = motif.getMatrixRowSums( );

    // Get motif degree max
    int64_t max_degree =  0;
    for(int u=0; u < motif.getNbrRow(); u++) {
      max_degree = MAX( degrees[u], max_degree);
    } 
    setMoments( max_degree );
    
    // Compute Occurence probability mu
    for(int u=0; u < motif.getNbrRow(); u++) {
      res = res * Moments_[ degrees[u] ];
    } 
    return( res );
  }
}


double
Mixnet::computeOccurenceProbability( const SparseAdjacency & motif ) {

  // Loop over all motif configurations  (c1, c2, ..., ck; 0 <= ci < Q )
    //  for(i=0; i < k; i++)
    //   term = term * alpha[conf[i]]

    // The matrix contains only the upper-triagular block
    //  for all edges (i, j)
    //  term = term * pi[conf[i), conf[j])
    // sum = sum + 
    
  int nnodes = motif.getNbrRow();

  // Create an configuratio iterator with vector of number nodes with Q_ states. 
  ConfigurationIterator ConfIter( Q_ , nnodes );

  // Reference to read the configuration content 
  const int *conf = ConfIter.restart();

  double mu   = 0.0;
  double term = 1.0;

  // Iterate configuration
  do {
    term = 1.0;
    for ( int i=0; i < nnodes; i++) {
      term = term * alpha_[conf[i]];
    }
    for ( int i=0; i < nnodes; i++) {
      // General case : directed graph
      if (directed_) {
	for ( int j=0; j < nnodes; j++) {
	  if( motif.getMatrix(i,j) ) {
	    term = term * pi_[ conf[i] + Q_ *conf[j]];
	  }
	}
      } else {
	for ( int j=i+1; j < nnodes; j++) {
	  if( motif.getMatrix(i,j) ) {
	    term = term * pi_[ conf[i] + Q_ *conf[j]];
	  }
	}
      }
    }
    mu += term;
  } while ( ConfIter.iterate() );
  return mu;
}

double
PseudoMixnet::computeOccurenceProbability( const SparseAdjacency & motif ) {

  // Loop over all motif configurations  (c1, c2, ..., ck; 0 <= ci < Q )
    //  for(i=0; i < k; i++)
    //   term = term * alpha[conf[i]]

    // The matrix contains only the upper-triagular block
    //  for all edges (i, j)
    //  term = term * pi[conf[i), conf[j])
    // sum = sum + 
    
  int nnodes = motif.getNbrRow();

  // Create an configuratio iterator with vector of number nodes with Q_ states. 
  ConfigurationIterator ConfIter( Q_ , nnodes );

  // Reference to read the configuration content 
  const int *conf = ConfIter.restart();

  double mu   = 0.0;
  double term = 1.0;

  // Counter to compute A(n,m) with R_[i] number of nodes belonging to class i
  int *d = new int[Q_]; 

  // init
  // Iterate configuration
  do {
    term = 1.0;
    // Init counters d
    for( int *it=d, *it_end=&d[Q_]; it != it_end; it++){
      *it = 0;
    }

    for ( int i=0; i < nnodes; i++) {
      // Equivalent to alpha product 
      term = term * ( R_[ conf[i] ] - d[ conf[i] ] );
      d[ conf[i] ]  =  d[ conf[i] ] + 1;
      if ( term <= 0 ) break;
    }
    
    // printList("Config. :", conf, &conf[nnodes], true);
    // cout << "  Number  : " << term << endl;

    if(  term > 0 ) {
      // General case : directed graph
      if (directed_) {
	for ( int i=0; i < nnodes; i++) {
	  for ( int j=0; j < nnodes; j++) {
	    //
	    // Loop not taken into account
	    if( i != j ) {
	      if( motif.getMatrix(i,j) ) {
		term = term * pi_[ conf[i] + Q_ * conf[j]];
	      } else {
		term = term * ( 1 - pi_[ conf[i] + Q_ * conf[j] ]);
	      }
	    }
	  }
	}
      } else {
	for ( int i=0; i < nnodes; i++) {
	  for ( int j=i+1; j < nnodes; j++) {
	    if( motif.getMatrix(i,j) ) {
	      term = term * pi_[ conf[i] + Q_ * conf[j]];
	    } else {
	      term = term * ( 1 - pi_[ conf[i] + Q_ * conf[j]]);
	    }
 	  }
	}
      }
      // cout << "  Contrib  : " << term << endl;
      mu += term;
    }
  } while ( ConfIter.iterate() );
  delete [] d;
  // return mu / factorial( nnodes, 1L );
  return mu / (BinomCoef( N_, nnodes ) * factorial( nnodes, 1L) );
}

double
RandomGraphModel::computeMean( const SparseAdjacency & motif, 
			       NonRedondantPermutation &NRPerm ) {

  double res;
  int64_t k = motif.getNbrRow();

  int64_t CardNRPermutation = NRPerm.getCardNRPermutation();

  double coef = BinomCoef( N_, k );
  double card = CardNRPermutation;
  double proba = computeOccurenceProbability( motif );
  res = coef * card * proba;
  
  // cout << "coef, card, proba, mean : " 
  //      << coef << " "
  //      << card << " "
  //      << proba << " "
  //      << res << " "
  //      << endl;

  return res ;  

}

double
RandomGraphModel::computeMoment2( const SparseAdjacency & motif, 
				  NonRedondantPermutation &NRPerm0) {

  SparseAdjacency motif0 (motif);
  SparseAdjacency motif1 (motif);

  NonRedondantPermutation NRPerm1(NRPerm0);
  double res;
  int k = motif.getNbrRow();
  int *inv_perm1 = new int[k];
  int *dperm1 = new int[k];

  int64_t CardNRPermutation = NRPerm0.getCardNRPermutation();

  // TODO OPT Already computed
  double var = CardNRPermutation * computeOccurenceProbability( motif );
  res = MultinomCoef( N_, (N_ - 2*k), k, k ) * var * var;

  /* 
  cout << "computeMoment2 " << MultinomCoef( N_, (N_ - 2*k), k, k ) << endl;
  cout << "N_, (N_ - 2*k), k, k )" 
       << N_ << " "
       << (N_ - 2*k) << " "
       << k << " "
       << k << " "
       << endl;
  */

  double coef;
  double MuSum;
  double OverlappingSum = 0;
  for ( int s=1; s <= k; s++) {
    if ( N_ >= (2*k -s) ) {
      coef = MultinomCoef( N_, k-s, s, k-s, N_ -2*k +s );
      /*  
	 cout << "computeMoment2 super " << coef 
	 << endl;
	 cout << "  N_, k-s, s, k-s, N_ -2*k +s" 
	 << N_ << " "
	 << (k-s) << " "
	 << s << " "
	 << (k-s) << " "
	 <<   N_ -2*k +s << " "
	 << endl;
      */
      
      MuSum = 0;
      
      //
      // Compute mu of the super-motif over R(m) x R(m)
      //
      
      // TODO (perf) could be done with only ij_index  
      const int *perm0 = NRPerm0.iterate("start");
      for ( ; perm0 != 0; ) {
	//
	// Build the lower part of the super-motif
	// with the motif
	// Called super-pattern
	// 
	SparseAdjacency super_pattern( motif, 2*k-s);
	super_pattern.permute( perm0, k );
	
	// Shift the nodes in the super-pattern
	// to build the upper part : overlap (s) and motif remaining (k-s)
	//                          start, size, shift_value
	super_pattern.shiftNumbering( 0, k, k-s );
	
	const int *perm1 = NRPerm1.iterate( "start" );
	// Determining inv( perm1 )
	for ( int i = 0; i < k; i++){
	  inv_perm1[ perm1[ i ]] = i;
	} 

	for ( ; perm1 != 0; ) {

	  // The resulting permutation is perm1(k) o inv_perm1(k_1)
	  for ( int i = 0; i < k; i++){
	    dperm1[ i ] = perm1[ inv_perm1[i] ];
	  } 
	  
	  // Determining inv( perm1 )
	  for ( int i = 0; i < k; i++){
	    inv_perm1[ perm1[ i ]] = i;
	  } 

	  SparseAdjacency super_motif( super_pattern );
	  
	  // Permute nodes
	  motif1.permute( dperm1 );
	  
	  // super_motif.print(" super-pattern (avant)");
	  // motif1.print("motif1 ");
	  
	  // Duplicate permuted matrix and finish building super-motif
	  super_motif.setMatrix( motif1 );
	  // super_motif.print(" super-motif");
	  
	  
	  // Compute mu( super-motif)  and update mu-sum
	  MuSum += computeOccurenceProbability( super_motif );
	  
	  perm1 = NRPerm1.iterate();
	}
	perm0 = NRPerm0.iterate();
      }
      OverlappingSum += MuSum * coef;
    }
  }

  delete [] inv_perm1;
  delete [] dperm1;

  res += OverlappingSum;
  return  res;  

}

double
PseudoMixnet::computeLambdaU( const SparseAdjacency & G, 
			      int *NodeList, 
			      int NodeListSize, 
			      int RemovedNode, 
			      int *NodeToClass ) {
  // s node removed 

  // Loop over all G configurations  (c1, c2, ..., ck; 0 <= ci < Q )
    //  for(i=0; i < k; i++)
    //   term = term * alpha[conf[i]]

    // The matrix contains only the upper-triagular block
    //  for all edges (i, j)
    //  term = term * pi[conf[i), conf[j])
    // sum = sum + 
  

  int s =RemovedNode;

  double mu   = 0.0;
  double term = 1.0;

  // Counter to compute A(n,m) with R_[i] number of nodes belonging to class i
  int *d = new int[Q_]; 

  // init


  // Init counters d : number of G nodes in classe i
  for( int *it=d, *it_end=&d[Q_]; it != it_end; it++){
    *it = 0;
  }

  for( int *it=NodeList, *it_end = &NodeList[NodeListSize]; 
       it != it_end; it++) {
    d[ NodeToClass[ *it ] ]++;
  }


  // int ClassOfs =  NodeToClass[s];
  // cout << endl;
  // cout << "    Suppressed node : " << setw(2) << s 
  //     << ", initial class : " << ClassOfs <<   endl;


  // Iterate configuration
  for ( int q=0; q < Q_; q++) {
    term = R_[ q ] - d[ q ];
    if ( term > 0 ) { 
    
      // printList("Config. :", conf, &conf[nnodes], true);
      // cout << "  Number  : " << term << endl;
    
      // General case : directed graph
      // cout << "    Class       : " << setw(3) << q << endl;
      // cout << "    Class Coef. : " << setiosflags(ios::fixed) 
      // << setprecision(2) << setw(6) << term << endl;
    
      int i, j ;
      if (directed_) {
	for ( int ii=0; ii < NodeListSize; ii++) {
	
	  i = NodeList[ii];
	  if( G.getMatrix(i,s) ) {
	    term = term * pi_[ NodeToClass[i] + Q_ * q];
	    // cout << "      node " << i 
	    //     << ":     Pi["  << setw(2) << NodeToClass[i] << "," 
	    //     <<  setw(2) << q << "]" 
	    //     << ", cumul: " << setw(6) << term << endl;
	  } else {
	    term = term * ( 1 - pi_[  NodeToClass[i] + Q_ * q ]);
	    // cout << "      node " << i << ": 1 - Pi["  << setw(2) << NodeToClass[i] << "," <<  setw(2) << q << "]" << endl;
	  }
	}
	for ( int jj=0; jj < NodeListSize; jj++) {
	  // Loop not taken into account
	  j = NodeList[jj];
	  if( G.getMatrix(s, j) ) {
	    term = term * pi_[ q + Q_ * NodeToClass[j] ];
	    // cout << "      node " << j 
	    //     << ":     Pi["  << setw(2) << NodeToClass[j] << "," 
	    //     <<  setw(2) << q << "]" 
	    //     << ", cumul: " << setw(6) << term << endl;
	  } else {
	    term = term * ( 1 - pi_[ q + Q_ * NodeToClass[j] ]);
	    // cout << "      node " << j << ": 1 - Pi["  << setw(2) 
	    // << NodeToClass[j] << "," <<  setw(2) << q << "]" << endl;
	  }
	}	  
      } else {
	for ( int ii=0; ii < NodeListSize; ii++) {
	  // Loop not taken into account
	  i = NodeList[ii];
	  if( G.getMatrix(i,s) ) {
	    term = term * pi_[ NodeToClass[i] + Q_ * q];
	    // cout << "      node " << i 
	    //     << ":     Pi["  << setw(2) << NodeToClass[i] << "," 
	    //     <<  setw(2) << q << "]" 
	    //     << ", cumul: " << setw(6) << term << endl;
	  } else {
	    term = term * ( 1 - pi_[  NodeToClass[i] + Q_ * q ]);
	    // cout << "      node " << i << ": 1 - Pi["  << setw(2) 
	    // << NodeToClass[i] << "," <<  setw(2) << q << "]" << endl;
	  }
	}
      }
    }
    // cout << "  Contrib  : " << term << endl;
    mu += term;
  }

  // cout << "    lamdaU = " << setw(6) << mu << endl;
  // cout << endl;

  delete [] d;

  // return mu / factorial( nnodes, 1L );
  return mu;
}

