# include "nr-permutation.h"

NonRedondantPermutation::NonRedondantPermutation( const SparseAdjacency & motif ) {
  
  int NbCol =  motif.getNbrCol();
  int NbRow =  motif.getNbrRow();
  bool MotifSymetry = motif.getSymmetry();

  int *CurrentIndexes = new int[ NbRow ];

  if (NbCol == NbRow) {
    Index_size_ =  NbRow;
    for ( int i= 0; i < Index_size_; i++) { CurrentIndexes[i] = i; }

    //
    // Generate all permutations and automorphisms;
    //
    Permutations_  = new list<int*>;
    Automorphisms_ = new list<int*>;

    // Put identity in Permutation list and Automorphism list
    int k = 0;
    k = k+1;
    copy_add_list_( Permutations_ , CurrentIndexes );
    copy_add_list_( Automorphisms_, CurrentIndexes );

    for( ; next_permutation( CurrentIndexes, CurrentIndexes + Index_size_ ); ) {
      k = k+1;
      copy_add_list_( Permutations_, CurrentIndexes );

      //
      // Test Automorphism
      //
      int *ij_index = motif.getIndexes();
      
      // Apply permutation
      for( int i=0;  ij_index[i] != NOT_INDEX; i++) {
	ij_index[i] = CurrentIndexes[ ij_index[i] ];  
      }
      
      // Build new Matrix
      SparseAdjacency newMat(ij_index, NbRow, MotifSymetry );
      delete [] ij_index;

      // Test new matrix
      if ( newMat.operator==( motif ) ) {
	copy_add_list_( Automorphisms_, CurrentIndexes );
      }
    };

    //
    // Get Non Redondant permutations
    //
    bool test;
    bool found = false;
    int *test_perm = new int[Index_size_];
    list<int*>::iterator it = Permutations_->begin();

    // Skip the identity and remove it from NR Permutations.
    // Permutations_->erase( it );
    // it = Permutations_->begin();

    for ( ; it != Permutations_->end(); it++ ) {
      //
      // Loop over all atomorphisms
      //
      list<int*>::iterator ita = Automorphisms_->begin();
      // Skip the identity automorphism
      // ?? ita++;
      for ( ; ita != Automorphisms_->end(); ita++ ) {

	// Apply permutation

	for( int i=0;  i< Index_size_ ; i++) {
	  test_perm[i] = (*it)[ (*ita) [i] ];  
	}
      	
	// Remove test_perm in permutation list
	list<int*>::iterator itx = it;
	itx++;
	found = false;
	for ( ; itx != Permutations_->end() && !(found) ; itx++ ) {

	  test = true;
	  for( int i=0;  (i< Index_size_) && test ; i++) {
	    test = ( test_perm[i] == (*itx)[i] );  
	  }
	  if (test) {  
	    delete [] *itx;
	    itx = Permutations_->erase( itx );
	    found = true;
	  }
	}
      }
    }
    CurentIt_ = Permutations_->begin();
    delete [] test_perm;
  } else {
    LOGMSG( 1, std::cout, "not squared matrix", "" ); 
  }
  delete [] CurrentIndexes;
}

NonRedondantPermutation::~NonRedondantPermutation( ) {


  // Delete NR permutations
  for ( list<int*>::iterator it = Permutations_->begin() ; 
	it != Permutations_->end(); it++ ) {
   delete [] *it;
  }

  delete  Permutations_;

  // Delete Automorphisms
  for ( list<int*>::iterator it = Automorphisms_->begin() ; 
	it != Automorphisms_->end(); it++ ) {
    delete [] *it;
  }

  delete  Automorphisms_;
 
}

const int *
NonRedondantPermutation::iterate( const char * mode  ) {

  int *perm;
  int start = 0;
  
  string str(mode);
  if( str.compare("start") == 0 )  
    start = 1;;

  if (start) {
    CurentIt_ = Permutations_->begin();
    if( CurentIt_ != Permutations_->end() ) {
      perm = *CurentIt_;
      CurentIt_++;
    } else {
      perm = 0;
    }
  } else {
    if( CurentIt_ != Permutations_->end()) {
      perm = *CurentIt_;
      CurentIt_++;
    } else {
      perm = 0;
    }
  }

  return perm ;
}


