# include <list>
# include <algorithm>
# include "sparse-adjacency.h"

class NonRedondantPermutation {

 protected:
  // Not used 
  int PermutationMode_;

  list<int*>::iterator CurentIt_;
  int Index_size_;

  list<int*> *Permutations_;
  // Does not contain the identity
  list<int*> *Automorphisms_;

  
  void copy_add_list_( list<int*> *perm , int *indexes ) {
    int *copy = new int[ Index_size_ ];
    memcpy( copy, indexes, sizeof( int ) * Index_size_ );
    perm->push_back( copy );

  }

  void copy_list_(list<int*> *SPerm , list<int*> *DPerm) {

    for ( list<int*>::iterator it = SPerm->begin() ; 
	  it != SPerm->end(); it++ ) {
      copy_add_list_( DPerm, *it);
    }
  }

 public:

  //! \brief Constructor from a motif.
  //! Generate all automorphisms and non redondant permutations
  //! of the specified motif. 
  //! \param motif .
  NonRedondantPermutation( const SparseAdjacency & motif );

  //! \brief Copy constructor
  //! \param A permutation to copy.
  NonRedondantPermutation( const NonRedondantPermutation &A) {

    PermutationMode_ = A.PermutationMode_;

    CurentIt_   = A.CurentIt_;
    Index_size_ = A.Index_size_;

    Permutations_  = new list<int*>;
    Automorphisms_ = new list<int*>;
    copy_list_( A.Permutations_,  Permutations_  );
    copy_list_( A.Automorphisms_, Automorphisms_ );
  }

  //! \brief Destructor
  ~NonRedondantPermutation();

  //! \brief Iterate permutation.
  //! \param mode :
  //!        "next" : iterate over the whole permutation set.
  //!        "restart" : restart the permutation iteration.
  //! \return The curent indexes are returned. 
  const int *iterate( const char * mode = "next");

  int64_t getCardNRPermutation() { return  Permutations_->size(); };
  int64_t getCardAutomorphism() { return  Automorphisms_->size(); };
  int64_t getCardPermutation() { return  
      (Automorphisms_->size() * Permutations_->size()); };
  void permuteMatrix( const int *);
};

