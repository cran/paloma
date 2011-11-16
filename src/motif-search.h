
# ifndef _SEARCH_MOTIF_H_
# define _SEARCH_MOTIF_H_

# include "sparse-adjacency.h"
# include "motif-tree.h"
# include "block-alloc.h"
# include <sys/time.h>
# include "shrink.h"
# include "util.h"

using namespace std;

// Used to display message Information
// # define MSG

# define TIMING

// Define output modes
const int CountMode=1,  WriteMode=2,  StoreMode=3, TCMode=4; 

// Define search modes
const int NonExactMode = 0; // or non-induiced
const int ShrinkMode = 1;   // Shrink mode, non-exact
const int ExactMode  = 2;   // Exact mode
const int DefaultSearchMode = NonExactMode; 

 
const int NullMode   = 10; // Don't search  

const int ComputeMotifs = 0; // Without shrink
const int ReadMotifs    = 1; // Read From a file 
const int WriteMotifs   = 2; // Read From a file 

// Chunk size for block allocation
const int ChunkSize = 2048;

class FindMotif {
  int64_t  _Count;
  int64_t  _CountAllMotifs;
  int            _OutputMode;
  ostream       *_OutputStream;
  int            _RealMotifSize;
  int           *_Constraints;
  bool          *_TopoClasses;
  int            _NumberOfChunks ;
  int           *_ChunkStorage; 
  vector< int*> *_FoundList;

  // Storage mode unused (up to now)
  BlockAllocation<int> *_FoundBList;

  // Additionnal counters
  int64_t    *_Counters;

  // Non-induced motif (non exact motif)
  // Undirected Case
  // k+1 step of the motif search.
  // Search the set all nodes G with have the same connectivity
  // that the k+1 node of the motif to detect. 
  void kpp_MatchUndirectedMotif(  
			       SparseAdjacency &motif,  
			       SparseAdjacency &G, 
			       int k_node, 
			       int *MapToG
			       );

  // Non-induced motif (non exact motif)
  // Directed Case
  // k+1 step of the motif search.
  // Search the set all nodes G with have the same connectivity
  // that the k+1 node of the motif to detect. 
  void kpp_MatchDirectedMotif(  
			     SparseAdjacency &motifSup,  
			     SparseAdjacency &motifInf,  
			     SparseAdjacency &G, 
			     int k_node, 
			     int *MapToG
			     );

  //
  //  Find colored motif
  //  Directed mode
  //  Use "fat" or the whole motif to map it to the colored Graph G
  //  Use the "Values" field, so that this field cannot be used for
  //  "shrinked" motifs.  
  void kpp_MatchUndirectedColoredMotif( 
				      SparseAdjacency &motif,
				      const int *MotifColor,
				      SparseAdjacency &G, 
				      const int *GColor, 
				      int k, 
				      int* MapToG 
				      );
  //
  //  Find colored motif
  //  Directed mode
  //  Use "fat" or the whole motif to map it to the colored Graph G
  //  Use the "Values" field, so that this field cannot be used for
  //  "shrinked" motifs. 
  void kpp_MatchDirectedColoredMotif( 
				      SparseAdjacency &motifSup, 
				      SparseAdjacency &motifInf,
				      const int *MotifColor,
				      SparseAdjacency &G, 
				      const int *GColor,
				      int k, 
				      int *MapToG
				    );

  //
  //  Find exact (induced) colored motif
  //  Directed mode
  //  Use "fat" or the whole motif to map it to the colored Graph G
  //  Use the "Values" field, so that this field cannot be used for
  //  "shrinked" motifs. 
  void kpp_MatchDirectedMotif( 
				      SparseAdjacency &motifSup, 
				      SparseAdjacency &motifInf,
				      SparseAdjacency &cmplSup, 
				      SparseAdjacency &cmplInf,
				      SparseAdjacency &G, 
				      int k, 
				      int *MapToG
				    );

  //
  //  Find exact (induced) colored motif
  //  Directed mode
  //  Use "fat" or the whole motif to map it to the colored Graph G
  //  Use the "Values" field, so that this field cannot be used for
  //  "shrinked" motifs. 
  void kpp_MatchDirectedColoredMotif( 
				      SparseAdjacency &motifSup, 
				      SparseAdjacency &motifInf,
				      SparseAdjacency &cmplSup, 
				      SparseAdjacency &cmplInf,
				      const int *MotifColor,
				      SparseAdjacency &G, 
				      const int *GColor,
				      int k, 
				      int *MapToG
				    );

  // Induced motif (Exact motif)
  // Undirected Case
  // k+1 step of the motif search.
  // Search the set all nodes G with have the same connectivity
  // that the k+1 node of the motif to detect. 
 void kpp_MatchUndirectedMotif( 
			      SparseAdjacency &motif,  
			      SparseAdjacency &cmpl,  
			      SparseAdjacency &G, 
			      int k, 
			      int* MapToG 
			       );

  // Induced motif (Exact motif)
  // Undirected Case
  // k+1 step of the motif search.
  // Search the set all nodes G with have the same connectivity
  // that the k+1 node of the motif to detect. 
 void kpp_MatchUndirectedColoredMotif( 
			      SparseAdjacency &motif,  
			      SparseAdjacency &cmpl,  
			      const int *MotifColor,
			      SparseAdjacency &G,
			      const int *GColor,
			      int k, 
			      int* MapToG 
			       );

 void kpp_AllUndirectedMotif( 
			    siblings_t *currents_sibs,
			    SparseAdjacency &G, 
			    int k, 
			    int id, 
			    int* MapToG,
			    int* constraints 
			     );

  //
  // Induced motif (Exact motif)
  // Undirected Case
  // FanMod algorithm
  // k+1 step of the motif search.
  // Search the set all nodes G with have the same connectivity
  // that the k+1 node of the motif to detect. 
 void kpp_MatchUndirectedMotifs( 
				int MotifSize,
				SparseAdjacency &G, 
				int k, 
				int* MapToG,
				int* StartNeighbours,
				int* EndNeighbours
				 );
  // Shrinked motif ??? version Michel
  void kpp_MatchUndirectedMotifShrink(  
			       SparseAdjacency &motif,  
			       SparseAdjacency &G, 
			       int k_node, 
			       int *MapToG,
			       Neibourgh *MultipleVerticesTwoSources,
			       Neibourgh *MultipleVerticesManySources,
			       int *ReverseOccurrence

			       );

  // Shrinked motif ???
  void kpp_MatchUndirectedMotifShrink(  
			       SparseAdjacency &motif,  
			       SparseAdjacency &G, 
			       int k_node, 
			       int *MapToG
			       );

// Shrinked motif ???
 void kpp_MatchDirectedMotifShrink(  
			     SparseAdjacency &motifSup,  
			     SparseAdjacency &motifInf,  
			     SparseAdjacency &G, 
			     int k_node, 
			     int *MapToG
			     );

  void _computeEquivalentNodes(  
			        SparseAdjacency &motif, 
			        int *Color
			      );
 public :

  FindMotif( ) { 
    _Count = 0;
    _CountAllMotifs = 0;
    _Counters = 0;
    _OutputMode   = CountMode;
# ifdef MSG
    _OutputStream = &cout;
# endif
    _RealMotifSize = 0;
    _Constraints   = 0;
    _TopoClasses   = 0;
    _FoundList     = 0;
    _NumberOfChunks = 0;
    _FoundBList     = 0;
    allocateFoundList();
  }

  FindMotif( int outputmode_, ostream *outputstream_  ) { 
    _Count = 0;
    _CountAllMotifs = 0;
    _Counters = 0;
    _OutputMode   = outputmode_;
    _OutputStream = outputstream_;
    _RealMotifSize = 0;
    _Constraints   = 0;
    _TopoClasses   = 0;
    _FoundList     = 0;
    allocateFoundList();
    _NumberOfChunks = 0;
    _FoundBList     = 0;
  }

  int64_t inline getCount() { return _Count; };
  const int64_t inline *getCounters() const { return _Counters; };
  int inline getOutputMode() { return _OutputMode; };
  void inline setOutputMode( int outputmode) { _OutputMode=outputmode; };

  void inline allocateCounters( int size ) { 
    if( _Counters != 0 ) 
      delete [] _Counters;
    _Counters = new int64_t[ size ];
    for( int64_t *it = _Counters, *it_end = _Counters + size ; 
	 it != it_end; it++) 
      *it = 0;
  }

  void inline deallocateCounters() { 
    if( _Counters ) {
      delete [] _Counters;
      _Counters = 0;
    }
  }
    
  void inline allocateConstraints( int size ) { 
    if( _Constraints != 0 ) 
      delete [] _Constraints;
    _Constraints = new int[ size ]; 
  }

  void inline deallocateConstraints() { 
    if( _Constraints ) {
      delete [] _Constraints;
      _Constraints = 0;
    }
  }
    
  const vector<int*> inline *getFoundList() { 
    return _FoundList ; 
  }

  void allocateFoundList() { 
    if( _FoundList ) 
      delete _FoundList;
    _FoundList = new vector<int*>;
    // cerr << _FoundList << endl;
  }

  void clearFoundList() { 
    if( _FoundList ) {
      _FoundList = deallocateFoundList(_FoundList); 
      delete _FoundList;
    }
    _FoundList = new vector<int*>;
    // cerr << _FoundList << endl;
  }  

  // Deallloc according to the list size and the ChunkSize
static  vector<int*> inline *deallocateFoundList( vector<int*> *foundlist) {
    int  NumberOfChunks =0 ;
    if( foundlist ) {
      NumberOfChunks = foundlist->size() / ChunkSize 
	+ ( (foundlist->size() % ChunkSize) != 0);
      for ( int i=0; i < NumberOfChunks ; i++) {
	int *ptr= (*foundlist)[i*ChunkSize];
	delete [] ptr;
      }
      delete foundlist;
      foundlist = 0;
    }
    return foundlist;
  }

void inline deallocateFoundList() { 
    if( _FoundList ) {
      _FoundList = deallocateFoundList(_FoundList );
    }
  }

  void inline allocateFoundBList( const int motif_size ) { 
    if( _FoundBList ) 
      delete _FoundBList;
    _FoundBList = new BlockAllocation<int>( motif_size, ChunkSize) ; 
  }

  void inline deallocateFoundBList() { 
    if( _FoundBList ) {
      delete _FoundBList;
      _FoundBList = 0;
    }
  }

  void inline allocateTopoClasses( int size ) { 
    if( _TopoClasses ) 
      delete [] _TopoClasses;
    _TopoClasses = new bool[ size * size ];
    for ( bool *it= _TopoClasses, 
	    *it_end = &_TopoClasses[ size * size ];
	  it != it_end; it++ ) {
      *it = false;
    } 
  }

  void inline deallocateTopoClasses() { 
    if( _TopoClasses ) {
      delete [] _TopoClasses;
      _TopoClasses = 0;
    }
  }

  vector<int*> *find(  SparseAdjacency &G, SparseAdjacency &motif,
		       bool use_constraints=true, bool color=false );
  vector<int*> *findExactly(  SparseAdjacency &G, SparseAdjacency &motif );

  vector<int*> *findShrink(  SparseAdjacency &G, SparseAdjacency &motif );

  // Find all motifs in a graph (generate all possible motifs)
  vector< int* > *findAll( SparseAdjacency &G, siblings_t &MotifTree,
			   int motif_size, bool motif_undirected );

  // Find all motifs in a graph (Fanmod type)
  vector< int* > *findAllMotifs(  SparseAdjacency &G, 
				  int MotifSize );


  int64_t writeOccurences(  
			     SparseAdjacency &G, 
			     SparseAdjacency &motifIn,  
			     SparseAdjacency &motifOut,  
			     int *MapToG
			   );

  vector< vector< int > > *getTopologicClasses( 
				   SparseAdjacency &motif, 
				   int *Color=0
				   );

  int computeNbrOfAutomorphism( 
			       SparseAdjacency &motif, 
			       int *Color=0
			      );

  vector< int > *getEquivalentNodes(  
				   SparseAdjacency &motif, 
				   int *Color
				   );

  // Test if motif1 is included non exactly (non-induced) in motif2
  bool isNonExactlyIncluded(  
		  SparseAdjacency &motif1, 
		  SparseAdjacency &motif2, 
		  int *Color1 = 0, 
		  int *Color2 = 0 
		  );

  // Test if motif1 is included exactly (induced) in motif2
  bool isExactlyIncluded(  
		  SparseAdjacency &motif1, 
		  SparseAdjacency &motif2, 
		  int *Color1 = 0, 
		  int *Color2 = 0 
		  );

  bool isNonExactIsomorph(  
		  SparseAdjacency &motif1, 
		  SparseAdjacency &motif2, 
		  int *Color1 = 0, 
		  int *Color2 = 0 
		  );

  ~FindMotif( ) {
    deallocateConstraints();
    deallocateTopoClasses();
    deallocateFoundList();
    /// ???    deallocateFoundList();
  } 

  inline void setConstraints( int *Constraints ) {
    if ( _Constraints ) delete [] _Constraints;
    _Constraints = Constraints;
  }

  inline void setConstraints( 
			     const int size, 
			     int value=NOT_INDEX 
			    ) {
    allocateConstraints( size );
    for( int *it=_Constraints, *it_end = &_Constraints[ size ];
	 it != it_end; it++) {
      *it = value;
    }
  }
  void writeList( const char *output_name, 
		  const vector<int*> *found, 
		  const int motif_size, 
		  bool block_alloc = true
		);
};
  
// High level method
// Change internal state of the object
void FindAllMotif( SparseAdjacency &G, int k, FindMotif &find);

void  buildAndfindAllMotifTree( SparseAdjacency &G, int k, FindMotif &find);
void  buildAndfindAllMotif( SparseAdjacency &G, int k, FindMotif &find, 
			    vector<int *> &motifs, vector<int64_t> &counts,
			    int search_mode = 0, int generate_mode = 0);

void  buildAndfindAllMotifShrink( SparseAdjacency &G, int k, FindMotif &find, 
			    vector<int *> &motifs, vector<int64_t> &counts);

SparseAdjacency ShrinkMotif( SparseAdjacency & motif, int* MapShrinkToOriginal);



// Motif database access
void setMotifDirectory( const char * motif_dir);

# endif

