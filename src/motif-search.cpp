# include <list>
# include "motif-search.h"
# include "heap.h"
# include "combinatorics.h"
# include "block-alloc.h"
# include "heap-union.h"

# ifndef _MOTIF_SEARCH_H_
# define _MOTIF_SEARCH_H_

// Motif file name
const char *MotifDir = ".";
const char *PatternMotifFileName = "%s/motif-%c-%d.db";
static char MotifFileName[1024];
static ofstream fmout;
static ifstream fmin;

static int64_t _CountAllMotifs;
static int64_t _NbrOfMotifs;

void getChildMotifs( vector< SparseAdjacency *> &AncestorMotifs, 
	       vector< vector< vector<int> > *> &AncestorCT, 
	       vector< SparseAdjacency *> &ChildMotifs,
		vector< vector< vector<int> > *> &ChildCT, 
		int motif_size, 
		int nb_edges,
		bool symmetry_, 
		FindMotif &find);

// Use it own FindMotif object
// Use high level method (use its find object)
int *getUniqPath(  SparseAdjacency & motif_, const int *Color_=0 ) {


  SparseAdjacency motif( motif_ );
  int MotifSize =  motif.getNbrRow();
  int *Constraints = new int[ MotifSize ];
  int *Color = new int[ MotifSize ];

  vector <int> *TCToNodes;

  // Init Constraints
  for( int i=0; i < MotifSize; i++) {
    Constraints[i] = NOT_INDEX;
  }

  // Init Color
  if( Color_ ) { 
    for( int i=0; i < MotifSize; i++) {
      Color[i] = Color_[i];
    }
  } else {
    for( int i=0; i < MotifSize; i++) {
      Color[i] = 0;
    }
  }
  int max_color = 0;

  FindMotif search( TCMode, 0);

  for( int i=0; i < MotifSize; i++) {

    // TCToNodes nodes in the same class
    TCToNodes = search.getEquivalentNodes( motif, Color);

    // Build TCToNodes mapping
    // TODO put inside nr-...h

    // Topological classes must be split to singleton
    //
    // Previous classes have been reduced to singleton
    if ( TCToNodes[ i ].size() > 1 ) {
      for( vector<int>::iterator it = TCToNodes[i].begin();  
	   it < TCToNodes[i].end(); it++) {
	if ( *it == i ) {
	  // Change the Color of considered node i
	  Color[*it] = ++max_color; 
	} else {
	  // Set a contraint ( Constraints[ i ] > k ) in relation 
	  // to the considered node i
	  Constraints[ *it ] = i;
	}
      }
    }

    delete [] TCToNodes;
  }
  delete [] Color;
  return Constraints;
}    



// High level method
// Change internal state of the object
vector< int* > *FindMotif::find(  
			   SparseAdjacency &G, 
			   SparseAdjacency &motif,
			   bool use_constraints,
			   bool color
			   ) {

  // TODO Choix entre les deux listes
  allocateFoundList();
  // ???
  allocateFoundBList( motif.getNbrRow() );
  _Count = 0;

  int G_NbrNodes = G.getNbrRow();

  // Remove poor neighbouring
  _RealMotifSize = motif.getNbrRow();


  // Motif
  int motif_NbrNodes = motif.getNbrRow();
  int *MapToG = new int[ motif_NbrNodes];
  int k_node = 0;

  // Change of type int64_t to int for color
  // TODO : Change in the SparseAdjacency Matrix the Value 
  //        (int64_t to a template) 
  int *MColor = 0 ;
  int *GColor = 0;
  const int64_t *p_motif = motif.getNodeValues(); 
  const int64_t *p_G = G.getNodeValues(); 
  if( color ) {
    MColor = new int [motif_NbrNodes];
    GColor = new int [G_NbrNodes];

    for ( int i=0; i < motif_NbrNodes; i++) {
      MColor[i] = *p_motif++;
    }
    for ( int i=0; i < G_NbrNodes; i++) {
      GColor[i] = *p_G++;
    }
  }

  // Get motif constraints
  if ( use_constraints ) {
    if ( color )
      setConstraints( getUniqPath( motif , MColor) );
    else
      setConstraints( getUniqPath( motif ));
  } else {
     // No constraints
    setConstraints( motif_NbrNodes );
  }
# ifdef MSG 
  printList( "Motif constraint ", _Constraints, 
	     &_Constraints[ motif.getNbrRow() ], true);
# endif

# ifdef MSG 
  motif.print("motif"); 
  G.print("Network"); 
# endif

  if ( motif.getSymmetry() &&  G.getSymmetry() ) {

    for ( int i = 0; i < G_NbrNodes; i++ ) {

      MapToG[0] = i;
      if( G.getAllRowSize( i ) != 0) {
 
	if( color ) {
	  if (MColor[0] == GColor[i]) 
	    kpp_MatchUndirectedColoredMotif( motif, MColor, G, GColor, 
					     k_node, MapToG
					     );
	} else {
	  kpp_MatchUndirectedMotif( motif, G, k_node, 
				    MapToG
				    );
	}
      }
    }
    
 
  } else if ( ! (motif.getSymmetry() ||  G.getSymmetry() ) ) {
    SparseAdjacency MotifUpper( motif.getTriangular( true ) );
    SparseAdjacency MotifLower( motif.getTriangular( false ) );

# ifdef MSG 
    MotifUpper.print( "Upper");
    MotifLower.print( "Lower");
# endif    

    for ( int i = 0; i < G_NbrNodes; i++ ) {
      MapToG[0] = i;
      if( color ) {
	if (MColor[0] == GColor[i]) 
	  kpp_MatchDirectedColoredMotif( MotifUpper,  MotifLower, MColor, 
					 G, GColor, 
					 k_node, MapToG
					 );
      } else {
	kpp_MatchDirectedMotif(  MotifUpper,  MotifLower, G, k_node, 
				 MapToG
				 );
      }
    }
  } else {
#   ifdef MSG 
    cerr << "The Graph and the motif not the same characteristics "
	 << "(directed/undirected)" << endl;
#   endif
  }
   vector< int*> *found_list=_FoundList;
  _FoundList = 0;

  if( color ) {
    delete [] MColor ;
    delete [] GColor ;
  }

  delete [] MapToG;

  return( found_list ); 
}
// High level method
// Change internal state of the object
// Find exact motif (induced)
vector< int* > *FindMotif::findExactly(  
			   SparseAdjacency &G, 
			   SparseAdjacency &motif 
			   ) {

  // TODO Choix entre les deux listes
  allocateFoundList();
  // ???
  allocateFoundBList( motif.getNbrRow() );
  _Count = 0;

  int G_NbrNodes = G.getNbrRow();

  // Remove poor neighbouring
  _RealMotifSize = motif.getNbrRow();


  // Motif
  int motif_NbrNodes = motif.getNbrRow();
  int *MapToG = new int[ motif_NbrNodes];
  int k_node = 0;

  // Get motif constraints
  setConstraints( getUniqPath( motif ));

# ifdef MSG 
  printList( "Motif constraint ", _Constraints, 
	     &_Constraints[ motif.getNbrRow() ], true);
# endif
  printList( "Motif constraint ", _Constraints, 
	     &_Constraints[ motif.getNbrRow() ], true);


# ifdef MSG 
  motif.print("motif"); 
  G.print("Network"); 
# endif

  int *ij_index_cmpl =  motif.getIndexes( true, false );

  if ( motif.getSymmetry() &&  G.getSymmetry() ) {
    // For exact motif
    SparseAdjacency cmpl( ij_index_cmpl, 
			  motif.getNbrRow(),   
			  motif.getSymmetry() );
    // cmpl.print("cmpl");
    
    for ( int i = 0; i < G_NbrNodes; i++ ) {

      MapToG[0] = i;
      if( G.getAllRowSize( i ) != 0)
	// kpp_MatchUndirectedMotif( motif, G, k_node, 
	//		  MapToG
	//		);

        // Exact
	kpp_MatchUndirectedMotif( motif, cmpl, G, k_node, 
        			  MapToG
       			  );
    }
    
 
  } else if ( ! (motif.getSymmetry() ||  G.getSymmetry() ) ) {
    SparseAdjacency MotifUpper( motif.getTriangular( true ) );
    SparseAdjacency MotifLower( motif.getTriangular( false ) );

# ifdef MSG 
    MotifUpper.print( "Upper");
    MotifLower.print( "Lower");
# endif    
    // For exact motif

    SparseAdjacency cmpl( ij_index_cmpl, 
			  motif.getNbrRow(),   
			  motif.getSymmetry() );

    SparseAdjacency cmplUpper( cmpl.getTriangular( true ) );
    SparseAdjacency cmplLower( cmpl.getTriangular( false ) );

    for ( int i = 0; i < G_NbrNodes; i++ ) {
      MapToG[0] = i;
      kpp_MatchDirectedMotif(  MotifUpper,  MotifLower,   
			       cmplUpper, cmplLower, G, k_node, 
    			       MapToG
			    );
    }
  } else {
#   ifdef MSG
    cerr << "The Graph and the motif not the same characteristics "
	 << "(directed/undirected)" << endl;
#   endif
  }
   vector< int*> *found_list=_FoundList;

  delete [] MapToG;
  _FoundList = 0;

  return( found_list ); 
}

// Low Level routine
// method
// Use the class internals
void FindMotif::_computeEquivalentNodes(  
			   SparseAdjacency &motif, 
			   int *Color
			   ) {

  _Count = 0;
  _OutputMode = TCMode;

  // Motif
  int motif_size = motif.getNbrRow();
  int *MapToG = new int[ motif_size];
  int k_node = 0;


  int *color_tmp_ = new int[ motif_size ];
  int *Color_;
  if (Color == 0) {
    for ( int *it= color_tmp_, 
	    *it_end = &color_tmp_[ motif_size ];
	  it != it_end; it++ ) {
      *it = 0;
    }
    Color_ = color_tmp_;
  } else {
    Color_ = Color;
  }

  // Set no constraint 
  setConstraints( motif_size );

  // Allocates and initializes array for TopologicalClasses
  allocateTopoClasses( motif_size );

  if ( motif.getSymmetry()) {
    for ( int i = 0; i < motif_size; i++ ) {
      MapToG[0] = i;
      if( motif.getAllRowSize( i ) != 0) 
	kpp_MatchUndirectedColoredMotif( motif, Color_, motif, Color_, 
					 k_node, MapToG
				       );
    }
    
 
  } else {
    SparseAdjacency motifUpper( motif.getTriangular( true ) );
    SparseAdjacency motifLower( motif.getTriangular( false ) );

# ifdef MSG 
    motifUpper.print( "Upper");
    motifLower.print( "Lower");
# endif    

    for ( int i = 0; i < motif_size; i++ ) {
      MapToG[0] = i;
      kpp_MatchDirectedColoredMotif(  motifUpper,  motifLower, Color_,
				      motif, Color_, 
				      k_node, MapToG
				   );
    }
  }
  delete [] MapToG;
  delete [] color_tmp_;
}

// Low Level routine
// method
// Use the class internals
// Return a array of topologic classes with contains node ID
vector< vector< int > > *FindMotif::getTopologicClasses(  
			   SparseAdjacency &motif, 
			   int *Color
			   ) {

  _computeEquivalentNodes( motif, Color );

  int motif_size = motif.getNbrRow();

  vector < vector< int > > *TCToNodes = new vector< vector< int > > ;
  TCToNodes -> resize( motif_size );
  bool *not_done = new bool[motif_size];
  for( int i = 0 ;  (i < motif_size); i++) {
    not_done[i] = true;
  }

  int CurClass = 0;
  for( int i = 0 ;  (i < motif_size) ; i++) {
    if ( not_done[i] ) {
      bool *it = _TopoClasses + ( motif_size * i);
      // Nodes are in the same classe
      for( int j = 0 ;  (j < motif_size) ; j++) {
	if( it[j] ) {
	  (*TCToNodes)[ CurClass ].push_back( j );
	  not_done[j] = false;
	}  
      }
      CurClass++;
    }
  }

  for( int i = motif_size-1 ;  i > -1 ; i--) {
    if(  (*TCToNodes)[i].size() == 0 ) {
      TCToNodes -> pop_back();
    }
  }
  delete [] not_done;
  return( TCToNodes ); 
}

int computeNbrOfAutomorphism(  SparseAdjacency & motif_ ) {

  SparseAdjacency motif( motif_ );
  // int *Color = new int[ MotifSize ];

  FindMotif searchAutoMorph( CountMode, 0);

  searchAutoMorph.find( motif_, motif_, false );
  return( searchAutoMorph.getCount() );
}

// Low Level routine
// method
// Use the class internals
vector< int > *FindMotif::getEquivalentNodes(  
			   SparseAdjacency &motif, 
			   int *Color
			   ) {

  _computeEquivalentNodes( motif, Color);

  int motif_size = motif.getNbrRow();

  vector< int > *TCToNodes = new vector< int > [motif_size];
 
  bool *not_done = new bool[motif_size];
  for( int i = 0 ;  (i < motif_size); i++) {
    not_done[i] = true;
  }
  bool again = true;

  int CurClass = 0;
  for( int i = 0 ;  (i < motif_size) && again ; i++) {
    bool *it = _TopoClasses + ( motif_size * i);
    again = true;
    bool NewClass;
    // Nodes are in the same classe
    for( int j = 0 ;  (j < motif_size) && again ; j++) {
      if( it[j] && not_done[j] ) {
	TCToNodes[ CurClass ].push_back( j );
	not_done[j] = false;
	NewClass = true;
      }  
      again = again || not_done[j];
    }
    if( NewClass)
      CurClass++;
  }

  delete [] not_done;
  return( TCToNodes ); 
}


// Low level method
// Set internal state of the FindMotif Class 
bool FindMotif::isNonExactIsomorph(  
			   SparseAdjacency &motif1, 
			   SparseAdjacency &motif2, 
			   int *Color1, 
			   int *Color2  
			   ) {

 // Get Topological Classes
  bool same_dim = ( motif1.getNbrRow() == motif2.getNbrRow() )
    && ( motif1.getSymmetry() == motif2.getSymmetry() )
    && ( motif1.getNbrOfValues() == motif2.getNbrOfValues() ) ;
  
  bool isomorph = false;
  if (same_dim ) {
    isomorph = isNonExactlyIncluded( motif1, motif2, Color1, Color2);
  }
  return isomorph;
}

// Low level method
// Set internal state of the FindMotif Class 
bool FindMotif::isNonExactlyIncluded(  
			   SparseAdjacency &motif1, 
			   SparseAdjacency &motif2, 
			   int *Color1, 
			   int *Color2  
			   ) {

  _Count = 0;
  _OutputMode = CountMode;

  int motif_size = motif1.getNbrRow();
  int graph_size = motif2.getNbrRow();
  int *MapToG = new int[ motif_size];


    if ( Color1 && Color2) {

      //
      // Colored graph
      //
      
      // Motif
       int k_node = 0;

      setConstraints( getUniqPath( motif1, Color1 ) );
      printList( "Motif constraint ", _Constraints, 
		 &_Constraints[ motif1.getNbrRow() ], true);
      
      if ( motif1.getSymmetry()) {
	for ( int i = 0; i < graph_size; i++ ) {
	  MapToG[0] = i;
	  if( motif2.getAllRowSize( i ) != 0 && (Color1[0] == Color2[i]) ) 
	    kpp_MatchUndirectedColoredMotif( motif1, Color1, motif2, Color2, 
					     k_node, MapToG
					     );
	}
	
      } else {
	SparseAdjacency motifUpper( motif1.getTriangular( true ) );
	SparseAdjacency motifLower( motif1.getTriangular( false ) );
	
# ifdef MSG 
	motifUpper.print( "Upper");
	motifLower.print( "Lower");
# endif    
      
	for ( int i = 0; i < graph_size; i++ ) {
	  MapToG[0] = i;
	  if( Color1[0] == Color2[i] )
	    kpp_MatchDirectedColoredMotif(  motifUpper,  motifLower, Color1,
					    motif2, Color2, 
					    k_node, MapToG
					    );
	}
      }
    } else {
      
      //
      //  Uncolored graph
      // 
      
      // Remove poor neighbouring
      _RealMotifSize = motif1.getNbrRow();

      // Motif
      // unused : int shrink_size = motif1.getNbrRow();
      int k_node = 0;
      
      // Get motif constraints
      setConstraints( getUniqPath( motif1 ) );
      
#   ifdef MSG 
      printList( "Motif constraint ", _Constraints, 
		 &_Constraints[ motif1.getNbrRow() ], true);
#   endif

      if ( motif1.getSymmetry() ) {
	
	for ( int i = 0; i < graph_size; i++ ) {
	  MapToG[0] = i;
	  if( motif2.getAllRowSize( i ) != 0) 
	    kpp_MatchUndirectedMotif( motif1, motif2, k_node, 
				      MapToG
				      ) ;
	}
	
      } else {
	SparseAdjacency MotifUpper( motif1.getTriangular( true ) );
	SparseAdjacency MotifLower( motif1.getTriangular( false ) );
	
#       ifdef MSG 
	MotifUpper.print( "Upper");
	MotifLower.print( "Lower");
#       endif    

	for ( int i = 0; i < motif_size; i++ ) {
	  MapToG[0] = i;
	  kpp_MatchDirectedMotif( MotifUpper,  MotifLower, motif2, k_node, 
				  MapToG
				  );
	}
      } 
    }
    delete [] MapToG ;

  return  (_Count != 0) ;
}

// Low level method
// Set internal state of the FindMotif Class 
bool FindMotif::isExactlyIncluded(  
			   SparseAdjacency &motif1, 
			   SparseAdjacency &motif2, 
			   int *Color1, 
			   int *Color2  
			   ) {

  _Count = 0;
  _OutputMode = CountMode;

  int motif_size = motif1.getNbrRow();
  int graph_size = motif2.getNbrRow();
  int *MapToG = new int[ motif_size];

    if ( Color1 && Color2) {

      //
      // Colored graph
      //
      
      // Motif
      int k_node = 0;

      setConstraints( getUniqPath( motif1, Color1 ) );
      // printList( "Motif constraint ", _Constraints, 
      //	 &_Constraints[ motif1.getNbrRow() ], true);

      int *ij_index_cmpl =  motif1.getIndexes( true, false );

      if ( motif1.getSymmetry()) {
	SparseAdjacency cmpl1( ij_index_cmpl, 
			      motif1.getNbrRow(),   
			      motif1.getSymmetry() );
	for ( int i = 0; i < graph_size; i++ ) {
	  MapToG[0] = i;
	  if( motif2.getAllRowSize( i ) != 0 && ( Color1[0] == Color2[i] ) ) 
	    kpp_MatchUndirectedColoredMotif( motif1, cmpl1, Color1, 
					     motif2, Color2, 
					     k_node, MapToG
					     );

	}
	
      } else {
	SparseAdjacency motifUpper( motif1.getTriangular( true ) );
	SparseAdjacency motifLower( motif1.getTriangular( false ) );

	SparseAdjacency cmpl( ij_index_cmpl, 
			      motif1.getNbrRow(),   
			      motif1.getSymmetry() );

	SparseAdjacency cmplUpper( cmpl.getTriangular( true ) );
	SparseAdjacency cmplLower( cmpl.getTriangular( false ) );
	
# ifdef MSG 
	motifUpper.print( "Upper");
	motifLower.print( "Lower");
# endif    

	for ( int i = 0; i < graph_size; i++ ) {
	  MapToG[0] = i;
	  if( Color1[0] == Color2[i] )
	    kpp_MatchDirectedColoredMotif(  motifUpper,  motifLower, 
					    cmplUpper, cmplLower,
					    Color1,
					    motif2, Color2, 
					    k_node, MapToG
					    );
	}
      }
      delete [] ij_index_cmpl;
    } else {
      
      //
      //  Uncolored graph
      // 
      
      // Remove poor neighbouring
      _RealMotifSize = motif1.getNbrRow();

      // Motif
      // unused : int shrink_size = motif1.getNbrRow();
      int k_node = 0;
      
      // Get motif constraints
      setConstraints( getUniqPath( motif1 ) );
      
#   ifdef MSG 
      printList( "Motif constraint ", _Constraints, 
		 &_Constraints[ motif1.getNbrRow() ], true);
#   endif

      int *ij_index_cmpl =  motif1.getIndexes( true, false );

      if ( motif1.getSymmetry() ) {
	SparseAdjacency cmpl1( ij_index_cmpl, 
			      motif1.getNbrRow(),   
			      motif1.getSymmetry() );
	
	for ( int i = 0; i < graph_size; i++ ) {
	  MapToG[0] = i;
	  if( motif2.getAllRowSize( i ) != 0) 
	    kpp_MatchUndirectedMotif( motif1, cmpl1, motif2, k_node, 
				      MapToG
				      ) ;
	}
	
      } else {
	SparseAdjacency MotifUpper( motif1.getTriangular( true ) );
	SparseAdjacency MotifLower( motif1.getTriangular( false ) );
	
#       ifdef MSG 
	MotifUpper.print( "Upper");
	MotifLower.print( "Lower");
#       endif    

	SparseAdjacency cmpl( ij_index_cmpl, 
			      motif1.getNbrRow(),   
			      motif1.getSymmetry() );

	SparseAdjacency cmplUpper( cmpl.getTriangular( true ) );
	SparseAdjacency cmplLower( cmpl.getTriangular( false ) );

	for ( int i = 0; i < motif_size; i++ ) {
	  MapToG[0] = i;
	  kpp_MatchDirectedMotif( MotifUpper,  MotifLower, 
				  cmplUpper, cmplLower, motif2, k_node, 
				  MapToG
				  );
	}
      }
      delete ij_index_cmpl;
    }

    delete [] MapToG ;

    return  (_Count != 0) ;
}

void FindMotif::kpp_MatchUndirectedMotif( 
			       SparseAdjacency &motif,  
			       SparseAdjacency &G, 
			       int k, 
			       int* MapToG 
			     ) {
  //
  // Undirected case
  //
    
  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif


  if ( (k) == motif.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

#   ifdef MSG
    cerr << "  New motif is stored " << endl;
#   endif

    if( _OutputMode == CountMode ) {
      _Count++;
      // cerr << "count  =" << _Count << endl;
      // printList( " MapToG : ", MapToG,  &MapToG[k+1], true);
      
    } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
      //_Count += writeOccurences( G, motif, motif, MapToG );   

      // TODO storage
      _FoundBList -> storeData( MapToG );
      _Count++;
    }

  } else {


    /////////////////////////////////////////////////////////
    //
    // Get neighbours u of k+1 in the motif
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    const int *V_motif_kpp      = motif.getCol(k+1);
    const int *V_motif_kpp_end  = motif.getColEnd(k+1);
    int V_motif_kpp_size = V_motif_kpp_end - V_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Neighbours  {u} = V_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", V_motif_kpp, V_motif_kpp_end, true );
#   endif

    //
    // Make the intersection :
    //    INTER{ V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ] ) }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
    //

    //
    // Intersection initialisation with V_G( MapToG (u=0) )
    //
    const int *it = V_motif_kpp;


    // Set to intersect
    const int **L = new const int*[V_motif_kpp_size];
    int       *TL = new int[V_motif_kpp_size];

    //
    // Store V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ]) set 
    // sorted heap.
    //
    int i = 0;
    for ( it = V_motif_kpp; (it != V_motif_kpp_end); it++) {
	L[i]    = G.getAllRow( MapToG[*it] );
	TL[i++] = G.getAllRowSize( MapToG[*it]);
    }
    Heap H(L, TL, V_motif_kpp_size );

#   ifdef MSG
    H.Debug();
#   endif

    //
    // Get the intersection list.
    //
    int a = H.NextCommonElement();
    while (a != NOT_INDEX) {

      // Apply constraints 
      if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	  || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {
	// if ( MapToG[ k+1 ] > MapToG[ constraint[ k+1 ] ]   ]
	if (!IsIn( MapToG, k+1, a)) {  

	  // New sub-motif found
	  MapToG[k+1] = a;
	  kpp_MatchUndirectedMotif( motif, G, k+1, MapToG );
	}
      }
      a = H.NextCommonElement();
    }
    delete[] L;
    delete[] TL;

  }
}


bool IsInNeighboring( SparseAdjacency &G, int *nodes_beg, int *nodes_end, int w) {
  int k = nodes_end - nodes_beg;
  const int **begin       = new const int*[k+1];
  const int **end    = new const int*[k+1];
  
  begin[0] = nodes_beg;
  end[0]   = nodes_end;
  int m = 1;
  
  for ( int *it = nodes_beg; it != nodes_end; it++ ) {
    begin[m]  = G.getAllRow( *it ) ;
    end[m]    = G.getAllRowEnd( *it );
    if( begin[m] != end[m] ) {
      m++;
    }
  }
  
  HeapUnion H( begin, end, m);
  // H.status();

  bool not_found = true;
  for( int a = H.nextElement();
       a != NOT_INDEX && not_found; a = H.nextElement() ) {
    //  H.status();
    if( a == w ) {
      not_found = false;
    }     
  } 
  
  delete [] begin;
  delete [] end;

  return ! not_found; 
}

// NbrNeighbours must be deallocated by the caller
int *getNumberOfNeighboursInSubSet( SparseAdjacency &G, 
				    int *nodes_beg, int *nodes_end) {
  int k = nodes_end - nodes_beg ;
  int *NbrNeighbours = new int[ k ];
  const int **begin  = new const int*[2];
  int *size  = new int[2];
  
  begin[0] = nodes_beg;
  size[0]  = nodes_end - nodes_beg;
  int i_k = 0;

  for ( int *it = nodes_beg; it != nodes_end; it++, i_k++ ) {
    begin[1]  = G.getAllRow( *it ) ;
    size[1]    = G.getAllRowSize( *it );

    if( size[1] != 0 ) {
      Heap H( begin, size, 2 );
      int count = 0;
      for( int a = H.NextCommonElement(); a != NOT_INDEX; 
	   a = H.NextCommonElement()) {
	count++;
      }
      NbrNeighbours[i_k] = count;

    } else {
      NbrNeighbours[i_k] = 0;
    }
  }
  
  delete [] begin;
  delete [] size;

  return NbrNeighbours; 
}

  
void FindMotif::kpp_MatchUndirectedMotif( 
			       SparseAdjacency &motif,  
			       SparseAdjacency &cmpl,  

			       SparseAdjacency &G, 
			       int k, 
			       int* MapToG 
			     ) {
  //
  // Undirected case
  //
    
  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif


  if ( (k) == motif.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

#   ifdef MSG
    cerr << "  New motif is stored " << endl;
#   endif

    if( _OutputMode == CountMode ) {
      _Count++;
    } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
      //      _Count += writeOccurences( G, motif, motif, MapToG );   

      // TODO storage
      _Count++;

    }

  } else {


    /////////////////////////////////////////////////////////
    //
    // Get neighbours u of k+1 in the motif
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    const int *V_motif_kpp      = motif.getCol(k+1);
    const int *V_motif_kpp_end  = motif.getColEnd(k+1);
    int V_motif_kpp_size = V_motif_kpp_end - V_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Neighbours  {u} = V_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", V_motif_kpp, V_motif_kpp_end, true );
#   endif

    //
    // Make the intersection :
    //    INTER{ V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ] ) }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
    //

    //
    // Intersection initialisation with V_G( MapToG (u=0) )
    //
    const int *it = V_motif_kpp;

    // Set to intersect
    const int **L = new const int*[V_motif_kpp_size];
    int       *TL = new int[V_motif_kpp_size];

    //
    // Store V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ]) set 
    // sorted heap.
    //
    int i = 0;
    for ( it = V_motif_kpp; (it != V_motif_kpp_end); it++) {
      L[i]    = G.getAllRow( MapToG[*it] );
      TL[i++] = G.getAllRowSize( MapToG[*it]);
    }
    Heap H(L, TL, V_motif_kpp_size );

#   ifdef MSG
    H.Debug();
#   endif

    //
    // Get the intersection list.
    //

    const int *V_cmpl_kpp      = cmpl.getCol(k+1);
    const int *V_cmpl_kpp_end  = cmpl.getColEnd(k+1);
    // int V_cmpl_kpp_size        = V_cmpl_kpp_end - V_cmpl_kpp;

    int a = H.NextCommonElement();
    while (a != NOT_INDEX) {

      // Apply constraints 
      if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	  || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {
	// if ( MapToG[ k+1 ] > MapToG[ constraint[ k+1 ] ]   ]
	if (!IsIn( MapToG, k+1, a)) {  

	  // beg
	  bool not_found = true;
	  for ( it = V_cmpl_kpp; (it != V_cmpl_kpp_end) && (not_found); it++) {
	    // G.getAllRow( MapToG[*it] );
	    not_found = ! G.isInAllRow( MapToG[*it], a );
	  }

	  if( not_found ) {
	    // end 
	    // New sub-motif found
	    MapToG[k+1] = a;
	    kpp_MatchUndirectedMotif( motif, cmpl, G, k+1, MapToG );
	  }
	}
      }
      a = H.NextCommonElement();
    }
    delete[] L;
    delete[] TL;

  }
}

void FindMotif::kpp_MatchUndirectedColoredMotif( 
			       SparseAdjacency &motif,  
			       SparseAdjacency &cmpl,  
			       const int *MotifColor,
			       SparseAdjacency &G, 
			       const int *GColor,
			       int k, 
			       int* MapToG 
			     ) {
  //
  // Undirected case
  //
    
  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif


  if ( (k) == motif.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

#   ifdef MSG
    cerr << "  New motif is stored " << endl;
#   endif

    if( _OutputMode == CountMode ) {
      _Count++;
    } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
      //      _Count += writeOccurences( G, motif, motif, MapToG );   

      // TODO storage
      _Count++;

    }

  } else {


    /////////////////////////////////////////////////////////
    //
    // Get neighbours u of k+1 in the motif
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    const int *V_motif_kpp      = motif.getCol(k+1);
    const int *V_motif_kpp_end  = motif.getColEnd(k+1);
    int V_motif_kpp_size = V_motif_kpp_end - V_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Neighbours  {u} = V_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", V_motif_kpp, V_motif_kpp_end, true );
#   endif

    //
    // Make the intersection :
    //    INTER{ V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ] ) }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
    //

    //
    // Intersection initialisation with V_G( MapToG (u=0) )
    //
    const int *it = V_motif_kpp;


    // Set to intersect
    const int **L = new const int*[V_motif_kpp_size];
    int       *TL = new int[V_motif_kpp_size];

    //
    // Store V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ]) set 
    // sorted heap.
    //
    int i = 0;
    for ( it = V_motif_kpp; (it != V_motif_kpp_end); it++) {
      L[i]    = G.getAllRow( MapToG[*it] );
      TL[i++] = G.getAllRowSize( MapToG[*it]);
    }
    Heap H(L, TL, V_motif_kpp_size );

#   ifdef MSG
    H.Debug();
#   endif

    //
    // Get the intersection list.
    //

    const int *V_cmpl_kpp      = cmpl.getCol(k+1);
    const int *V_cmpl_kpp_end  = cmpl.getColEnd(k+1);
    // int V_cmpl_kpp_size        = V_cmpl_kpp_end - V_cmpl_kpp;

    int a = H.NextCommonElement();
    while (a != NOT_INDEX) {

      // Apply constraints 
      if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	  || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {
	// if ( MapToG[ k+1 ] > MapToG[ constraint[ k+1 ] ]   ]
	if (!IsIn( MapToG, k+1, a)) {  

	  // Same Color
	  if ( GColor[ a ] == MotifColor[k+1] ) {
	    // beg
	    bool not_found = true;
	    for ( it = V_cmpl_kpp; (it != V_cmpl_kpp_end) && (not_found); it++) {
	      // G.getAllRow( MapToG[*it] );
	      not_found = ! G.isInAllRow( MapToG[*it], a );
	    }

	    if( not_found ) {
	      // end 
	      // New sub-motif found
	      MapToG[k+1] = a;
	      kpp_MatchUndirectedColoredMotif( motif, cmpl, MotifColor, 
					       G, GColor, k+1, MapToG );
	    }
	  }
	}
      }
      a = H.NextCommonElement();
    }
    delete[] L;
    delete[] TL;

  }
}


// Non-induced motif (non exact motif)
// k+1 step of the motif search.
// Search the set all nodes G with have the same connectivity
// that the k+1 node of the motif to detect. 
void FindMotif::kpp_MatchDirectedMotif( 
			     SparseAdjacency &motifSup, 
			     SparseAdjacency &motifInf, 
			     SparseAdjacency &G, 
			     int k, 
			     int *MapToG
			   ) {

  // Remat: the kth node is connected at least once to 
  // the k-1 preceeding nodes (the motif is connex)

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif

  if ( (k) == motifSup.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //
# ifdef MSG
    cerr << "  New motif is stored " << endl;
# endif
    if( _OutputMode == CountMode ) {
      _Count++;

      // 
      // printList( " MapToG : ", MapToG,  &MapToG[k+1], true);
      
    } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
      //_Count += writeOccurences( G, motif, motif, MapToG );   

      // TODO storage
      _FoundBList -> storeData( MapToG );
      _Count++;
    }
  } else {
    

    int MapToG_u;

    /////////////////////////////////////////////////////////
    //
    // Get Incoming and Outcomming edges u: u -> k+1 or u <- k+1
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    
    const int *Vin_motif_kpp      = motifInf.getRow(k+1);
    const int *Vin_motif_kpp_end  = motifInf.getRowEnd(k+1);
          int Vin_motif_kpp_size  = Vin_motif_kpp_end - Vin_motif_kpp;
    const int *Vout_motif_kpp     = motifSup.getCol(k+1);
    const int *Vout_motif_kpp_end = motifSup.getColEnd(k+1);
          int Vout_motif_kpp_size = Vout_motif_kpp_end - Vout_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Incoming {u} = Vin_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vin_motif_kpp, Vin_motif_kpp_end, true );

    cerr << "  Motif-Outcoming {u} = Vout_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vout_motif_kpp, Vout_motif_kpp_end, true );
#   endif
    
    // Make the intersection :
    //    INTER{ Vout_G( MapToG[ u=Vin_motif(k+1), with (u<=k+1) ] ) 
    //       and Vin_G( MapToG[ u=Vout_motif(k+1), with (u<=k+1) ] }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
      
    //
    // Intersection initialisation with Vout_G( MapToG[ u=0] )
    //
    const int *it;

    // Set to intersect
    const int **L = new 
      const int*[Vin_motif_kpp_size + Vout_motif_kpp_size];
    int       *TL = 
      new int[Vin_motif_kpp_size + Vout_motif_kpp_size];

    //
    // Store Vout_G( MapToG[ u=Vin_motif(k+1) with (u<=k+1) ]) set 
    // in the sorted heap.
    //
    int i = 0;
    for ( it = Vin_motif_kpp; (it != Vin_motif_kpp_end); it++) {
      MapToG_u =  MapToG[*it];
      TL[i] = G.getColSize( MapToG_u );
      if( TL[i] != 0) {
	L[i++]    = G.getCol( MapToG_u );
      } else {
	// Empty set
	i = 0;
	break;
      }
    }
    
    // Test if incoming arcs are required in the pattern (  Vin_motif_kpp_size )
    // and no outcomming arcs have found ( G.getCol( MapToG_u ) )
    if ( (i != 0) || ( Vin_motif_kpp_size == 0)) {
      //
      // Store Vin_G( MapToG[ u=Vout_motif(k+1) with (u<=k+1) ]) set 
      // in the sorted heap.
      //
      int j = i;
      for ( it = Vout_motif_kpp; (it != Vout_motif_kpp_end); it++) {
	MapToG_u =  MapToG[*it];
	TL[j] = G.getRowSize( MapToG_u );
	if( TL[j] != 0) {
	  L[j++]    = G.getRow( MapToG_u );
	} else {
	  // Empty set
	  j = i;
	}
      }

      // Build the Heap
      // TODO : Heap deal with no list (i == 0)  to avoid the test
      // Test if outcoming arcs are required in the pattern (  Vout_motif_kpp_size )
      // and no incomming arcs have found ( G.getRow( MapToG_u ) )
      if( (j > 0) && ((i != j) || ( Vout_motif_kpp_size == 0)) ) {
	Heap H(L, TL, j );

#     ifdef MSG
	H.Debug();
#     endif

	//
	// Get the intersection list.
	//
	int a = H.NextCommonElement();
	while (a != NOT_INDEX) {
	  
	  // Apply constraints 
	  if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	      || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {
	  
	    if (!IsIn( MapToG, k+1, a)) { 

	      // New sub-motif found
	
	      MapToG[k+1] = a;
	      kpp_MatchDirectedMotif( motifSup, motifInf, G, k+1, 
				      MapToG );
	    }
	  }
	  a = H.NextCommonElement();
	}
      }
    }
    delete[] L;
    delete[] TL;

    
  }
}

void FindMotif::kpp_MatchUndirectedColoredMotif( SparseAdjacency &motif,
				      const int *MotifColor,
				      SparseAdjacency &G, 
				      const int *GColor, 
				      int k, 
				      int* MapToG 
				    ) {
  //
  // Undirected case
  //
    
  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif

  // TODO  change  motif.getNbrRow() by a value
  if ( (k) == motif.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

#   ifdef MSG
    cerr << "  New motif is stored " << endl;
#   endif

    switch( _OutputMode ) {

    case WriteMode : 
      for ( int i = 0; i < motif.getNbrRow(); i++){
	*(_OutputStream) <<  MapToG[i] << " " ;
      }
      *(_OutputStream) << endl;
      break;

    case StoreMode : 
      {
      
	// Allocate and Copy motif
	int *new_motif = new int[ motif.getNbrRow() ] ;
	for ( int i = 0; i < motif.getNbrRow(); i++){
	  new_motif[i] = MapToG[i];
	}
	// Add in list
	_FoundList->push_back( new_motif );
      }
      break;

    case TCMode :
      {
	int motif_size = motif.getNbrRow();	
	for ( int i = 0; i < motif_size; i++){
	  _TopoClasses[ i*motif_size +  MapToG[i] ] = true; 
	}
      }
      break;

    default :
      // Count Mode
      break;
    }  
    _Count++;

  } else {


    /////////////////////////////////////////////////////////
    //
    // Get neighbours u of k+1 in the motif
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    const int *V_motif_kpp      = motif.getCol(k+1);
    const int *V_motif_kpp_end  = motif.getColEnd(k+1);
    int V_motif_kpp_size = V_motif_kpp_end - V_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Neighbours  {u} = V_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", V_motif_kpp, V_motif_kpp_end, true );
#   endif

    //
    // Make the intersection :
    //    INTER{ V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ] ) }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
    //

    //
    // Intersection initialisation with V_G( MapToG (u=0) )
    //
    const int *it = V_motif_kpp;

    // Node ID in G
    // int  MapToG_u = MapToG[ *it ];

    // Set to intersect
    const int **L = new const int*[V_motif_kpp_size];
    int       *TL = new int[V_motif_kpp_size];

    //
    // Store V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ]) set 
    // sorted heap.
    //
    int i = 0;
    for ( it = V_motif_kpp; (it != V_motif_kpp_end); it++) {
	L[i]    = G.getAllRow( MapToG[*it] );
	TL[i++] = G.getAllRowSize( MapToG[*it]);
    }
    Heap H(L, TL, V_motif_kpp_size );

#   ifdef MSG
    H.Debug();
#   endif

    //
    // Get the intersection list.
    //
    int a = H.NextCommonElement();
    while (a != NOT_INDEX) {

      // Apply constraints 
      if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	  || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {

	// Same Color
	if ( GColor[ a ] == MotifColor[k+1] ) {

	  if (!IsIn( MapToG, k+1, a)) {  
	    
	    // New sub-motif found
	    MapToG[k+1] = a;
	    kpp_MatchUndirectedColoredMotif( motif, MotifColor, G, GColor, 
					     k+1, MapToG );
	  }
	}
      }
      a = H.NextCommonElement();
    }
    delete[] L;
    delete[] TL;

  }
}

void FindMotif::kpp_MatchDirectedColoredMotif( 
				    SparseAdjacency &motifSup, 
				    SparseAdjacency &motifInf,
				    const int *MotifColor,
				    SparseAdjacency &G, 
				    const int *GColor,
				    int k, 
				    int *MapToG
				  ) {
  //
  // Directed part
  //
    

  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif

  if ( (k) == motifSup.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

# ifdef MSG
    cerr << "  New motif is stored " << endl;
# endif

    switch( _OutputMode ) {

    case WriteMode : 
      for ( int i = 0; i < motifSup.getNbrRow(); i++){
	*(_OutputStream) <<  MapToG[i] << " " ;
      }
      *(_OutputStream) << endl;
      break;

    case StoreMode : 
      {
	// Allocate and Copy motif
	int *new_motif = new int[ motifSup.getNbrRow() ] ;
	for ( int i = 0; i < motifSup.getNbrRow(); i++){
	  new_motif[i] = MapToG[i];
	}
	// Add in list
	_FoundList->push_back( new_motif );
      }
      break;

    case TCMode : 
      {
	int motif_size = motifSup.getNbrRow();	
	for ( int i = 0; i < motif_size; i++){
	  _TopoClasses[ i*motif_size +  MapToG[i] ] = true; 
	}
      }
      break;

    default :
      // Count Mode
      break;
    }  
    _Count++;
    
  } else {
    

    int MapToG_u;

    /////////////////////////////////////////////////////////
    //
    // Get Incoming and Outcomming edges u: u -> k+1 or u <- k+1
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    
    const int *Vin_motif_kpp      = motifInf.getRow(k+1);
    const int *Vin_motif_kpp_end  = motifInf.getRowEnd(k+1);
          int Vin_motif_kpp_size  = Vin_motif_kpp_end - Vin_motif_kpp;
    const int *Vout_motif_kpp     = motifSup.getCol(k+1);
    const int *Vout_motif_kpp_end = motifSup.getColEnd(k+1);
          int Vout_motif_kpp_size = Vout_motif_kpp_end - Vout_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Incoming {u} = Vin_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vin_motif_kpp, Vin_motif_kpp_end, true );

    cerr << "  Motif-Outcoming {u} = Vout_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vout_motif_kpp, Vout_motif_kpp_end, true );
#   endif
    
    // Make the intersection :
    //    INTER{ Vout_G( MapToG[ u=Vin_motif(k+1), with (u<=k+1) ] ) 
    //       and Vin_G( MapToG[ u=Vout_motif(k+1), with (u<=k+1) ] }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
      
    //
    // Intersection initialisation with Vout_G( MapToG[ u=0] )
    //
    const int *it;

    // Set to intersect
    const int **L = new 
      const int*[Vin_motif_kpp_size + Vout_motif_kpp_size];
    int       *TL = 
      new int[Vin_motif_kpp_size + Vout_motif_kpp_size];

    //
    // Store Vout_G( MapToG[ u=Vin_motif(k+1) with (u<=k+1) ]) set 
    // in the sorted heap.
    //
    int i = 0;
    for ( it = Vin_motif_kpp; (it != Vin_motif_kpp_end); it++) {
      MapToG_u =  MapToG[*it];
      TL[i] = G.getColSize( MapToG_u );
      if( TL[i] != 0) {
	L[i++]    = G.getCol( MapToG_u );
      } else {
	// Empty set
	i = 0;
	break;
      }
    }
    
    // Test if incoming arcs are required in the pattern (  Vin_motif_kpp_size )
    // and no outcomming arcs have found ( G.getCol( MapToG_u ) )
    if ( (i != 0) || ( Vin_motif_kpp_size == 0)) {
      //
      // Store Vin_G( MapToG[ u=Vout_motif(k+1) with (u<=k+1) ]) set 
      // in the sorted heap.
      //
      int j = i;
      for ( it = Vout_motif_kpp; (it != Vout_motif_kpp_end); it++) {
	MapToG_u =  MapToG[*it];
	TL[j] = G.getRowSize( MapToG_u );
	if( TL[j] != 0) {
	  L[j++]    = G.getRow( MapToG_u );
	} else {
	  // Empty set
	  j = i;
	}
      }

      // Build the Heap
      // TODO : Heap deal with no list (i == 0)  to avoid the test
      // Test if outcoming arcs are required in the pattern (  Vout_motif_kpp_size )
      // and no incomming arcs have found ( G.getRow( MapToG_u ) )
      if( (j > 0) && ((i != j) || ( Vout_motif_kpp_size == 0)) ) {
	Heap H(L, TL, j );

#     ifdef MSG
	H.Debug();
#     endif

	//
	// Get the intersection list.
	//
	int a = H.NextCommonElement();
	while (a != NOT_INDEX) {
	  
	  // Apply constraints 
	  if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	      || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {

	    // Same Color
	    if ( GColor[ a ] == MotifColor[k+1] ) {

	      if (!IsIn( MapToG, k+1, a)) { 

		// New sub-motif found
		
		MapToG[k+1] = a;
		kpp_MatchDirectedColoredMotif(  motifSup, motifInf, MotifColor, 
						G, GColor,
						k+1, MapToG );
	      }
	    }
	  }
	  a = H.NextCommonElement();
	}
      }
    }
    delete[] L;
    delete[] TL;
    
  }
}


void FindMotif::kpp_MatchDirectedColoredMotif( 
				    SparseAdjacency &motifSup, 
				    SparseAdjacency &motifInf,
				    SparseAdjacency &cmplSup, 
				    SparseAdjacency &cmplInf,
				    const int *MotifColor,
				    SparseAdjacency &G, 
				    const int *GColor,
				    int k, 
				    int *MapToG
				  ) {
  //
  // Directed part
  //
    

  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif

  if ( (k) == motifSup.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

# ifdef MSG
    cerr << "  New motif is stored " << endl;
    printList("MaptoG ok: ", MapToG, &MapToG[k+1], true );    
# endif

    switch( _OutputMode ) {

    case WriteMode : 
      for ( int i = 0; i < motifSup.getNbrRow(); i++){
	*(_OutputStream) <<  MapToG[i] << " " ;
      }
      *(_OutputStream) << endl;
      break;

    case StoreMode : 
      {
	// Allocate and Copy motif
	int *new_motif = new int[ motifSup.getNbrRow() ] ;
	for ( int i = 0; i < motifSup.getNbrRow(); i++){
	  new_motif[i] = MapToG[i];
	}
	// Add in list
	_FoundList->push_back( new_motif );
      }
      break;

    case TCMode : 
      {
	int motif_size = motifSup.getNbrRow();	
	for ( int i = 0; i < motif_size; i++){
	  _TopoClasses[ i*motif_size +  MapToG[i] ] = true; 
	}
      }
      break;

    default :
      // Count Mode
      break;
    }  
    _Count++;
    
  } else {
    

    int MapToG_u;

    /////////////////////////////////////////////////////////
    //
    // Get Incoming and Outcomming edges u: u -> k+1 or u <- k+1
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    
    const int *Vin_motif_kpp      = motifInf.getRow(k+1);
    const int *Vin_motif_kpp_end  = motifInf.getRowEnd(k+1);
          int Vin_motif_kpp_size  = Vin_motif_kpp_end - Vin_motif_kpp;
    const int *Vout_motif_kpp     = motifSup.getCol(k+1);
    const int *Vout_motif_kpp_end = motifSup.getColEnd(k+1);
          int Vout_motif_kpp_size = Vout_motif_kpp_end - Vout_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Incoming {u} = Vin_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vin_motif_kpp, Vin_motif_kpp_end, true );

    cerr << "  Motif-Outcoming {u} = Vout_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vout_motif_kpp, Vout_motif_kpp_end, true );
#   endif

    // printList("MaptoG : ", MapToG, &MapToG[k+1], true );    

    // Make the intersection :
    //    INTER{ Vout_G( MapToG[ u=Vin_motif(k+1), with (u<=k+1) ] ) 
    //       and Vin_G( MapToG[ u=Vout_motif(k+1), with (u<=k+1) ] }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
      
    //
    // Intersection initialisation with Vout_G( MapToG[ u=0] )
    //
    const int *it;

    // Set to intersect
    const int **L = new 
      const int*[Vin_motif_kpp_size + Vout_motif_kpp_size];
    int       *TL = 
      new int[Vin_motif_kpp_size + Vout_motif_kpp_size];

    //
    // Store Vout_G( MapToG[ u=Vin_motif(k+1) with (u<=k+1) ]) set 
    // in the sorted heap.
    //
    int i = 0;
    for ( it = Vin_motif_kpp; (it != Vin_motif_kpp_end); it++) {
      MapToG_u =  MapToG[*it];
      TL[i] = G.getColSize( MapToG_u );
      if( TL[i] != 0) {
	L[i++]    = G.getCol( MapToG_u );
      } else {
	// Empty set
	i = 0;
	break;
      }
    }
    
    // Test if incoming arcs are required in the pattern (  Vin_motif_kpp_size )
    // and no outcomming arcs have found ( G.getCol( MapToG_u ) )
    if ( (i != 0) || ( Vin_motif_kpp_size == 0)) {
      //
      // Store Vin_G( MapToG[ u=Vout_motif(k+1) with (u<=k+1) ]) set 
      // in the sorted heap.
      //
      int j = i;
      for ( it = Vout_motif_kpp; (it != Vout_motif_kpp_end); it++) {
	MapToG_u =  MapToG[*it];
	TL[j] = G.getRowSize( MapToG_u );
	if( TL[j] != 0) {
	  L[j++]    = G.getRow( MapToG_u );
	} else {
	  // Empty set
	  j = i;
	}
      }

      // Build the Heap
      // TODO : Heap deal with no list (i == 0)  to avoid the test
      // Test if outcoming arcs are required in the pattern (  Vout_motif_kpp_size )
      // and no incomming arcs have found ( G.getRow( MapToG_u ) )
      if( (j > 0) && ((i != j) || ( Vout_motif_kpp_size == 0)) ) {
	Heap H(L, TL, j );

#     ifdef MSG
	H.Debug();
#     endif

	//
	// TODO with heap
	//

	const int *V_cmpl_in_kpp      = cmplInf.getRow(k+1);
	const int *V_cmpl_in_kpp_end  = cmplInf.getRowEnd(k+1);

	const int *V_cmpl_out_kpp      = cmplSup.getCol(k+1);
	const int *V_cmpl_out_kpp_end  = cmplSup.getColEnd(k+1);

	//
	// Get the intersection list.
	//
	int a = H.NextCommonElement();
	while (a != NOT_INDEX) {
	  
	  // Apply constraints 
	  if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	      || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {

	    // Same Color
	    if ( GColor[ a ] == MotifColor[k+1] ) {

	      if (!IsIn( MapToG, k+1, a)) { 

		// ?? cout << " node selected " << a << ", k=" << k << endl;
		bool not_found = true;
		for ( it = V_cmpl_in_kpp; (it != V_cmpl_in_kpp_end) && (not_found); 
		      it++) {
		  // G.getAllRow( MapToG[*it] );
		  not_found = ! G.isInCol( MapToG[*it], a );
		  // ??? cout << "   A. rentrantes : " <<  MapToG[*it] << endl;
		  // ??? printList("   G sortantes ", G.getCol( MapToG[*it] ),
		  // ???	    G.getColEnd( MapToG[*it]) , true);
		}
		for ( it = V_cmpl_out_kpp; (it != V_cmpl_out_kpp_end) && (not_found); 
		      it++) {
		  // G.getAllRow( MapToG[*it] );
		  not_found = ! G.isInRow( MapToG[*it], a );
		   // ??? cout << cout << "   A. sortantes : " <<  MapToG[*it] << endl;
		   // ??? cout << printList("   G rentrantes ", G.getRow( MapToG[*it] ),
		  // ??? cout << G.getRowEnd( MapToG[*it]) , true);
		}

		// New sub-motif found
		if( not_found ) {	
		  MapToG[k+1] = a;
		  kpp_MatchDirectedColoredMotif(  motifSup, motifInf,
						  cmplSup, cmplInf,
						  MotifColor, 
						  G, GColor,
						  k+1, MapToG );
		}
	      }
	    }
	  }
	  a = H.NextCommonElement();
	}
      }
    }
    delete[] L;
    delete[] TL;
    
  }
}

void FindMotif::kpp_MatchDirectedMotif( 
				    SparseAdjacency &motifSup, 
				    SparseAdjacency &motifInf,
				    SparseAdjacency &cmplSup, 
				    SparseAdjacency &cmplInf,
				    SparseAdjacency &G, 
				    int k, 
				    int *MapToG
				  ) {
  //
  // Directed part
  //
    

  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif

  if ( (k) == motifSup.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

# ifdef MSG
    cerr << "  New motif is stored " << endl;
# endif

    switch( _OutputMode ) {

    case WriteMode : 
      for ( int i = 0; i < motifSup.getNbrRow(); i++){
	*(_OutputStream) <<  MapToG[i] << " " ;
      }
      *(_OutputStream) << endl;
      break;

    case StoreMode : 
      {
	// Allocate and Copy motif
	int *new_motif = new int[ motifSup.getNbrRow() ] ;
	for ( int i = 0; i < motifSup.getNbrRow(); i++){
	  new_motif[i] = MapToG[i];
	}
	// Add in list
	_FoundList->push_back( new_motif );
      }
      break;

    case TCMode : 
      {
	int motif_size = motifSup.getNbrRow();	
	for ( int i = 0; i < motif_size; i++){
	  _TopoClasses[ i*motif_size +  MapToG[i] ] = true; 
	}
      }
      break;

    default :
      // Count Mode
      break;
    }  
    _Count++;
    
  } else {
    

    int MapToG_u;

    /////////////////////////////////////////////////////////
    //
    // Get Incoming and Outcomming edges u: u -> k+1 or u <- k+1
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    
    const int *Vin_motif_kpp      = motifInf.getRow(k+1);
    const int *Vin_motif_kpp_end  = motifInf.getRowEnd(k+1);
          int Vin_motif_kpp_size  = Vin_motif_kpp_end - Vin_motif_kpp;
    const int *Vout_motif_kpp     = motifSup.getCol(k+1);
    const int *Vout_motif_kpp_end = motifSup.getColEnd(k+1);
          int Vout_motif_kpp_size = Vout_motif_kpp_end - Vout_motif_kpp;

#   ifdef MSG
    cerr << "  Motif-Incoming {u} = Vin_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vin_motif_kpp, Vin_motif_kpp_end, true );

    cerr << "  Motif-Outcoming {u} = Vout_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vout_motif_kpp, Vout_motif_kpp_end, true );
#   endif
    
    // Make the intersection :
    //    INTER{ Vout_G( MapToG[ u=Vin_motif(k+1), with (u<=k+1) ] ) 
    //       and Vin_G( MapToG[ u=Vout_motif(k+1), with (u<=k+1) ] }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
      
    //
    // Intersection initialisation with Vout_G( MapToG[ u=0] )
    //
    const int *it;

    // Set to intersect
    const int **L = new 
      const int*[Vin_motif_kpp_size + Vout_motif_kpp_size];
    int       *TL = 
      new int[Vin_motif_kpp_size + Vout_motif_kpp_size];

    //
    // Store Vout_G( MapToG[ u=Vin_motif(k+1) with (u<=k+1) ]) set 
    // in the sorted heap.
    //
    int i = 0;
    for ( it = Vin_motif_kpp; (it != Vin_motif_kpp_end); it++) {
      MapToG_u =  MapToG[*it];
      TL[i] = G.getColSize( MapToG_u );
      if( TL[i] != 0) {
	L[i++]    = G.getCol( MapToG_u );
      } else {
	// Empty set
	i = 0;
	break;
      }
    }
    
    // Test if incoming arcs are required in the pattern (  Vin_motif_kpp_size )
    // and no outcomming arcs have found ( G.getCol( MapToG_u ) )
    if ( (i != 0) || ( Vin_motif_kpp_size == 0)) {
      //
      // Store Vin_G( MapToG[ u=Vout_motif(k+1) with (u<=k+1) ]) set 
      // in the sorted heap.
      //
      int j = i;
      for ( it = Vout_motif_kpp; (it != Vout_motif_kpp_end); it++) {
	MapToG_u =  MapToG[*it];
	TL[j] = G.getRowSize( MapToG_u );
	if( TL[j] != 0) {
	  L[j++]    = G.getRow( MapToG_u );
	} else {
	  // Empty set
	  j = i;
	}
      }

      // Build the Heap
      // TODO : Heap deal with no list (i == 0)  to avoid the test
      // Test if outcoming arcs are required in the pattern (  Vout_motif_kpp_size )
      // and no incomming arcs have found ( G.getRow( MapToG_u ) )
      if( (j > 0) && ((i != j) || ( Vout_motif_kpp_size == 0)) ) {
	Heap H(L, TL, j );

#     ifdef MSG
	H.Debug();
#     endif

	//
	// TODO with heap
	//

	const int *V_cmpl_in_kpp      = cmplInf.getRow(k+1);
	const int *V_cmpl_in_kpp_end  = cmplInf.getRowEnd(k+1);

	const int *V_cmpl_out_kpp      = cmplSup.getCol(k+1);
	const int *V_cmpl_out_kpp_end  = cmplSup.getColEnd(k+1);

	//
	// Get the intersection list.
	//
	int a = H.NextCommonElement();
	while (a != NOT_INDEX) {
	  
	  // Apply constraints 
	  if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	      || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {


	      if (!IsIn( MapToG, k+1, a)) { 

		bool not_found = true;
		for ( it = V_cmpl_in_kpp; (it != V_cmpl_in_kpp_end) && (not_found); 
		      it++) {
		  // G.getAllRow( MapToG[*it] );
		  not_found = ! G.isInCol( MapToG[*it], a );
		}
		for ( it = V_cmpl_out_kpp; (it != V_cmpl_out_kpp_end) && (not_found); 
		      it++) {
		  // G.getAllRow( MapToG[*it] );
		  not_found = ! G.isInRow( MapToG[*it], a );
		}

		// New sub-motif found
		if( not_found ) {	
		  MapToG[k+1] = a;
		  kpp_MatchDirectedMotif(  motifSup, motifInf,
					   cmplSup, cmplInf,
					   G, 
					   k+1, MapToG );
		}
	      }
	  }
	  a = H.NextCommonElement();
	}
      }
    }
    delete[] L;
    delete[] TL;
    
  }
}

// TODO store in memory lists
void findNextGeneration( SparseAdjacency &G,    
			 vector< SparseAdjacency *> &AncestorMotifs, 
			 vector< vector< vector<int> > *> &AncestorCT,
			 FindMotif &find, int motif_size, int nb_edges,
			 vector <int *> &motifs, vector <int64_t> &counts, 
			 int generation_mode=0, 
			 int search_mode = 0 ) {



  // motifs, counts, motif list found with the number of occurence  
  vector< SparseAdjacency *> ChildMotifs;
  vector< vector< vector<int> > *> ChildCT;

  string buf[3];

  if ( AncestorMotifs.size() == 0 ) {
    
    // Version Michel ???
    /*
    if( (search_mode == ShrinkMode) && (G.getNodeValues() == 0) ) {
      cerr << __FILE__ << " Graph not shrunk" << endl;
      return ;
    } 
    */

    if( generation_mode == ComputeMotifs || generation_mode == WriteMotifs ) {
      //
      // Generate clique
      //
      SparseAdjacency *motif = new SparseAdjacency( motif_size, "clique", 
						    G.getSymmetry() );
      vector< vector<int> > *CT = find.getTopologicClasses( *motif );
      
      // Print automorphisms
      // cout << " Nbr automorphisms " << computeNbrOfAutomorphism( *motif ) << endl;

      ChildCT.push_back( CT );
      ChildMotifs.push_back( motif );
      nb_edges =  motif_size * (motif_size-1);
      if ( G.getSymmetry() )
	nb_edges = nb_edges / 2;
      
      if( generation_mode == WriteMotifs ) {
	char sym = 'u';
	if( ! G.getSymmetry() ) sym = 'd'; 
	
	sprintf ( MotifFileName, PatternMotifFileName, 
		  MotifDir, sym, motif_size);
	fmout.open(MotifFileName);
	if( ! fmout.is_open() ) {
#   ifdef MSG
	  cerr << "Can't open motif DB file : " << MotifFileName << endl; 
#   endif
	  return;
	}
	fmout << "Motif size: " << motif_size << endl; 
	fmout << endl; 
      }
    } else if( generation_mode == ReadMotifs ) {
      char sym = 'u';
      if( ! G.getSymmetry() ) sym = 'd'; 
      
      sprintf ( MotifFileName, PatternMotifFileName, 
		MotifDir, sym, motif_size);
      fmin.open(MotifFileName);
      if( ! fmin.good() ) {
#   ifdef MSG 
	cerr << "Can't open motif DB file : " << MotifFileName << endl; 
#   endif
	return;
      }
      fmin >> buf[0] >> ws >> buf[1] >> ws >> buf[2];
      // cerr << buf[0] << " " <<  buf[1] << " " <<  buf[2] << endl;
      fmin >> buf[0] >> ws >> buf[1] >> ws;
      // cerr << buf[0] << " " <<  buf[1] << endl;

      SparseAdjacency *CMotif = new SparseAdjacency(fmin, "", G.getSymmetry());
      ChildMotifs.push_back( CMotif );

    }
  } else {

    if( generation_mode == ReadMotifs ) {

      fmin >> buf[0] >> ws >> buf[1] >> ws;
      // cerr << buf[0] << " " <<  buf[1] << endl;

      SparseAdjacency *CMotif = new SparseAdjacency(fmin, "", G.getSymmetry());
      if( CMotif -> getNbrRow() )
	ChildMotifs.push_back( CMotif );

    } else {

      getChildMotifs( AncestorMotifs, AncestorCT, 
		      ChildMotifs, ChildCT, 
		      motif_size, nb_edges, G.getSymmetry(), 
		      find );
    }
  }
  AncestorMotifs.clear();
  AncestorCT.clear();
  
  // Configure FindMotif
  find.setOutputMode( CountMode );

  for( size_t i=0; i < ChildMotifs.size(); i++) {
    // TO Optimize pass ChildCT[i]
    // ChildMotifs[i]->print( "new generation ", "all" );

    // unused 1
    // if ( root ) {
    //  addMotifConstraints( root, root, *ChildMotifs[i], getUniqPath( *ChildMotifs[i] ), 1);
    // } else {
    getTime( 0 );
      
    // ChildMotifs[i]->print("motif ");

    // Print automorphisms
    // cout << " Nbr automorphisms " << computeNbrOfAutomorphism( *ChildMotifs[i] ) << endl;

    if ( search_mode == NonExactMode ){
      find.find( G, *ChildMotifs[i] );
    } else if ( search_mode == ShrinkMode ) {
      find.findShrink( G, *ChildMotifs[i] );
    } else if( search_mode == ExactMode ) {
      find.findExactly( G, *ChildMotifs[i] );
    } else if( search_mode == NullMode ) {
      
    }
    
    if( generation_mode == WriteMode ) {
      fmout << "Motif number written : " <<  _NbrOfMotifs << endl;
      ChildMotifs[i]->write( fmout );
    }
    // Unused 1 (see before)
    // } 
    // Store counts and motif
    int64_t c = find.getCount();
    // cout << " count " << c << endl;
    
    // Used to select only motif found
    // if ( c ) {
    counts.push_back( c );
    motifs.push_back( ChildMotifs[i] -> getIndexes( ) );
    // }
    getTime( 1 );
  
    _CountAllMotifs += find.getCount() ;
    
    _NbrOfMotifs = _NbrOfMotifs + 1;
  
#   ifdef MSG
    cout << " -> Motif " << _NbrOfMotifs << " : " << find.getCount() << endl;
#   endif
  
    // cout << " -> Motif " << _NbrOfMotifs << " : " << find.getCount() << endl;
    
    find.deallocateFoundBList();
  }

  find.setOutputMode( CountMode );

  if( ChildMotifs.size() != 0 ) {
    findNextGeneration( G, ChildMotifs, ChildCT, 
			find, motif_size, nb_edges -1, 
			motifs, counts, 
			generation_mode, search_mode
			);
  } else if ( generation_mode == WriteMotifs ) {
    fmout << endl;
    fmout.close();
  } else if ( generation_mode == ReadMotifs ) {
    fmin.close();
  }
# ifdef MSG
  cerr << " delete motif"  << endl;
# endif

}


void buildAndfindAllMotif( SparseAdjacency &G, int k, FindMotif &find,
			   vector< int * > &motifs, 
			   vector<int64_t> &counts,
			   int generate_mode,
			   int search_mode ) {

  vector< SparseAdjacency *> motifList;
  vector< vector< vector<int> >*> CTList;
  
  int nb_edges = 0;
  findNextGeneration( G, motifList, CTList, find, k, 
		      nb_edges, motifs, counts, generate_mode, search_mode );

# ifdef MSG
  cout << "Count all motifs " << _CountAllMotifs << endl;
  cout << "Nbr of motifs " << _NbrOfMotifs << endl;
  cout << "Nbr of motifs " << counts.size() << " " << motifs.size() << endl;
# endif
  
}

void FindMotif::writeList( const char *output_name, 
			   const vector<int*> *found, 
			   const int motif_size, 
			   bool block_alloc ) {
  
    
  // save output

# ifdef MSG 
  bool use_cout = false;
  ofstream *file = 0;

  streambuf* cout_buffer = cout.rdbuf();
  
  // Open 
  if ( strlen( output_name ) == 0) {
    use_cout =true;
  } else {
    file = new ofstream( output_name );
    // redirect file to cout
    cout.rdbuf ( file->rdbuf() );
  }

  if ( block_alloc ) {

    int *ptr;
    for( vector<int*>::const_iterator it = found->begin();
	 (it != found->end() ); it++) {
      ptr = *it;
      for(; *ptr != 0 ; ){
	for(int i=0; i < motif_size ; i++){
	  cout << " " << *ptr;
	  ptr++;
	}
	cout << endl;
      }
    }
  } else {

    int *ptr = 0;
    for( vector<int*>::const_iterator it = found->begin();
	 (it != found->end() ); it++) {
      for(int i=0; i < motif_size ; i++){
	cout << " " << *ptr;
	ptr++;
      }
      cout << endl;
    }
  }


    if( ! use_cout ) {
    file->close();
  
    // restore old output
    cout.rdbuf (cout_buffer);
  }
# endif

}

// High level method
// Change internal state of the object
// Find all motifs in a graph (Fanmod type)
// In test
vector< int* > *FindMotif::findAllMotifs(  
			   SparseAdjacency &G, 
			   int MotifSize
			   ) {

  _Count = 0;

  // Motif
  int *MapToG = new int[ MotifSize];
  int GraphSize = G.getNbrRow(); 
# ifdef MSG 
  G.print("Network"); 
# endif

  int *Neighbours = new int[ GraphSize ];
  int *StartNeighbours = Neighbours;
  // int *EndNeighbours = Neighbours;

  if ( G.getSymmetry() ) {

    for ( int i = 0; i < GraphSize; i++ ) {

      MapToG[0] = i;

      if( G.getAllRowSize( i ) != 0) {
	  // int *Start_n_ = StartNeighbours;
	  int *End_n_   = StartNeighbours;
	  int a =  MapToG[0];
	  // Add new neighbours if > a and update EndNeighbours
	  for( const int *it_n = G.getAllRow( a ), 
		 *it_n_end = G.getAllRowEnd( a ); 
	       it_n != it_n_end; it_n++) {
	    if( *it_n > a ) {
	      *End_n_ = *it_n;
	      End_n_++;
	    }
	  } 

#         ifdef MSG
	  printList( " MapToG", MapToG, &MapToG[1], true);

 	  printList( " Neigh", Start_n_, End_n_, true);
#         endif

	  // FanMod algorithm
	  kpp_MatchUndirectedMotifs( MotifSize, G, 0,  
				     MapToG, Neighbours, End_n_
			);

      }
    }
  } else {
#   ifdef MSG     
     cerr << "The graph must be symmetric " << endl;
#   endif
  } 
   vector< int*> *found_list=_FoundList;

   delete [] MapToG ;
   delete [] Neighbours;

  _FoundList = 0;
  return( found_list ); 
}

// Low level routine
// Fanmod type algorithm
// In test
void FindMotif::kpp_MatchUndirectedMotifs( 
					  int MotifSize,
					  SparseAdjacency &G, 
					  int k, 
					  int* MapToG,
					  int* StartNeighbours,
					  int* EndNeighbours
					   ) {
  //
  // Undirected case
  //
    
  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes
    
# ifdef MSG
    cerr << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif


    if ( (k) == (MotifSize-1) ) {
      //
      // A new motif is found
      //
      
#   ifdef MSG
      cerr << "  New motif is stored " << endl;
#   endif
      
      if( _OutputMode == CountMode ) {
	_Count++;
# ifdef MSG
	printList( "motif ", MapToG, &MapToG[k+1], true) ;
# endif 
	
      } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
	//      _Count += writeOccurences( G, motif, motif, MapToG );   
	
	// TODO storage
	//_FoundBList -> storeData( MapToG );
	//_FoundList->push_back( Map );

	// Allocate a new block if necessary
	// in motif size unit
	if( (_Count % ChunkSize ) == 0) {
	  _ChunkStorage= new int[ChunkSize*MotifSize];
	}

	// Copy Map
	int *save_ptr=_ChunkStorage;
	for( int *it=MapToG, *it_end=&MapToG[MotifSize]; it != it_end;) {
	  *_ChunkStorage++ = *it++;
	}
	_FoundList->push_back( save_ptr );
	_Count++;
      }
      
    } else {
      
      
      /////////////////////////////////////////////////////////
      //
      // Get neighbours u of k+1 in the motif
      // with u <= k+1 node of the  pattern
      //
      /////////////////////////////////////////////////////////
      int *End_;

      if( (StartNeighbours != EndNeighbours) || (k == 0) ) {
	
	// Neighbors could be the all the graph
	int *Start_ = new int[ G.getNbrRow() ];

	int a = MapToG[0];

	int* it = StartNeighbours, *it_end =  EndNeighbours;
	// printList( "  Vext", it, it_end, true);

	// 
	// For all v in V_ext
	for( ; it != it_end ; it++) {

	  // cout    << "  k = " << k << "test : " <<  *it << endl;
	  MapToG[k+1] = *it;

	  // Copy in Vext' the vertices not done
	  int *it_dst = Start_;
	  for ( int *it_c = (it+1) ;
		it_c != it_end; it_c++, it_dst++){
	    *it_dst = *it_c;
	  }

	  End_ = it_dst;

#         ifdef MSG
	  printList( "  Vext'", Start_, End_, true);
#	  endif

#   ifdef MSG
	  printList( "  MapToG", MapToG, &MapToG[k+2], true);
#   endif

	  int *End_n_ = End_;

	  // Prepare New Vextensions 
	  // Add new neighbours if > a and update EndNeighbours
	  for( const int *it_n = G.getAllRow( MapToG[k+1] ), 
		 *it_n_end = G.getAllRowEnd(  MapToG[k+1] ); 
	       it_n != it_n_end; it_n++) {

	    // TODO : Optimize
	    if( *it_n > a  && 
		( !IsInNeighboring( G, MapToG, &MapToG[k+1], *it_n)) ) {
	      *End_n_ = *it_n;
	      End_n_++;
	    }
	  } 

	  // Sort them
	  sort( Start_, End_n_ );
	  End_n_ = unique (Start_, End_n_);

#         ifdef MSG
	  printList( "  Vext'", Start_, End_n_, true);
#	  endif
	  
	  // if ( (it+1) != End_ )
	  kpp_MatchUndirectedMotifs( MotifSize, G, k+1, MapToG, Start_, End_n_);
	  
	}
	delete [] Start_ ;
      }
    }
}

// ??? motif is a connexe graph.
// High level method
// Change internal state of the object

vector< int* > *FindMotif::findShrink(  
			   SparseAdjacency &G, 
			   SparseAdjacency &motif 
			   ) {

  allocateFoundList();

  _Count = 0;

  int G_NbrNodes = G.getNbrRow();

  // Remove poor neighbouring
  _RealMotifSize = motif.getNbrRow();

  // Mapping from shrunk motif node to original node
  // use to keep the original motif constraints
  int *map = new int[_RealMotifSize];
  int *map_1 = new int[_RealMotifSize];

  SparseAdjacency ShrinkM = ShrinkMotif( motif, map ); 
  for( int *ptr = map_1, *ptr_end = map_1 + _RealMotifSize; 
        ptr < ptr_end; ptr++) {
    *ptr = -1;
  }
  for( int i=0, i_end = ShrinkM.getNbrRow(); i < i_end; i++) {
    map_1[ map[i] ] = i;
  }

# ifdef MSG 
  motif.print("ShrinkM");
# endif

  // Motif
  int motif_NbrNodes = ShrinkM.getNbrRow();
  int *MapToG = new int[ motif_NbrNodes];
  int k_node = 0;

  // Get Topological Classes
  // ans constraints

  // Set shrunk motif constraints
  int *constrM = getUniqPath( motif );
  // keep the motif (not shrunk mot) contraints
  int *constr = new int[ ShrinkM.getNbrRow() ];
  for( int *ptr = constr, *ptr_map=map, *ptr_end = constr + ShrinkM.getNbrRow();
       ptr < ptr_end; ptr++, ptr_map++) {
    *ptr = constrM[ *ptr_map ];
    if ( *ptr != -1 ) {
      // Transform from Original to New indexes
      *ptr = map_1[ *ptr ];
    }
  }

  setConstraints( constr);

  printList( "Motif constraint ", _Constraints, &_Constraints[ ShrinkM.getNbrRow() ], true);

# ifdef MSG 
  printList( "Motif constraint ", _Constraints, &_Constraints[ ShrinkM.getNbrRow() ], true);
# endif

# ifdef MSG 
  motif.print("ShrinkM"); 
  G.print("Network"); 
# endif


  if ( ShrinkM.getSymmetry() &&  G.getSymmetry() ) {
    getTime( 0 );
    SparseAdjacency FShrinkM( ShrinkM );

    // Version Michel
 
    int n =  _RealMotifSize;
    int NbVertices = G_NbrNodes;
    int NbM = (n * (n - 1)) / 2;
    int NbMultipleVertices = 0;
    // Find the NbMultipleVertices
    int *NbAdj       = new int[ NbVertices ]; 
    int *SortedNbAdj = new int[ NbVertices ];
    for (int i = 0; i < NbVertices; i++)
      NbAdj[i] = 0;
    for (int i = 0; i < NbVertices; i++)
      SortedNbAdj[i] = G.getAllRowSize( i );

    sort( SortedNbAdj, SortedNbAdj + NbVertices);

    for (int i = 0; i < NbM; i++)
      NbMultipleVertices += SortedNbAdj[NbVertices - 1 - i];
    if (NbMultipleVertices > NbVertices)
      NbMultipleVertices = NbVertices;

    Neibourgh *MultipleVerticesTwoSources = new Neibourgh[NbMultipleVertices];

    for (int i = 0; i < NbMultipleVertices; i++)
      MultipleVerticesTwoSources[i].Sources = new int[2];
    Neibourgh *MultipleVerticesManySources = new Neibourgh[NbMultipleVertices];
    for (int i = 0; i < NbMultipleVertices; i++)
      MultipleVerticesManySources[i].Sources = new int[n];

     int *ReverseOccurrence = new int[NbVertices];
     for (int i = 0; i < NbVertices; i++)
       ReverseOccurrence[i] = -1;

    // end version Michel

    for ( int i = 0; i < G_NbrNodes; i++ ) {
      MapToG[0] = i;
      if( G.getAllRowSize( i ) != 0) {
	// Version Michel ??
	// kpp_MatchUndirectedMotifShrink( FShrinkM, G, k_node, 
	//			  MapToG
	//			);

	ReverseOccurrence[ i ] = 0;

	kpp_MatchUndirectedMotifShrink( FShrinkM, G, k_node, 
					MapToG,
					MultipleVerticesTwoSources,
					MultipleVerticesManySources,
					ReverseOccurrence
					);

        ReverseOccurrence[ i ] = -1;
      }
    }

    // Version Michel
    for (int i = 0; i < NbMultipleVertices; i++) {
	MultipleVerticesTwoSources[i].DestroyMe();
	MultipleVerticesManySources[i].DestroyMe();
    }

    delete[] MultipleVerticesTwoSources;
    delete[] MultipleVerticesManySources;
    delete[] ReverseOccurrence;
    delete [] NbAdj;
    delete [] SortedNbAdj;
    getTime( 1 );
 
  } else if ( ! (ShrinkM.getSymmetry() ||  G.getSymmetry() ) ) {
    getTime( 0 );
    SparseAdjacency FShrinkUpper( ShrinkM.getTriangular( true ) );
    SparseAdjacency FShrinkLower( ShrinkM.getTriangular( false ) );

# ifdef MSG 
    FShrinkUpper.print( "Upper");
    FShrinkLower.print( "Lower");
# endif    

    for ( int i = 0; i < G_NbrNodes; i++ ) {
      MapToG[0] = i;
      kpp_MatchDirectedMotifShrink(  FShrinkUpper,  FShrinkLower, G, k_node, 
    			       MapToG
			    );
    }
    getTime( 1 );
  } else {
#   ifdef MSG 
    cerr << "The Graph and the motif not the same characteristics "
	 << "(directed/undirected)" << endl;
#   endif
  }
  
  delete [] map;
  delete [] map_1;
  delete [] MapToG;

   vector< int*> *found_list=_FoundList;
  _FoundList = 0;
  return( found_list ); 
}

// Version Michel
void FindMotif::kpp_MatchUndirectedMotifShrink( 
			       SparseAdjacency &motif,  
			       SparseAdjacency &G, 
			       int k, 
			       int* MapToG,
			        Neibourgh *MultipleVerticesTwoSources,
			       Neibourgh *MultipleVerticesManySources,
			       int *ReverseOccurrence
			     ) {
  //
  // Undirected case
  //
    
  // ??? Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cout << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif


  if ( (k) == motif.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

#   ifdef MSG
    cout << "  New motif is stored " << endl;
#   endif

    if( _OutputMode == CountMode ) {
      int OccurenceSize = motif.getNbrRow();
      int *OccurrenceVertices = MapToG;
      /*
      int *NbNeibourghsInOccurrence = new int[OccurenceSize];
      for (int i = 0; i < OccurenceSize; i++)
	NbNeibourghsInOccurrence[i] = 0;
      */

      /*
      for (int N = 0; N < OccurenceSize; N++)
	for (int O = 0; O < OccurenceSize; O++)
	  if (BooleanAdjacencies[OccurrenceVertices[O]][OccurrenceVertices[N]])
	    NbNeibourghsInOccurrence[O]++;
      */
      // TODO : To optimize
      // Get the the number of neighbours inside the sub-graph
      // matching the motif, extracted from the graph G
      int *NbNeibourghsInOccurrence = getNumberOfNeighboursInSubSet( G, 
						MapToG, 
						&MapToG[OccurenceSize] );

      unsigned long int PR = ComputeNbOccurrences(G, OccurrenceVertices, OccurenceSize, &motif , NbNeibourghsInOccurrence, MultipleVerticesTwoSources, MultipleVerticesManySources, ReverseOccurrence );
  
      printList("NbNeibourghsInOccurrence : ", NbNeibourghsInOccurrence, 
		&NbNeibourghsInOccurrence[OccurenceSize], true );
    // NbResults += PR;
      /*
      if (PR > 0)
	{
	  int NS = OccurrenceVertices[0];
	  cout << NS << ", " << NbAdj[NS] << ", " << PR << endl;
	  if (PR != NbAdj[NS] * (NbAdj[NS] - 1))
	    {
	      cout << NS << ", " << NbAdj[NS] << endl;
	      cerr << NS << ", " << NbAdj[NS] << endl;
	    }
	}
      */

    /*std::cout << "Solution partielle : " << std::endl;
     std::cout << "Nombre de sommets : " << M->NbVertices << std::endl;
     for (int i = 0; i < M->NbVertices; i++)
     std::cout << OccurrenceVertices[i] + 1 << " ";
     std::cout << std::endl;
     std::cout << "Comptage partiel : " << PR << std::endl;*/
      delete[] NbNeibourghsInOccurrence;
#   ifdef MSG
      cout << "count  = " << _Count << ", adjust = " << PR << endl;
#   endif
      _Count += PR;
    } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
      // _Count += writeOccurences( G, motif, motif, MapToG );   
    }

  } else {


    /////////////////////////////////////////////////////////
    //
    // Get neighbours u of k+1 in the motif
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    const int *V_motif_kpp      = motif.getCol(k+1);
    const int *V_motif_kpp_end  = motif.getColEnd(k+1);
    int V_motif_kpp_size = V_motif_kpp_end - V_motif_kpp;

#   ifdef MSG
    cout << "  Motif-Neighbours  {u} = V_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", V_motif_kpp, V_motif_kpp_end, true );
#   endif

    //
    // Make the intersection :
    //    INTER{ V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ] ) }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
    //

    //
    // Intersection initialisation with V_G( MapToG (u=0) )
    //
    const int *it = V_motif_kpp;

    getTime( 2 );

    // Node ID in G
    // int  MapToG_u = MapToG[ *it ];

    // Ordered lists of neigbours of MapToG_u

    // Set to intersect
    const int **L = new const int*[V_motif_kpp_size];
    int       *TL = new int[V_motif_kpp_size];

    //
    // Store V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ]) set 
    // sorted heap.
    //
    int i = 0;
    for ( it = V_motif_kpp; (it != V_motif_kpp_end); it++) {
	L[i]    = G.getAllRow( MapToG[*it] );
	TL[i++] = G.getAllRowSize( MapToG[*it]);
    }
    Heap H(L, TL, V_motif_kpp_size );

#   ifdef MSG
    H.Debug();
#   endif

    //
    // Get the intersection list.
    //
    int a = H.NextCommonElement();
    while (a != NOT_INDEX) {

      // Apply constraints 
      if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	  || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {

	// if ( MapToG[ k+1 ] > MapToG[ constraint[ k+1 ] ]   ]
	if (!IsIn( MapToG, k+1, a)) {  

	  // New sub-motif found
	  MapToG[k+1] = a;
	  // Version Michel

	  ReverseOccurrence[a] = k;

	  kpp_MatchUndirectedMotifShrink( motif, G, k+1, MapToG, 
					  MultipleVerticesTwoSources,
					  MultipleVerticesManySources,
					  ReverseOccurrence
					  );
	  ReverseOccurrence[a] = -1;

	}
      }
      a = H.NextCommonElement();
    }
    delete[] L;
    delete[] TL;

    getTime( 3 );
  }
}

void FindMotif::kpp_MatchUndirectedMotifShrink( 
			       SparseAdjacency &motif,  
			       SparseAdjacency &G, 
			       int k, 
			       int* MapToG 
			     ) {
  //
  // Undirected case
  //
    
  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cout << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif


  if ( (k) == motif.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //

#   ifdef MSG
    cout << "  New motif is stored " << endl;
#   endif

    if( _OutputMode == CountMode ) {
      int64_t cumul = 1;
      const int64_t *G_nb_poor_neighbours = G.getNodeValues();
      /*
      const int64_t *nb_poor_neighbours = motif.getNodeValues();
      for (int i = 0; i <  motif.getNbrRow(); i++) {
	int64_t G_degree = G.getAllRowSize( MapToG[i] );
	int64_t M_degree = motif.getAllRowSize( i );
	int64_t Diff_degree = G_degree - M_degree;

	cout << "nb_poor_neighbours  :"   <<   nb_poor_neighbours[i] << endl;
	cout << "G_nb_poor_neighbours:" << G_nb_poor_neighbours[MapToG[i]] << endl;
	cout << " G_degree :" <<  G_degree << endl;
	cout << " M_degree :" <<  M_degree << endl;
	// if ( Diff_degree >= nb_poor_neighbours ) {
	cumul *= BinomCoef(Diff_degree, 
			   static_cast<int> (nb_poor_neighbours[i]));
	  // } else {
	  // Not enough motif neighbours
	  //cumul = 0;
	  //}
      }
*/     


      // TODO : Optimize
      
      unsigned long int choose = 1;
      unsigned long int fact = 1;

      for (int i = 0; i <  motif.getNbrRow(); i++) {
	int64_t G_degree = G.getAllRowSize( MapToG[i] );
	int64_t M_degree = motif.getAllRowSize( i );
	int64_t Diff_degree = G_nb_poor_neighbours[MapToG[i]]+ G_degree - M_degree;
	int64_t nb_poor_neighbours =  motif.getNodeValues()[i];
#   ifdef MSG	
	cout << "G_nb_poor_neighbours:" << G_nb_poor_neighbours[MapToG[i]] << endl;
	cout << "nb_poor_neighbours  :"   <<   nb_poor_neighbours << endl;

	cout << " G_degree :" <<  G_degree << endl;
	cout << " M_degree :" <<  M_degree << endl;
	cout << " Diff_degree :" <<  Diff_degree << endl;
#   endif	
	if ( Diff_degree >= nb_poor_neighbours ){
	  for (int j = 0 ; j < nb_poor_neighbours ; j++){
	    choose *= ( Diff_degree - j);
	    fact *= (j+1);
	  }
	} else { 
	  choose = 0;
	}
 
      }
      cumul= static_cast<int64_t> (double( choose)/double(fact));
#   ifdef MSG
      cout << " map : " ;
      for (int i = 0; i <  motif.getNbrRow(); i++) {
	cout << MapToG[i] << " " ;
      }
      cout << endl;

      cout << " cumul : " << cumul << endl;
#   endif

      _Count += cumul;
      // cout << "count ??? =" << Count << endl;
    } else if ( (_OutputMode == WriteMode) ||  (_OutputMode == StoreMode) ) {
      // _Count += writeOccurences( G, motif, motif, MapToG );   
    }

  } else {


    /////////////////////////////////////////////////////////
    //
    // Get neighbours u of k+1 in the motif
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    const int *V_motif_kpp      = motif.getCol(k+1);
    const int *V_motif_kpp_end  = motif.getColEnd(k+1);
    int V_motif_kpp_size = V_motif_kpp_end - V_motif_kpp;

#   ifdef MSG
    cout << "  Motif-Neighbours  {u} = V_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", V_motif_kpp, V_motif_kpp_end, true );
#   endif

    //
    // Make the intersection :
    //    INTER{ V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ] ) }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
    //

    //
    // Intersection initialisation with V_G( MapToG (u=0) )
    //
    const int *it = V_motif_kpp;

    getTime( 2 );

    // Node ID in G
    // int  MapToG_u = MapToG[ *it ];

    // Ordered lists of neigbours of MapToG_u

    // Set to intersect
    const int **L = new const int*[V_motif_kpp_size];
    int       *TL = new int[V_motif_kpp_size];

    //
    // Store V_G( MapToG[ u=V_motif(k+1) with (u<=k+1) ]) set 
    // sorted heap.
    //
    int i = 0;
    for ( it = V_motif_kpp; (it != V_motif_kpp_end); it++) {
	L[i]    = G.getAllRow( MapToG[*it] );
	TL[i++] = G.getAllRowSize( MapToG[*it]);
    }
    Heap H(L, TL, V_motif_kpp_size );

#   ifdef MSG
    H.Debug();
#   endif

    //
    // Get the intersection list.
    //
    int a = H.NextCommonElement();
    while (a != NOT_INDEX) {

      // Apply constraints 
      if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	  || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {

	// if ( MapToG[ k+1 ] > MapToG[ constraint[ k+1 ] ]   ]
	if (!IsIn( MapToG, k+1, a)) {  

	  // New sub-motif found
	  MapToG[k+1] = a;
	  kpp_MatchUndirectedMotifShrink( motif, G, k+1, MapToG );
	}
      }
      a = H.NextCommonElement();
    }
    delete[] L;
    delete[] TL;

    getTime( 3 );
  }
}


void FindMotif::kpp_MatchDirectedMotifShrink( 
			     SparseAdjacency &motifSup, 
			     SparseAdjacency &motifInf, 
			     SparseAdjacency &G, 
			     int k, 
			     int *MapToG
			   ) {
  //
  // Directed part
  //
    

  // Warning: the kth node must be connected at least once to 
  // the k-1 preceeding nodes

# ifdef MSG
  cout << "k-motif " << k << " -> " << MapToG[k] << endl;
# endif

  if ( (k) == motifSup.getNbrRow() - 1 ) {
    //
    // A new motif is found
    //
# ifdef MSG
    cout << "  New motif is stored " << endl;
# endif
    if( _OutputMode == CountMode ) {
      int64_t cumul = 1;   
      // Incomming
      const int64_t *nb_poor_neighbours_in  = motifInf.getNodeValues();
      for (int i = 0; i <  motifInf.getNbrRow(); i++) {
	int G_degree = G.getRowSize( MapToG[i] );
	int M_degree = motifInf.getRowSize( i );
	int Diff_degree = G_degree - M_degree;
	cumul *= BinomCoef(Diff_degree, 
			      static_cast<int> (nb_poor_neighbours_in[i*2]) );
      }
      if ( cumul != 0 ) {
	// Outcoming
	const int64_t *nb_poor_neighbours_out = motifSup.getNodeValues();
	for (int i = 0; i <  motifSup.getNbrCol(); i++) {
	  int G_degree = G.getColSize( MapToG[i] );
	  int M_degree = motifSup.getColSize( i );
	  int Diff_degree = G_degree - M_degree;
	  cumul *= BinomCoef( Diff_degree, 
			      static_cast<int> (nb_poor_neighbours_out[i*2+1]));
	}
      }
      
      _Count += cumul;
    } else  if( _OutputMode == WriteMode ) {
      // _Count += writeOccurences( G, motifInf, motifSup, MapToG );
    } 
  } else {
    

    int MapToG_u;

    /////////////////////////////////////////////////////////
    //
    // Get Incoming and Outcomming edges u: u -> k+1 or u <- k+1
    // with u <= k+1 node of the  pattern
    //
    /////////////////////////////////////////////////////////
    
    const int *Vin_motif_kpp      = motifInf.getRow(k+1);
    const int *Vin_motif_kpp_end  = motifInf.getRowEnd(k+1);
          int Vin_motif_kpp_size  = Vin_motif_kpp_end - Vin_motif_kpp;
    const int *Vout_motif_kpp     = motifSup.getCol(k+1);
    const int *Vout_motif_kpp_end = motifSup.getColEnd(k+1);
          int Vout_motif_kpp_size = Vout_motif_kpp_end - Vout_motif_kpp;

#   ifdef MSG
    cout << "  Motif-Incoming {u} = Vin_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vin_motif_kpp, Vin_motif_kpp_end, true );

    cout << "  Motif-Outcoming {u} = Vout_m(k+1), with u <= k+1 : "
	 << endl; 
    printList("      ", Vout_motif_kpp, Vout_motif_kpp_end, true );
#   endif
    
    // Make the intersection :
    //    INTER{ Vout_G( MapToG[ u=Vin_motif(k+1), with (u<=k+1) ] ) 
    //       and Vin_G( MapToG[ u=Vout_motif(k+1), with (u<=k+1) ] }
    // 
    // Remark : for Exact matching do something more
    // Possibles Nodes must have the same interaction pattern
      
    //
    // Intersection initialisation with Vout_G( MapToG[ u=0] )
    //
    const int *it;

    // Set to intersect
    const int **L = new 
      const int*[Vin_motif_kpp_size + Vout_motif_kpp_size];
    int       *TL = 
      new int[Vin_motif_kpp_size + Vout_motif_kpp_size];

    //
    // Store Vout_G( MapToG[ u=Vin_motif(k+1) with (u<=k+1) ]) set 
    // in the sorted heap.
    //
    int i = 0;
    for ( it = Vin_motif_kpp; (it != Vin_motif_kpp_end); it++) {
      MapToG_u =  MapToG[*it];
      TL[i] = G.getColSize( MapToG_u );
      if( TL[i] != 0) {
	L[i++]    = G.getCol( MapToG_u );
      } else {
	// Empty set
	i = 0;
	break;
      }
    }
    
    // Test if incoming arcs are required in the pattern (  Vin_motif_kpp_size )
    // and no outcomming arcs have found ( G.getCol( MapToG_u ) )
    if ( (i != 0) || ( Vin_motif_kpp_size == 0)) {
      //
      // Store Vin_G( MapToG[ u=Vout_motif(k+1) with (u<=k+1) ]) set 
      // in the sorted heap.
      //
      int j = i;
      for ( it = Vout_motif_kpp; (it != Vout_motif_kpp_end); it++) {
	MapToG_u =  MapToG[*it];
	TL[j] = G.getRowSize( MapToG_u );
	if( TL[j] != 0) {
	  L[j++]    = G.getRow( MapToG_u );
	} else {
	  // Empty set
	  j = i;
	}
      }

      // Build the Heap
      // TODO : Heap deal with no list (i == 0)  to avoid the test
      // Test if outcoming arcs are required in the pattern (  Vout_motif_kpp_size )
      // and no incomming arcs have found ( G.getRow( MapToG_u ) )
      if( (j > 0) && ((i != j) || ( Vout_motif_kpp_size == 0)) ) {
	Heap H(L, TL, j );

#     ifdef MSG
	H.Debug();
#     endif

	//
	// Get the intersection list.
	//
	int a = H.NextCommonElement();
	while (a != NOT_INDEX) {
	  
	  // Apply constraints 
	  if( ( _Constraints[ k+1 ] == NOT_INDEX ) 
	      || ( MapToG[ _Constraints[ k+1 ] ] < a ) ) {
	  
	    if (!IsIn( MapToG, k+1, a)) { 

	      // New sub-motif found
	
	      MapToG[k+1] = a;
	      kpp_MatchDirectedMotifShrink( motifSup, motifInf, G, k+1, 
				      MapToG );
	    }
	  }
	  a = H.NextCommonElement();
	}
      }
    }
    delete[] L;
    delete[] TL;

    
  }
}

void buildAndfindAllMotifShrink( SparseAdjacency &G, int k, FindMotif &find,
			   vector< int * > &motifs, vector<int64_t> &counts ) {

#   ifdef MSG
  cout << " Network nodes: " <<  G.getNbrRow() << endl;
  cout << " Network edges: " <<  G.getNbrOfValues() << endl;
#   endif

  int *map = new int[k];
  SparseAdjacency ShrinkG = ShrinkMotif( G, map );
 
#   ifdef MSG
  cout << " Shrink Network nodes: " <<  ShrinkG.getNbrRow() << endl;
  cout << " Shrink Network edges: " <<  ShrinkG.getNbrOfValues() << endl;
#   endif

  // Build Clique
  
  SparseAdjacency *motif = new SparseAdjacency( k, "clique", G.getSymmetry() );
  vector< SparseAdjacency *> motifList;

  // A pointer is requested
  motifList.push_back( motif );

  vector< int > CT;
  for( int i=0; i<k; i++){
    CT.push_back( i );
  }
  vector < vector< int > > CTmotif;
  CTmotif.push_back( CT );

  vector< vector< vector<int> >*> CTList;
  CTList.push_back( &CTmotif ) ;


  // Find occurences
  // int * map;
  SparseAdjacency ShrinkM = ShrinkMotif( *motif, map ); 


  find.findShrink( ShrinkG, ShrinkM );
  _CountAllMotifs = find.getCount();

# ifdef MSG
  cout << " -> Motif " << _NbrOfMotifs << " : " << find.getCount() << endl;
# endif

  _NbrOfMotifs = 1;

  // Store counts and motif
  counts.push_back( find.getCount() );
  motifs.push_back( motif->getIndexes( ) );


  // siblings_t *root = new siblings_t;
  // addMotif( root, root, *motif, 1);
  // root -> push_back( new MotifTree(0, 0 ) );
 
  // Init TC
  int nb_edges =  k*(k-1);
  if ( G.getSymmetry() )
    nb_edges = nb_edges / 2;

  findNextGeneration( ShrinkG, motifList, CTList, find, k, nb_edges, 
		      motifs, counts, ComputeMotifs, ShrinkMode );

# ifdef MSG
  cout << "Count all motifs " << _CountAllMotifs << endl;
  cout << "Nbr of motifs " << _NbrOfMotifs << endl;
  cout << "Nbr of motifs " << counts.size() << " " << motifs.size() << endl;
# endif
  
  /*
  const int64_t *counters = find.getCounters(); 
  for (int i=0; i < (*root)[0]->_id; i++) {
    cout << i << " " << counters[i] << endl;
  }
  */
  delete [] map; 
}

// The MapShrinkToOriginal storage must be allocate/delete by the calling block
// It must be have the size of the (complete) motif
SparseAdjacency ShrinkMotif( SparseAdjacency & motif, int* MapShrinkToOriginal)
{
  int NbrNodes =  motif.getNbrRow();

  int     NbrNewNodes = 0;
  int    *NewEdge;
  int64_t *NodeValue;
  size_t     NbrNodeValues;
  int *NewLabel = new int[NbrNodes];
  int *PoorVertex = new int[NbrNodes];

  if( motif.getSymmetry() ) {

    
    // Determine Vertices with only one neighbour (called PoorVertex)
     for (int i = 0; i < NbrNodes; i++)
      // Doute ??? normalement la matrice est symetrique
      if ( (motif.getRowSize(i) + motif.getColSize(i)) < 2)
	PoorVertex[i] = true;
      else
	PoorVertex[i] = false;
    
    // For each vertex, count "PoorNeighbour"s
    int *NbrPoorNeighbours = new int[NbrNodes];
    for (int i = 0; i < NbrNodes; i++)
      NbrPoorNeighbours[i] = 0;
    for (int i = 0; i < NbrNodes; i++){
      for (const int *it = motif.getRow(i); 
	   it !=  motif.getRowEnd(i);
	   it++) {
	if (PoorVertex[*it])
	  NbrPoorNeighbours[i]++;
      }
      // Doute ??? normalement la matrice est symetrique
      for ( const int *it = motif.getCol(i); 
	   it !=  motif.getColEnd(i);
	   it++) {
	if (PoorVertex[*it])
	  NbrPoorNeighbours[i]++;
      }
    }
    
    // Renumbering nodes
    
    for (int i = 0; i < NbrNodes; i++)
      if (!PoorVertex[i]) {
	MapShrinkToOriginal[NbrNewNodes] = i;
	NewLabel[i] = NbrNewNodes++;
      } else {
	NewLabel[i] = NOT_INDEX;
      }
    
    NewEdge = new int[(motif.getNbrOfValues() + 1) * 2];
    int edge_index = 0;
    NodeValue = new int64_t[NbrNewNodes];
    for (int i = 0; i < NbrNodes; i++) {
      if (!PoorVertex[i]) {
	for (const int *it = motif.getRow(i); 
	     it !=  motif.getRowEnd(i);
	     it++) {
	  if(  NewLabel[ *it ] != NOT_INDEX ) { 
	    NewEdge[ 2 * edge_index ] = NewLabel[ i ];
	    NewEdge[ 2 * edge_index + 1 ] = NewLabel[ *it ];
	    edge_index++;
	  }
	}
	NodeValue[ NewLabel[i] ] =  NbrPoorNeighbours[i];
      }
    }

    NewEdge[ 2 * edge_index ]     = NOT_INDEX;
    NewEdge[ 2 * edge_index + 1 ] = NOT_INDEX;

    NbrNodeValues = 1;

    delete [] NbrPoorNeighbours;

  } else {

    // Determine Vertices with only one neighbour (called PoorVertex)
    // For directed motif PoorVertex must have exactly one neighbour
    // (Incomming or Outcoming edge)
    for (int i = 0; i < NbrNodes; i++)
      if ( (motif.getRowSize(i) + motif.getColSize(i)) < 2)
	PoorVertex[i] = true;
      else
	PoorVertex[i] = false;
    
    // For each vertex, count "PoorNeighbour"s
    //  NbrPoorNeighbours[k][0] : incoming edges
    //  NbrPoorNeighbours[k][1] : outcoming edges
    int (*NbrPoorNeighbours)[2] = new int [NbrNodes][2];

    for (int i = 0; i < NbrNodes; i++) {
      NbrPoorNeighbours[i][0] = 0;
      NbrPoorNeighbours[i][1] = 0;
    }
    for (int i = 0; i < NbrNodes; i++){
      // Incoming arcs
      for (const int *it = motif.getRow(i); 
	   it !=  motif.getRowEnd(i);
	   it++) {
	if (PoorVertex[*it])
	  NbrPoorNeighbours[i][0]++;
      }
      // Outcoming arcs
      for (const int *it = motif.getCol(i); 
	   it !=  motif.getColEnd(i);
	   it++) {
	if (PoorVertex[*it])
	  NbrPoorNeighbours[i][1]++;
      }
    }
    
    // Renumbering nodes
    
    for (int i = 0; i < NbrNodes; i++)
      if (!PoorVertex[i]) { 
	MapShrinkToOriginal[NbrNewNodes] = i;
	NewLabel[i] = NbrNewNodes++;
      } else
	NewLabel[i] = NOT_INDEX;
    
    NewEdge = new int[(motif.getNbrOfValues() + 1) * 2];
    int edge_index = 0;
    // TODO : to change
    NodeValue = new int64_t[NbrNewNodes * 2];
    for (int i = 0; i < NbrNodes; i++) {
      if (!PoorVertex[i]) {
	// Outcoming arcs
	for (const int *it = motif.getCol(i); 
	     it !=  motif.getColEnd(i);
	     it++) {
	  if(  NewLabel[ *it ] != NOT_INDEX ) { 
	    NewEdge[ 2 * edge_index     ] = NewLabel[ i ];
	    NewEdge[ 2 * edge_index + 1 ] = NewLabel[ *it ];
	    edge_index++;
	  }
	}
	// Store Incoming, outcoming arcs.
	NodeValue[ NewLabel[i]*2    ] =  NbrPoorNeighbours[i][0];
	NodeValue[ NewLabel[i]*2 + 1] =  NbrPoorNeighbours[i][1];
      }
    }

    NewEdge[ 2 * edge_index ]     = NOT_INDEX;
    NewEdge[ 2 * edge_index + 1 ] = NOT_INDEX;
 
    NbrNodeValues = 2;
    delete [] NbrPoorNeighbours ;
  }

  SparseAdjacency ShrinkedMotif ( NewEdge, 
				  NbrNewNodes, 
				  motif.getSymmetry(), 
				  NodeValue, NbrNodeValues );  

  // ShrinkedMotif.print("shrinked");

  delete [] NewEdge;
  delete [] NodeValue;
  delete [] NewLabel;
  delete [] PoorVertex;

  return ShrinkedMotif ;
}

void getChildMotifs( vector< SparseAdjacency *> &AncestorMotifs, 
	       vector< vector< vector<int> > *> &AncestorCT, 
	       vector< SparseAdjacency *> &ChildMotifs,
		vector< vector< vector<int> > *> &ChildCT, 
		int motif_size, 
		int nb_edges,
		bool symmetry_, 
		FindMotif &find) {

  int *perm = new int[ motif_size ];
  int k=0;

  // For all parent motifs
  for( vector< SparseAdjacency *>::iterator it=AncestorMotifs.begin() ;
       it != AncestorMotifs.end(); it++, k++) {
    vector< vector<int> > *CT = AncestorCT[k];
    
#   ifdef MSG
    cerr << "Topo Classes of motif " << k 
	 << ", number of classes " << CT-> size()<< endl;
    // for( int i = 0; i < CT->size(); i++) {
    for( vector< vector<int> >::iterator it_tc = CT->begin();
	 it_tc != CT->end(); it_tc++ ) {
      cerr << "  ";
      // for( int j = 0; j < (*CT)[i].size(); j++) {
      for( vector<int>::iterator it = it_tc->begin(); 
	   it != it_tc->end(); it++) {
	cerr << *it << ", ";
      }
      cerr << endl;
    }
#   endif

#   ifdef MSG
    cerr << "GENERATION " << nb_edges - 1 
	 << " ancetre " << k << " number of CT " << CT->size() << endl;
#   endif
    // cout << "GENERATION " << nb_edges - 1 
    //	 << " ancetre " << k << " number of CT " << CT->size() << endl;

    if( nb_edges >= (motif_size-1) ) {
      //
      // Select to nodes in all CT configurations 
      //
      for( size_t i = 0; i < CT->size(); i++) {
	for( size_t j = 0; j <= i; j++) {
	  int src = 0;
	  int dst = 0;
	  
#         ifdef MSG
	  cerr << " CT id : " << i << ", " << j << endl;
#         endif

	  // Same TC case
	  if( i == j) {
	    // Need two nodes in the same class
	    if ( (*CT)[i].size() > 1 ) {
	      dst = 1;
	    } else
	      dst = NOT_INDEX;
	  }

	  if (dst != NOT_INDEX) {
	    
	    // Two arcs to test
	    int u_max = 2;
	    if ( symmetry_ || ( i==j ) )
	      u_max = 1;

	    // Selected nodes
	    int src_node = (*CT)[i][src];
	    int dst_node = (*CT)[j][dst];

	    // Two arcs to remove
	    for ( int u=0; u < u_max ; u++) {

	      // Change src and dst
	      if ( u == 1 ) {
		// Selected nodes
		int tmp = src_node;
 		src_node = dst_node;
		dst_node = tmp;
	      }

	      SparseAdjacency *motif = new SparseAdjacency( *AncestorMotifs[k] );
#             ifdef MSG	      	      
	      cerr << " remove edge : " << src_node << ", " << dst_node << endl;
#             endif

	      // Remove edge/arc
	      //   If no edge/arc is present ( returned 
	      //   by removedArc) then
	      //   it the same for all nodes in the same CTs
	      if ( motif -> removeArc( src_node, dst_node ) 
		   && motif -> optimizeConnexity( perm ) ) {
		// printList( "permutation", 
		//	 perm, &perm[motif_size], true);
		
		bool no_isomorph =  true;
		// Remove the motif if an isomorphim exists with one 
		// in the ChildMotifs list
		for( vector< SparseAdjacency *>::iterator jt=ChildMotifs.begin() ;
		     jt != ChildMotifs.end() && no_isomorph; jt++) {
		  // (*jt)->print( "test iso", "all");
		  no_isomorph = ! find.isNonExactIsomorph( *motif, *(*jt) );
		}
		if( no_isomorph ) {
		  // Store the new motif
		  motif->print( "New motif ", "all" );
		  vector< vector<int> > *CT = find.getTopologicClasses( *motif );

		  ChildCT.push_back( CT );
		  ChildMotifs.push_back( motif );
		} else {
		  // motif->print( "delete isomorphe ", "all" );
		  delete motif;
		}
	      } else {
		// motif->print( "delete optimize/ nor edge ", "all" );
		delete motif;
	      }
	    }
	  }
	}
      }
    }
  }
  delete [] perm;
}

void setMotifDirectory( const char * motif_dir) {
  MotifDir = motif_dir;
}

# endif
