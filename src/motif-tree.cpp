# include <vector>
# include "motif-tree.h"


void printTList( const char *ptr, 
		const int *begin, const int *end, 
		bool endl_ ) {

# ifdef MSG 
  cout << ptr;
  for( const int *it = begin; (it != end ); it++) {
      cout << " " << *it;
  }
  if( endl_ )
    cout << endl;
# endif
}


void printTree(siblings_t *current_sibs, string &indent, int motif_size ) {

  indent.append( "  " );

# ifdef MSG 

  for ( siblings_t::iterator it = current_sibs->begin(); 
	it != current_sibs->end(); it++) {
    cout << indent ;
    printTList(" Col", (*it)->_col, (*it)->_col_end, true);
    cout << indent << " child " << ((*it)->_child) << endl; 
    cout << indent << " ID " << (*it)->_id << endl; 
    cout << indent << " Constraints " << endl; 
    for (vector<int *>::iterator cit=(*it)->_constraint.begin(); 
	 cit != (*it)->_constraint.end(); cit++) {
      if ( (*cit) ) {
	cout << indent << " -> " ; 
	for( int l = 0; l < motif_size; l++) {
	  cout << (*cit)[l] << " ";
	}
	cout << endl; 
      } 
    }
    cout << indent << " Flux In" << endl; 
    cout << indent << " -> " ; 
    for (vector<int >::iterator fit=(*it)->_flux_in.begin(); 
	 fit != (*it)->_flux_in.end(); fit++) {
	  cout << (*fit) << " ";
    }
    cout <<endl;  ; 

    cout << indent << " Flux Out" << endl; 
    cout << indent << " -> " ; 
    for (vector<int >::iterator fit=(*it)->_flux_out.begin(); 
	 fit != (*it)->_flux_out.end(); fit++) {
	  cout << (*fit) << " ";
    }
    cout <<endl;  ; 
      
    if( (*it)->_child ) {
      printTree( (*it)->_child, indent, motif_size );
    }  
  }
# endif
  indent.erase( indent.end() - 2 , indent.end() );
}


void addMotif( siblings_t *root, siblings_t *current_sibs, 
	       SparseAdjacency &motif, int k ) {
  

  if( motif.getSymmetry() ) {

    bool not_found = true;
    for ( siblings_t::iterator it = current_sibs->begin(); 
	  it != current_sibs->end(); it++) {
      if( ( (*it)->_col_size == motif.getColSize(k) ) 
	  && equal( (*it) -> _col,  (*it) -> _col_end, motif.getCol(k) ) 
	  ) {
	// Next generation
	not_found = false;
	addMotif( root, (*it)-> _child, motif, k+1 );
	// if ( k != (motif.getNbrRow()-1) )
	(*it)-> _id++;
      }  
      // Next brother/sister
    }

    // New motif
    // Add a new brother/sister
    if ( not_found ) {
      // Number of motifs
      
      // Compute constraints
      int *constraints = 0;
      siblings_t *child = 0;
      int id ;
      if ( root -> size( ) != 0 ) {
	id = (*root)[0] -> _id;
      } else {
	id = 0;
      }

      int counter;
      for (int m = motif.getNbrRow()-1; m > k; m--) {
	siblings_t *siblings = new siblings_t;
	if ( m ==  motif.getNbrRow()-1 ) {
	  counter = id;
	} else {
	  counter = 1;
	}
	// Constraint
	siblings -> push_back( 
			      new MotifTree( motif.getCol(m), 
					     motif.getColSize(m), 
					     constraints, constraints[m], 
					     child, counter ) 
			       );
	child = siblings;
      }
      if ( k ==  motif.getNbrRow()-1 ) {
	counter = id;
      } else {
	counter = 1;
      }
      
      current_sibs -> push_back( 
				new MotifTree( motif.getCol(k), 
					       motif.getColSize(k), 
					       constraints, constraints[k], 
					       child, counter ) 
				 );
      
      string str("new generate");
      printTree( current_sibs, str, motif.getNbrCol() );

      // if ( k ==  (motif.getNbrRow() - 1))
      //	(*root)[0] -> _id++;
    }



  } else {
# ifdef MSG 
    cout << "Not implemented" << endl;
# endif
  }

}

void addMotifConstraints( siblings_t *root, siblings_t *current_sibs, 
			  SparseAdjacency &motif, int *constraints, int k ) {
  

  if( motif.getSymmetry() ) {
    bool not_found = true;
    for ( siblings_t::iterator it = current_sibs->begin(); 
	  it != current_sibs->end(); it++) {

      if( ( (*it)->_col_size == motif.getColSize(k) ) 
	  && equal( (*it) -> _col,  (*it) -> _col_end, motif.getCol(k)) 
	  ) {

	(*it) -> _constraint.push_back( constraints );
	// Next generation
	not_found = false;
	addMotifConstraints( root, (*it)-> _child, motif, constraints, k+1 );
	// if ( k != (motif.getNbrRow()-1) )
	(*it)-> _id++;
      }  
      // Next brother/sister
    }

    // New motif
    // Add a new brother/sister
    if ( not_found ) {
      // Number of motifs
      
      siblings_t *child = 0;
      int id ;
      if ( root -> size( ) != 0 ) {
	id = (*root)[0] -> _id;
      } else {
	id = 0;
      }

      int counter;
      for (int m = motif.getNbrRow()-1; m > k; m--) {
	siblings_t *siblings = new siblings_t;
	if ( m ==  motif.getNbrRow()-1 ) {
	  counter = id;
	} else {
	  counter = 1;
	}
	// Constraint
	siblings -> push_back( 
			      new MotifTree( motif.getCol(m), 
					     motif.getColSize(m), 
					     constraints, constraints[m], 
					     child, counter ) 
			       );
	child = siblings;
      }
      if ( k ==  motif.getNbrRow()-1 ) {
	counter = id;
      } else {
	counter = 1;
      }
      
      current_sibs -> push_back( 
				new MotifTree( motif.getCol(k), 
					       motif.getColSize(k), 
					       constraints, constraints[k], 
					       child, counter ) );
      
      string str("new generate");
      printTree( current_sibs, str,  motif.getNbrRow() );
      
      //if ( k ==  (motif.getNbrRow() - 1))
      //	(*root)[0] -> _id++;
    }

  } else {
# ifdef MSG 
    cout << "Not implemented" << endl;
# endif
  }
  
}

void optimizeConstraints( siblings_t *current_sibs, int k ) {

  /// if( motif.getSymmetry() ) {
  for ( siblings_t::iterator it = current_sibs->begin(); 
	it != current_sibs->end(); it++) {
    bool no_constraint = false;
    // Flux in
    for (vector<int *>::iterator cit=(*it)->_constraint.begin(); 
	 cit != (*it)->_constraint.end(); cit++) {
      if ( (*cit)[k-1] == NOT_INDEX ) {
	no_constraint = true;
      }
      (*it)->_flux_in.push_back((*cit)[k-1]);
    }
    if ( no_constraint ) {
      ((*it)->_flux_in)[0] = NOT_INDEX;
      (*it)->_flux_in.resize(1);
    } else {
      sort( (*it)->_flux_in.begin(), (*it)->_flux_in.end() );
      vector<int >::iterator 
	fit = 
	unique( (*it)->_flux_in.begin(), (*it)->_flux_in.end() );
      (*it)->_flux_in.resize( fit - (*it)->_flux_in.begin() );
    }
    // Copy
    (*it)->_flux_in_t = new int[ (*it)->_flux_in.size() ];
    (*it)->_flux_in_end_t = (*it)->_flux_in_t + (*it)->_flux_in.size() ;
    
    int l=0;
    for( vector<int >::iterator fit = (*it)->_flux_in.begin(); 
	 fit != ((*it)->_flux_in.end());  fit++, l++ ){
      (*it)->_flux_in_t[l] = *fit;
    }
    
    
    // Flux out
    for (vector<int *>::iterator cit=(*it)->_constraint.begin(); 
	 cit != (*it)->_constraint.end(); cit++) {
      //	if ( (*cit)[k] == NOT_INDEX ) {
      //  no_constraint = true;
      //}
      (*it)->_flux_out.push_back((*cit)[k]);
    }
    sort( (*it)->_flux_out.begin(), (*it)->_flux_out.end() );
    vector<int >::iterator cit = 
      unique( (*it)->_flux_out.begin(), (*it)->_flux_out.end() );
    (*it)->_flux_out.resize( cit - (*it)->_flux_out.begin() );
    
    
    // Copy
    (*it)->_flux_out_t = new int[ (*it)->_flux_out.size() ];
    (*it)->_flux_out_end_t = (*it)->_flux_out_t + (*it)->_flux_out.size() ;
    
    l=0;
    for( vector<int >::iterator fit = (*it)->_flux_out.begin(); 
	 fit != ((*it)->_flux_out.end());  fit++, l++ ){
      (*it)->_flux_out_t[l] = *fit;
    }
    
    if( (*it)->_child ) 
      optimizeConstraints( (*it)->_child, k+1 );
    
  }
}
      

