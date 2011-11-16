//  sparse-adjacency.h
// 
//  Copyright (C) 2009 Laboratoire Statistique & Genome
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or (at
//  your option) any later version.
// 
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
//

//! \file    sparse-adjacency.h
//! \author  Gilles Grasseau
//! \brief   Using adjacency matrix
//! \version 0.1
//! \date    September 2009


# ifndef _FS_ADJACENCY_H_
# define _FS_ADJACENCY_H_

# include <iostream>
# include <fstream>
# include <vector>
# include <algorithm>
# include <string.h>
# include <stdint.h>

# include "util.h"
// # include "common-matrix.h"

using namespace std;

# define NOT_INDEX -1
# define OPT_NOCHECK_BOUNDS 0

// TODO : 
//  + Complete the method tests
//  + End the comments
//  + Reading from file method

//! \brief Permute array row and column containers and set
//!        the egde list
//! \param sigma permutation
//! \param array NNodes size array to permute
//! \param inc group of values to permute together

template <typename T> 
void _permute( const int *sigma, T *array, int array_size, int inc=1);

//! \class ??? SparseAdjacency 
//! \brief ???
//! \brief The SparseAdjacency class deals with using
//!  binary like adjacency matrix. SparseAdjacencyMatrix use a sparse 
//!  squared matrix representation.


class SparseAdjacency {

 protected :

  int _MatrixSize;

  //! Specify if a symetric representation is used.
  //! In this case only the upper-triangular matrix is stored (i <= j).
  //! remark : if this flag is false a SparseAdjacencyMatrix can be 
  //! symetric 
  bool symmetric_;

  // Node names
  vector<const char *> NodeNames;

  // Values related to nodes
  // ??? type
  int64_t *NodeValues_;
  int NodeValuesInc_;

  int64_t nbr_values_;
  int64_t nbr_diagonal_values_;

  int **row_;
  int **col_;

  int *row_size;
  int *col_size;

  // End pointer of the arrays
  int **row_end;
  int **col_end;

  // To have all the indexes of the row for a symmetric matrix
  int **all_row_;
  int **all_row_end;
  int *all_row_size;

  //! \brief Remove the (i,j) matrix element in row arrays
  //! \param i Row index.
  //! \param val value to remove.
  bool inline _removeInRow( int i, int val ) {
    
    int *it = remove( row_[i], row_end[i], val );
    bool found = false;

    if ( it != row_end[i] ) {
      row_size[i]--;
      row_end[i] = it; 
      found = true;
      it = remove( all_row_[i], all_row_end[i], val );
      all_row_size[i]--;
      all_row_end[i] = it; 
      
      it = remove( all_row_[val], all_row_end[val], i );
      all_row_size[val]--;
      all_row_end[val] = it; 
    }
    return found;
  }

  //! \brief Remove the (i,j) matrix element in colunn arrays
  //! \param i row index.
  //! \param val value to remove.
  bool inline _removeInCol( int i, int val ) {
    
    int *it = remove( col_[i], col_end[i], val );
    bool found = false;

    if ( it != col_end[i] ) {
      col_size[i]--;
      col_end[i] = it; 
      found = true;
    }
    return found;
  }

  //! \brief Allocate row and column descriptors
  //! \param size size of the graph
  void _allocateRowAndColumnDescriptors (int mat_size ) {
 
    _MatrixSize = mat_size;
    row_     = new int*[_MatrixSize]; 
    row_size = new int [_MatrixSize];
    row_end  = new int*[_MatrixSize];
    int **it1 = row_, **it2 = row_end;
    for ( int *it = row_size, *it_end = &row_size[_MatrixSize]; 
	  it != it_end; it++, it1++, it2++) {
      *it = 0;
      *it1 = 0; *it2 = 0;
    }  

    all_row_     = new int*[_MatrixSize]; 
    all_row_size = new int [_MatrixSize];
    all_row_end  = new int*[_MatrixSize];
    it1 = all_row_; it2 = all_row_end;
    for ( int *it = all_row_size, *it_end = &all_row_size[_MatrixSize]; 
	  it != it_end; it++, it1++, it2++) {
      *it  = 0;
      *it1 = 0; *it2 = 0;
    }  

    col_     = new int*[_MatrixSize];
    col_size = new int [_MatrixSize];
    col_end  = new int*[_MatrixSize];

    it1 = col_; it2 = col_end;
    for ( int *it = col_size, *it_end = &col_size[_MatrixSize]; 
	  it != it_end; it++, it1++, it2++ ) {
      *it = 0;
      *it1 = 0; *it2 = 0;
    }  
  }

  //! \brief Fill the Sparse matrix with a dense adjacency Matrix.
  //! ??? Attention les boucles sont mal gérées pour le non-dirigé : 
  //! elles sont présentes à la fois dans row_ et dans col_. 
  //! Il faut voir sir l'algo de bactracking a besoin de ce doublon
  //! il y a aussi un pb avec la definition i <= j upper et 
  //! i > j lower 
  void _fillDenseMatrix( const bool mat[], const int  matrix_size, 
			 const bool symmetric = false ) {
    
    if( symmetric_ ) {
      for( int i=0; i < _MatrixSize; i++) {
	int sum_row = 0; 
	int sum_col = 0; 
	const bool *line = &mat[i*_MatrixSize];
	
	// Lower part
	for( const bool *ptr=&line[0], *ptr_end=&line[MIN(i+1,_MatrixSize)]; 
	     ptr != ptr_end; ptr++) {
	  if( *ptr ) sum_col++;
	}

	// Upper part
	for( const bool *ptr=&line[i], *ptr_end=&line[_MatrixSize]; 
	     ptr != ptr_end; ptr++) {
	  if( *ptr ) sum_row++;
	}
	row_size[i] = sum_row;
	row_[i]     = new int[ sum_row ];
	row_end[i]  = row_[i] + sum_row;
	
	col_size[i] = sum_col;
	col_[i]     = new int[ sum_col ];
	col_end[i]  = col_[i] + sum_col;
	
	all_row_size[i] = sum_row + sum_col;
	all_row_[i]     = new int[ sum_row + sum_col ];
	all_row_end[i]  = all_row_[i] + sum_row + sum_col;
	
	// Lower triangular part
	int *store     = col_[i];
	int *all_store = all_row_[i]; 
	int index = 0;
	for( const bool *ptr=&line[0], *ptr_end=&line[MIN(i+1,_MatrixSize)]; 
	     ptr != ptr_end; ptr++, index++) {
	  if( *ptr ) {
	    *store++     = index;
	    *all_store++ = index;
	  }
	}

	// Upper triangular part
	store     = row_[i];
	// Redo diagonal
	index--;
	for( const bool *ptr=&line[i], *ptr_end=&line[_MatrixSize]; 
	     ptr != ptr_end; ptr++, index++) {
	  if( *ptr ) {
	    *store++     = index;
	    *all_store++ = index;
	  }
	}
      }
    } else {
      for( int i=0; i < _MatrixSize; i++) {
	// ???
	int sum_row = 0; 
	int sum_col = 0;
	const bool *line = &mat[i];
	const bool *column = &mat[i*_MatrixSize];
	for( const bool *ptr=&line[0], *ptr_end=line+_MatrixSize*_MatrixSize; 
	     ptr != ptr_end; ptr+=_MatrixSize) {
	  if( *ptr ) sum_row++;
	}
	for( const bool *ptr=&column[0], *ptr_end=&column[_MatrixSize]; 
	     ptr != ptr_end; ptr ++) {
	  if( *ptr ) sum_col++;
	}
	row_size[i] = sum_row;
	row_[i]     = new int[ sum_row ];
	row_end[i]  = row_[i] + sum_row;

	col_size[i] = sum_col;
	col_[i]     = new int[ sum_col ];
	col_end[i]  = col_[i] + sum_col;

	int *store     = row_[i];
	int index = 0;
	for( const bool *ptr=&line[0], *ptr_end=line+_MatrixSize*_MatrixSize; 
	     ptr != ptr_end; ptr+=_MatrixSize, index++) {
	  if( *ptr ) {
	    *store++     = index;
	  }
	}

	store  = col_[i];
	index  = 0;
	for( const bool *ptr=&column[0], 
	       *ptr_end=&column[_MatrixSize]; 
	     ptr != ptr_end; ptr++, index++) {
	  if( *ptr ) {
	    *store++     = index;
	  }
	}
      }
    }
    
  }

  //! \brief Delete row and column descriptors
  void _deleteRowAndColumnDescriptors () {
 
    if( row_ ) {
      delete [] row_;
      row_ = 0;
    }
    if( row_end ) {
      delete [] row_end;
      row_end = 0;
    }
    if( row_size ) {
      delete [] row_size;
      row_size = 0;
    }
    // Inv    if ( symmetric_ ) {
      if( all_row_ ) {
	delete [] all_row_;
	all_row_ = 0;
      }
      if( all_row_end ) {
	delete [] all_row_end;
	all_row_end = 0;
      }
      if( all_row_size ) {
	delete [] all_row_size;
	all_row_size = 0;
      }
      // Inv }

    if( col_ ) {
      delete [] col_;
      col_ = 0;
    }
    if( col_end ) {
      delete [] col_end;
      col_end = 0;
    }
    if( col_size ) {
      delete [] col_size;
      col_size = 0;
    }
  }

  //! \brief Delete row and column containers 
  void _deleteRowAndColumnContents( ) {


    for( int **it = row_, **it_end = &row_[_MatrixSize];
	 it != it_end; it++) {
      delete [] *it;
    }
    for( int **it = col_, **it_end = &col_[_MatrixSize];
	 it != it_end; it++) {
      delete [] *it;
    }
    if ( symmetric_ ) {
      for( int **it = all_row_, **it_end = &all_row_[_MatrixSize];
	   it != it_end; it++) {
	delete [] *it;
      }
    }
  }

  void _deleteNodeValues() {
    if( NodeValues_) {
      delete [] NodeValues_;
      NodeValues_ = 0;
    }
  }

  //! \brief Allocate row and column containers and set
  //! the egde list
  //! \param size size of the graph
  void _copyRowAndColumnContents( const int *ij_index ) {

    int i_ind, j_ind;

    // ISO C++
    // vector<int> row_v[_MatrixSize];
    // vector<int> col_v[_MatrixSize];
    vector<int> *row_v = new vector<int>[_MatrixSize];
    vector<int> *col_v = new vector<int>[_MatrixSize];

    if( symmetric_ ) {
      for( int i=0; ij_index[2*i] != NOT_INDEX; i++) {

	// Transform to an upper triangular matrix
	//
	if( ij_index[2*i] <= ij_index[2*i+1] ) {
	  i_ind = ij_index[2*i];
	  j_ind = ij_index[2*i+1];
	} else {
	  i_ind = ij_index[2*i+1];
	  j_ind = ij_index[2*i];
	}
	row_v[ i_ind ].push_back( j_ind );
	col_v[ j_ind ].push_back( i_ind );
      }
    } else {
      for( int i=0; ij_index[2*i] != NOT_INDEX; i++) {
	row_v[ ij_index[2*i+1]   ].push_back( ij_index[2*i] );
	col_v[ ij_index[2*i]     ].push_back( ij_index[2*i+1] );
      }
    }

    // Sort and uniq rows
    vector<int>::iterator it_;
    for( int i=0; i < _MatrixSize; i++) {
      sort( row_v[i].begin(), row_v[i].end() );
      it_ = unique( row_v[i].begin(), row_v[i].end() );
      row_v[i].resize( it_ - row_v[i].begin());
    }

    // Sort and uniq columns
    for( int j=0; j < _MatrixSize; j++) {
      sort( col_v[j].begin(), col_v[j].end() );
      it_ = unique( col_v[j].begin(), col_v[j].end() );
      col_v[j].resize( it_ -  col_v[j].begin());
    }

    for ( int i = 0; i < _MatrixSize; i++) {
      row_size[i] = row_v[i].size();
      if ( row_size[i] ) {
	row_[i]    = new int[ row_size[i] ]; 
	row_end[i] = row_[i] + row_size[i];
	int *it_dest = row_[i];
	for( vector<int>::iterator it = row_v[i].begin(), it_end = row_v[i].end();
	       it != it_end; it++, it_dest++) {
	  *it_dest = * it;
	}
      } else {
	row_[i]    = 0;
	row_end[i] = 0;
      }
    }
    for ( int j = 0; j < _MatrixSize; j++) {
      col_size[j] = col_v[j].size();
      if ( col_size[j] ) {
	col_[j]    = new int[ col_size[j] ]; 
	col_end[j] = col_[j] + col_size[j];
	int *it_dest = col_[j];
	for( vector<int>::iterator it = col_v[j].begin(), 
	       it_end = col_v[j].end();
	       it != it_end; it++, it_dest++) {
	  *it_dest = * it;
	}
      } else {
	col_[j]    = 0;
	col_end[j] = 0;
      }
    }

    // Update all_row
    if( symmetric_ ) {
      for( int i=0; i < _MatrixSize; i++) {
	all_row_size[i] = row_size[i] + col_size[i] ;
	all_row_[i]     = new int [ all_row_size[i]  ];
	all_row_end[i]  = all_row_[i] + all_row_size[i] ;
	
	int *ptr = all_row_[i] ;
	int *it;
	for ( it = col_[i]; 
	      it != col_end[i]; it++ ) {
	  *ptr = *it; ptr++;
	}
	for ( it = row_[i]; 
	      it != row_end[i]; it++ ) {
	  *ptr = *it; ptr++;
	}
      }
    }

    delete [] row_v;
    delete [] col_v;
  }

  //! \brief Copy the array in NodeValues field
  //! \param NodeValues values to copy. The array size must be equal 
  //!    to  NodeValuesInc * MatrixSize_.
  //! \param NodeValuesInc number of values affected to one node
  void setNodeValues( const int64_t *NodeValues = 0,
		      int NodeValuesInc = 1 ) { 

    NodeValuesInc_ = NodeValuesInc;
    if ( NodeValues != 0 ) {
      int array_size = _MatrixSize;
      NodeValues_ = new int64_t[ array_size * NodeValuesInc_];
      const int64_t *src = NodeValues;
      for( int64_t *it = NodeValues_; 
	   it <  (NodeValues_ + array_size * NodeValuesInc_); ) 
	*it++ =  *src++;
    } else {
      NodeValues_ = 0;
    }
  }

  //! \brief Check if the pair of indexes contained in 'array' are 
  //!        consistent with matrix dimensions.
  //!        Generate an error if one index is out of bounds.
  //! \param array indexes to check.
  //! \return Return the numer of pair indexes. 
  size_t inline _checkEdgeBounds( const int *array ) {
    bool outofbounds = false;
    int k = 0;

    for( int i=0;  array[i] != NOT_INDEX; ) {
      k++;
      if ( (array[i] >= _MatrixSize ) || (array[i+1] >= _MatrixSize) )
	outofbounds = true;
      i++;
      i++;
    }
    if ( outofbounds ) {
#   ifdef MSG 
      LOGMSG( 1, std::cout, "Index greater than one matrix dimension", "" );
#   endif
    }
    return k;
  }


  //! \brief Method filling the adjacency matrix with
  //!        the edje list (i_index, j_index) (same as the constructor)
  //! \param ij_index  Array of index pairs (i,j) (array[2, ...], i and j >= 0).
  //!                  This array must be ended with the values (NOT_INDEX, NOT_INDEX)
  //! \param symmetric Reduce the matrix to a upper-triangular matrix (i <= j).
  //!                 If necessary a matrix element (if i>j) could be transposed
  //!                 to obtain an upper-triangular form.
  // TODO 
  //             else "ptr" is directly assign to the internal matrix data.
  //             In "false" case do not use the pointer for other 
  //              computations 
  void  initMatrix( const int ij_index[], 
		    const int matrix_size=0, 
		    const bool symmetric = false, 
		    const int64_t *NodeValues = 0,
		    size_t NodeValuesInc = 1 ) { 


    bool outofbounds = false;

    symmetric_ = symmetric;

    _MatrixSize = matrix_size ;

    //
    // Verify indexes
    //
    _checkEdgeBounds( ij_index );
 
    if ( outofbounds ) { 
#   ifdef MSG   
      LOGMSG( 1, std::cout, "Index greater than one matrix dimension", "" );
#   endif
    }

    _allocateRowAndColumnDescriptors( _MatrixSize );

    _copyRowAndColumnContents( ij_index);

    updateNbrOfValues_();   

    symmetric_ = symmetric;

    setNodeValues( NodeValues, NodeValuesInc);

  }

  // \brief Copy method
  // \param SparseAdjacency matrix
  void 
  copyMatrix_( const SparseAdjacency &A) {

    _MatrixSize = A._MatrixSize; 
    symmetric_ = A.symmetric_;

    row_     = new int* [_MatrixSize];
    col_     = new int* [_MatrixSize];
    row_size = new int  [_MatrixSize];
    col_size = new int  [_MatrixSize];
    row_end  = new int* [_MatrixSize];
    col_end  = new int* [_MatrixSize];
    all_row_     = new int* [_MatrixSize];
    all_row_end  = new int* [_MatrixSize];
    all_row_size = new int [_MatrixSize];
  
    for( int i=0; i < _MatrixSize; i++) {
      row_[i]     = new int [ A.row_size[i] ];
      row_size[i] = A.row_size[i] ;
      row_end[i]  = row_[i] + row_size[i];

      int *ptr = row_[i] ;
      for ( int* it = A.row_[i]; it != A.row_end[i]; it++) {
	*ptr = *it; ptr++;
      }

      if( ! symmetric_ ) {
	// all_row not used for unsymmetic matrix
	all_row_[i] = 0;
	all_row_end[i] = 0;
	all_row_size[i] = 0;

      } else {
	all_row_[i]     = new int [ A.all_row_size[i] ];
	all_row_size[i] = A.all_row_size[i] ;
	all_row_end[i]  = all_row_[i] + all_row_size[i];
	ptr = all_row_[i] ;
	for ( int* it = A.all_row_[i]; it != A.all_row_end[i]; it++) {
	  *ptr = *it; ptr++;
	}
      }
    }

    for( int i=0; i < _MatrixSize; i++) {
      col_[i]     = new int [  A.col_size[i] ];
      col_size[i] = A.col_size[i];
      col_end[i]  = col_[i] + col_size[i];

      int *ptr = col_[i] ;
      for (  int* it = A.col_[i]; it != A.col_end[i]; it++) {
	*ptr = *it; ptr++;
      }
    }
    symmetric_           = A.symmetric_;
    nbr_values_          = A.nbr_values_;
    nbr_diagonal_values_ = A.nbr_diagonal_values_; 

    int u = A.getNodeValuesInc();

    const int64_t *v = A.getNodeValues();

    setNodeValues( v, u );

    NodeNames = A.NodeNames ;
  } 

  // \brief Copy-submatrix method
  // \param SparseAdjacency matrix
  // Values in the matrix are not copied
  void 
    copyMatrix_( const SparseAdjacency &A, int matrix_size) {
    _MatrixSize =  matrix_size;
    symmetric_ = A.symmetric_;

    row_     = new int* [_MatrixSize];
    col_     = new int* [_MatrixSize];
    row_size = new int  [_MatrixSize];
    col_size = new int  [_MatrixSize];
    row_end  = new int* [_MatrixSize];
    col_end  = new int* [_MatrixSize];
    all_row_     = new int* [_MatrixSize];
    all_row_end  = new int* [_MatrixSize];
    all_row_size = new int [_MatrixSize];
  
    for( int i=0; i < _MatrixSize; i++) {

      if ( i < A.getNbrRow() ) {
	row_[i]     = new int [ A.row_size[i] ];
	
	int *ptr = row_[i] ;
	int nb_row = 0;
	for ( int* it = A.row_[i]; (it != A.row_end[i]) && (*it <  _MatrixSize); 
	      it++) {
	  *ptr = *it; ptr++;
	  nb_row++;
	}
	row_size[i] = nb_row ;
	row_end[i]  = row_[i] + nb_row;
	
	if( ! symmetric_ ) {
	  // all_row not used for unsymmetic matrix
	  all_row_[i] = 0;
	  all_row_end[i] = 0;
	  all_row_size[i] =0 ;
	  
	} else {
	  all_row_[i]     = new int [ A.all_row_size[i] ];
	  ptr = all_row_[i] ;
	  int nb_row = 0;
	  for ( int* it = A.all_row_[i]; 
		it != (A.all_row_end[i]) && (*it <  _MatrixSize); it++) {
	    *ptr = *it; ptr++;
	    nb_row++;
	  }
	  all_row_size[i] = nb_row ;
	  all_row_end[i]  = all_row_[i] + nb_row;
	}
      
	col_[i]     = new int [  A.col_size[i] ];
	int nb_col =0;
	ptr = col_[i] ;
	for (  int* it = A.col_[i]; it != A.col_end[i] && (*it < _MatrixSize); 
		 it++) {
	    *ptr = *it; ptr++;
	    nb_col++;
	  }
	  col_size[i] = nb_col;
	  col_end[i]  = col_[i] + nb_col;
      } else {
	  row_[i] = 0;
	  row_end[i] = 0;
	  row_size[i] =0 ;
	  col_[i] = 0;
	  col_end[i] = 0;
	  col_size[i] =0 ;
	  all_row_[i] = 0;
	  all_row_end[i] = 0;
	  all_row_size[i] =0 ;
      }
    }
    symmetric_            = A.symmetric_;
    updateNbrOfValues_();
      
      
    // Values not copied
    // int u = A.getNodeValuesInc();
    // const int64_t *v = A.getNodeValues();
    // setNodeValues( v, u );
    if( A.getNodeValues() != 0)

#   ifdef MSG
      cerr << "copyMatrix : Sub matrix with values : not implemented" << endl;
#   endif

    NodeValues_ =0;
  } 


  void updateNbrOfValues_() {
    int64_t n_val = 0;

    for( int i=0; i < _MatrixSize; i++) {
      n_val += row_size[i];
    }
    nbr_values_ = n_val;

    // The index values are sorted
    int n=0;
    for( int i=0; i < _MatrixSize; i++) {
      for ( int *it = row_[i] ; 
	    it != row_end[i]; it++ ) {
	if( *it == i ) n++;
      }
    }
    nbr_diagonal_values_ = n;
  } 

  // initiallize (allocation and copy)
  void initRow( int i, int *it_start, int *it_end ) {
    row_[i]     = new int [ it_end  - it_start ];
    row_size[i] = it_end  - it_start;
    row_end[i]  = row_[i] + row_size[i] ;

    int *it_dest = row_[i];
    for (  int *it = it_start; 
	   it != it_end; it++, it_dest++) {
      *it_dest = *it;
    }
  }

  // initiallize (allocation and copy)
  void initColumn( int j, int *it_start, int *it_end ) {
    col_[j] = new int [ it_end  - it_start ];
    col_size[j] = it_end  - it_start;
    col_end[j]  = col_[j] + col_size[j] ;

    int *it_dest = col_[j];
    for (  int *it = it_start; 
	   it != it_end; it++, it_dest++) {
      *it_dest = *it;
    }
  }

  void inline nextNeighbours( int i , bool *connected, 
		       bool *ToDo )  { 

      ToDo[ i ] = false;
      for( int *it = all_row_[i]; it != all_row_end[i]; it++ ) {
	connected[*it] = true;
      }

      bool all = true;
      for( int k=0; (all) && (k < _MatrixSize); k++) {
	all &= connected[i];
      }
      if ( ! all ) {
	for( int *it=all_row_[i]; it != all_row_end[i]; it++ ) {
	  if( ToDo [*it] )
	    nextNeighbours( *it, connected, ToDo );
	} 
      }
  }



 public :

  //! \brief Constructor - Empty matrix.
  SparseAdjacency() { 
    // LOGInFct( std::cout );
    _MatrixSize = 0; 
    col_= 0;
    row_=0;
    col_size = 0;
    row_size = 0;
    col_end = 0;
    row_end = 0;
    symmetric_ = false;
    // LOGOutFct( std::cout );
  }  

  SparseAdjacency(const int matrix_size, bool symmetric=true ) { 
    // LOGInFct( std::cout );

    _allocateRowAndColumnDescriptors ( matrix_size );

    symmetric_ = symmetric;

    updateNbrOfValues_();
    NodeValues_ = 0;
    // LOGOutFct( std::cout );
  }
  //! \brief Constructor -  The adjacency matrix is read from a file
  //! \param filename name of the file containiNg the edge list
  //! \param symmetric the matrix is symmetric (if true unsymmetric otherwise)
  //! \param diag take into account diagonal terms.

  // ??? envoyer un flux
  SparseAdjacency( const char *filename, 
		   const bool symmetric = true, 
		   const bool diag = false, const char *format=0 ) { 
 
    ifstream infile;
    string str;
    size_t n = 0;
    string Format;
    bool error = false;

    if (format) 
      Format = format;

    //
    // Read Color file
    //
    int NbrColor = 0;
    int64_t *color = 0;

    if ( Format == "col" ) {
      string FileColorName( filename );
      size_t pos = FileColorName.find_last_of( "." );
      FileColorName.erase( FileColorName.begin()+pos, FileColorName.end()) ;
      FileColorName.append(".");
      FileColorName.append(Format);

      // Number of lines
      infile.open (FileColorName.c_str(), ifstream::in);
      while ( infile.good() ) {
	infile >> str ;
	n++;
      }

      NbrColor = n;
      color = new int64_t[NbrColor];

      // Go to the begining of the file
      infile.close();
      infile.open (FileColorName.c_str(), ifstream::in);

      string buf[2];
      int index = 0;

      while ( infile.good() ) {
	infile >> ws >> buf[0] >> ws ;

	color[ index++ ] = atoi( buf[0].c_str() );
      }
      // Ajust Number of colors
      NbrColor = index;

      // Simple check
      for ( index=0; index < NbrColor; index++) {
	if ( (color[index] < 0) &&  (color[index] >= NbrColor) ) {
	  error = true;
	}
      }
      if (error) {
#   ifdef MSG
	cerr << "Bad color file: " <<  FileColorName << endl;
#   endif
      }
      infile.close();
    }

    //
    // Read edge list file
    //
    n = 0;

    infile.open (filename, ifstream::in);
    while ( infile.good() ) {
      infile >> str ;
      n++;
    }

    size_t  NbrEdges = (n / 2);
    int *edge = new int[2*NbrEdges+2];

    // Go to the begining of the file
    infile.close();
    infile.open (filename, ifstream::in);

    int NbrLoops = 0;
    string buf[2];
    int NodeID[2];
    int index = 0;
    int max_ind1 = 0, max_ind2 = 0;

    while ( infile.good() ) {
      infile >> ws >> buf[0] >> ws >> buf[1] >> ws;

      // TODO Test lignes vides
      for (int i=0; i < 2; i++) {
 
	bool notFound=true;

	int k = 0;
	int pos = NOT_INDEX;
	while( (((size_t) k) < NodeNames.size()) && notFound ) {

	  // for( int j = 0; j < NodeNames.size(); j++) {
	  //   cout << NodeNames[j] << "_" ;
	  // }
	  // cout << endl;

	  if ( buf[i].compare( NodeNames[k]) == 0) {
	    notFound = false;
	    pos = k;
	    // cout << "Found node : " <<  buf[i].c_str() << "==" 
	    // <<  NodeNames[k] << " " <<  pos << " " << NodeNames.size() 
	    // << endl;
	  }
	  k++;
	}
	if( notFound ) { 
	  char *ptr = new char[ buf[i].length() +1 ];
	  memcpy ( ptr,  buf[i].c_str(),  buf[i].length() + 1  );
	  NodeNames.push_back( ptr );
	  pos = NodeNames.size() -1;
	  /// cout << "New node : " <<  buf[i].c_str() 
	  //     << " " <<  pos << " " << NodeNames.size() 
	  //     << " > "  << NodeNames[pos] << " " 
	  //     << buf[i].length() << endl;
	}
	NodeID[i] = pos;
      }

      // cout << NodeID[0] << " " << NodeID[1] << " " 
      //      << NodeNames.size() << endl;

      if ( ( diag ) || (NodeID[0] != NodeID[1]) ) {
	edge[index*2]   = NodeID[0];
	edge[index*2+1] = NodeID[1];
	max_ind1 = MAX( max_ind1, NodeID[0] );
	max_ind2 = MAX( max_ind2, NodeID[1]);

	// cout << "test  " << index << " "
	//      <<  NodeID[0]<< "-" <<  NodeID[1] << endl;

	index++;
      } else {
	NbrLoops++;
	// cout << "loop removed " << NodeID[0] << "==" << NodeID[1] << endl;
      }

    }

    // cout << "test " << index << endl;
    edge[index*2]   = NOT_INDEX;
    edge[index*2+1] = NOT_INDEX;
    infile.close();

    // square_matrix
    max_ind1 = MAX( max_ind1, max_ind2);
    max_ind2 = max_ind1;

#   ifdef MSG
    cout << "File " << filename << " read :" << endl;  
    cout << "  Number of edges : " << NbrEdges << endl; 
    cout << "  Number of nodes : " << NodeNames.size() << endl;
    if (symmetric) {
      cout << "  Undirected matrix" <<  
	"(" <<  max_ind1+1 << ","  << max_ind2+1 << ")" << endl;
      if( diag ) 
	cout << "  Loops taking into account" << endl;
      else {
	cout << "  Loops removed (" << NbrLoops << " loops)" << endl;
	cout << "  Remain edges : " << NbrEdges - NbrLoops  << endl;
      }
    }
    else {
      cout << "  Directed matrix" <<  
	"(" <<  max_ind1+1 << ","  << max_ind2+1 << ")" << endl;
      if( diag ) 
	cout << "  Loops taking into account" << endl;
      else {
	cout << "  Loops removed (" << NbrLoops << " loops)" << endl;
	cout << "  Remain edges : " << NbrEdges - NbrLoops  << endl;
      }
    }
#   endif

    // Print node names
    // for (int it = 0 ; it < NodeNames.size(); it++) {
    //  cout << it << "->" <<   NodeNames[it] << endl;
    // } 

    if (!error)
      initMatrix( edge, max_ind1+1 , symmetric, color);

    
    delete [] color;
    delete [] edge;

  }

  //
  SparseAdjacency(const int matrix_size, const char *name, bool symmetric=true ) {
 
    // LOGInFct( std::cout );

    _allocateRowAndColumnDescriptors ( matrix_size );

    symmetric_ = symmetric;

    string motifName( name );
    // Clique without loop
    if( motifName.compare("clique") == 0 ) {
      if( symmetric_ ) {
	for ( int i = 0; i < (_MatrixSize - 1); i++) {
	  row_[i]     = new int [ _MatrixSize - i -1 ];
	  row_size[i] = _MatrixSize - i - 1 ;
	  row_end[i]  = row_[i] + (_MatrixSize - i - 1);
	  int *ptr = row_[i];
	  for ( int j = i+1 ; j < _MatrixSize ; j++, ptr++) {
	    *ptr = j;
	  }
	}
	row_    [ _MatrixSize - 1 ] = 0;
	row_end [ _MatrixSize - 1 ] = 0;
	row_size[ _MatrixSize - 1 ] = 0;
	
	for ( int i = 0; i < _MatrixSize ; i++) {
	  col_[i]     = new int [ i ];
	  col_size[i] = i;
	  col_end[i]  = col_[i] + (i);
	  int *ptr = col_[i];
	  for ( int j = 0 ; j < i; j++, ptr++) {
	    *ptr = j;
	  }
	}
	for ( int i = 0; i < _MatrixSize; i++) {
	  all_row_[i]     = new int [ _MatrixSize -1 ];
	  all_row_size[i] = _MatrixSize - 1 ;
	  all_row_end[i]  = all_row_[i] + (_MatrixSize - 1);
	  int *ptr = all_row_[i];
	  for ( int j = 0 ; j < _MatrixSize ; j++) {
	    if( i != j ) {
	      *ptr++ = j;
	    }
	  }
	}
      } else {
	for ( int i = 0; i < _MatrixSize; i++) {
	  row_[i]     = new int [ _MatrixSize -1 ];
	  row_size[i] = _MatrixSize -1 ;
	  row_end[i]  = row_[i] + (_MatrixSize - 1);
	  int *ptr = row_[i];
	  for ( int j = 0 ; j < _MatrixSize ; j++) {
	    if( i != j ) {
	      *ptr++ = j;
	    }
	  }
	}
	// Indexes must be exchanges
	for ( int i = 0; i < _MatrixSize; i++) {
	  col_[i]     = new int [ _MatrixSize -1 ];
	  col_size[i] = _MatrixSize -1 ;
	  col_end[i]  = col_[i] + (_MatrixSize - 1);
	  int *ptr = col_[i];
	  for ( int j = 0 ; j < _MatrixSize ; j++) {
	    if( i != j ) {
	      *ptr++ = j;
	    }
	  }
	}

      } 
    }
    
    updateNbrOfValues_();
    NodeValues_ = 0;
    // LOGOutFct( std::cout );
  }

  // \brief Copy constructor
  SparseAdjacency( const SparseAdjacency &A) {
    copyMatrix_( A );
  } 

  //! \brief Copy sub_matrix or larger matrix
  //! Values are not copied
  SparseAdjacency( const SparseAdjacency &A, int matrix_size ) {
    if ( (matrix_size <= 0) )
      copyMatrix_( A );
    else
      copyMatrix_( A, matrix_size );
  } 

  //! \brief Constructor -  The adjacency matrix is filled with
  //!        the edje list (i_index, j_index).
  //!  WARNING : matrix_size must be greater than max(index) + 1
  //! \param ij_index  Array of index pairs (i,j) (array[2, ...], i and j >= 0).
  //!                  This array must be ended with the values (NOT_INDEX, NOT_INDEX)
  //! \param symmetric Reduce the matrix to a upper-triangular matrix (i <= j).
  //!                 If necessary a matrix element (if i>j) could be transposed
  //!                 to obtain an upper-triangular form.
  //! 
  // TODO 
  //             else "ptr" is directly assign to the internal matrix data.
  //             In "false" case do not use the pointer for other 
  //              computations 
  SparseAdjacency(   const int ij_index[], 
		     const int matrix_size, 
		     const bool symmetric = false,
		     const int64_t *NodeValues = 0,
		     size_t NodeValuesInc = 1) { 


    // LOGInFct( std::cout );
    
    initMatrix( ij_index, matrix_size, symmetric, NodeValues, NodeValuesInc );

    // LOGInFct( std::cout );
  }

  //! \brief Constructor -  The adjacency matrix is filled with
  //!        a dense matrix.
  //! \param matrix  boolean adjacency matrix 
  //! \param symmetric Reduce the matrix to a upper-triangular matrix (i <= j).
  //!                 If necessary a matrix element (if i>j) could be transposed
  //!                 to obtain an upper-triangular form.
  //! 
  SparseAdjacency(   const bool mat[],
		     const int  matrix_size, 
		     const bool symmetric = false,
		     const int64_t *NodeValues = 0,
		     size_t NodeValuesInc = 1) { 

    // LOGInFct( std::cout );

    symmetric_ = symmetric;
    _MatrixSize = matrix_size ;

    _allocateRowAndColumnDescriptors( _MatrixSize );

    _fillDenseMatrix( mat, _MatrixSize, symmetric);

    updateNbrOfValues_();  

    setNodeValues( NodeValues, NodeValuesInc );

    // LOGInFct( std::cout );

  }

  //! \brief Constructor -  The adjacency matrix is filled with
  //!        a dense matrix.
  //! \param matrix  integer adjacency matrix 
  //! \param symmetric Reduce the matrix to a upper-triangular matrix (i <= j).
  //!                 If necessary a matrix element (if i>j) could be transposed
  //!                 to obtain an upper-triangular form.
  //! 
  SparseAdjacency(   const int mat[],
		     const int  matrix_size, 
		     const char *type, 
		     const bool symmetric = false,
		     const int64_t *NodeValues = 0,
		     size_t NodeValuesInc = 1) { 

    string str( type );
    // LOGInFct( std::cout );

    if( str.compare("dense") == 0 ) {
      symmetric_ = symmetric;
      _MatrixSize = matrix_size ;
      
      _allocateRowAndColumnDescriptors( _MatrixSize );
      
      bool *bool_mat = new bool[_MatrixSize * _MatrixSize];
      bool *it_bool=bool_mat;
      for( const int *it = mat, *it_end = mat +  _MatrixSize * _MatrixSize; 
	   it != it_end; it++, it_bool++) {
	*it_bool = (*it != 0);
      }
      
      _fillDenseMatrix( bool_mat, _MatrixSize, symmetric);
      
      delete [] bool_mat;

      updateNbrOfValues_();  

      setNodeValues( NodeValues, NodeValuesInc );
    }

    // LOGOutFct( std::cout );

  }

  //! \brief read the matrix (dense format) in a stream
  //!  the file is filed with "0" and "1" 
  //! \param infile ifstream
  //! \param sep sperator betwen to values
  //! \param buffer_size maximum size of a file line.
  SparseAdjacency( ifstream &infile, const char *sep, 
		   const bool symmetric=true, 
		   int buffer_size=256 ) {

    char *buffer = new char[ buffer_size ];
    int inc = strlen( sep ) + 1;
    int mat_size = 0;
    int i = 0;
    bool GoOn = true;
    bool *mat = 0;

    while( GoOn ) { 
      infile.getline( buffer, buffer_size ); 
      if (mat == 0) {
	mat_size = strlen( buffer ) / inc;
	mat      = new bool[mat_size * mat_size];
      }
      if (mat != 0) {
	char *ptr_r =  buffer;
	for( bool *ptr = &mat[i], *ptr_end = ptr + mat_size*mat_size; 
	     ptr != ptr_end; ptr += mat_size, ptr_r += inc) {
	  *ptr = (*ptr_r == '1');
	}
	i++;
	GoOn = (i < mat_size ) && infile.good() ;
      }
    }

    if( i == mat_size ) {
      symmetric_  = symmetric;
      _MatrixSize = mat_size ;

      _allocateRowAndColumnDescriptors( _MatrixSize );
      
      _fillDenseMatrix( mat, _MatrixSize, symmetric);

      updateNbrOfValues_();

      // No values on nodes
      setNodeValues();

      delete [] mat;
      delete [] buffer;
    } else {
      symmetric_  = symmetric;
      _MatrixSize = 0 ;

      delete [] mat;
      delete [] buffer;
    }
  }


  // \brief Destructor
  ~SparseAdjacency() { 
    _deleteRowAndColumnContents( ); 
    _deleteRowAndColumnDescriptors();
    _deleteNodeValues() ; 
  };



  //! \brief Get the number of rows in Matrix.
  //! \return Return the row number.
  int inline getNbrRow( ) const { return ( _MatrixSize ); } 

  //! \brief Get the number of rows in Matrix.
  //! \return Return the row number.  
  int inline getNbrCol( ) const { return ( _MatrixSize ); }

  //! \brief Get the Matrix symmetry structure.
  //! \return Return true if symmetric data structure.  
  bool inline getSymmetry( ) const { return ( symmetric_ ); }

  //! \brief Get the number of non-zero values.
  //! If the matrix is symmetric, return the 
  //! number of non-zero values of the upper triangular 
  //! part (see getNbrOfMatrixValues to get the total number).
  //! \return Return the total number of pairs (i,j) .  
  // ??? integer  
  int64_t inline getNbrOfValues( ) const {
    return( nbr_values_ ); 
  }

  //! \brief Get the number of non-zero values.
  //! If the matrix is symmetric, return the 
  //! number of non-zero values of the whole matrix 
  //! (see getNbrOfValues to get the non-zero values 
  //! of the upper triangular part).
  //! \return Return the total number of pairs (i,j) .  
  int inline getNbrOfMatrixValues( ) const {
    return( (static_cast<int>(symmetric_) ) * (nbr_values_ - nbr_diagonal_values_) 
	    + nbr_values_ ); 
  }

  //! \brief Get a list of non-zero element in Matrix (edge list).
  //! \return Return a list array of row indexes (list<int>[nbr_of_columns]) .  
  //! To Test ??? for unsymmetric with diag = true
  int inline *getIndexes( bool cmpl = false, bool diag = true ) const { 

    int *ret_array;
    int k = 0;

    if (  !cmpl ) {
      if( ! diag ) {
#   ifdef MSG
	cerr << "Not implemented" << endl;
#   endif
      }
      ret_array = new int [ 2 *(nbr_values_ + 1)];
      
      for( int i=0; i < _MatrixSize; i++) {
	for ( int* it = row_[i] ; it != row_end[i]; it++ ) {
	  ret_array[k++] = *it;
	  ret_array[k++] = i;
	}
      }
    } else {

      int nbr_edges = _MatrixSize * (_MatrixSize - 1) / 2 ;
      if( ! symmetric_ ) {
	nbr_edges = nbr_edges * 2;
      }
      if( diag ) {
	nbr_edges += _MatrixSize;
      }

      nbr_edges = nbr_edges - nbr_values_;

      ret_array = new int [ 2 *( nbr_edges + 1)];
      
      int NoIndex = NOT_INDEX;
      for( int i=0; i < _MatrixSize; i++) {
	int* it = row_[i];
	if(  it == row_end[i] ) 
	  it = &NoIndex; 

	int j=0;
	if (symmetric_) {
	  j = i;
	}
	for( ; j < _MatrixSize; j++) {
	  if( *it != j  ) {
	    if( (diag) || (i != j) ) {
	      // Must be switched (arc: j -> i)
	      ret_array[k++] = j;
	      ret_array[k++] = i;
	    }
	  } else {
	    // Next element
	    if( it != row_end[i] ) 
	      it++; 
	  }
	}
      }
    }
    ret_array[k] = NOT_INDEX;
    ret_array[k] = NOT_INDEX;

    return( ret_array );
  }

  //! \brief Get the adjacency Matrix.
  //! \return a bool array[nbr_element * nbr_element] .  
  //! To Test ??? for unsymmetric with diag = true
  bool inline *getDenseMatrix( ) const { 
    
    bool *mat = new bool[_MatrixSize * _MatrixSize];
    
    for( bool *ptr=mat, *ptr_end=&mat[_MatrixSize * _MatrixSize]; 
	 ptr != ptr_end; ptr++) {
      *ptr = false;
    }

    for( int i=0; i < _MatrixSize; i++) {
      for ( int* it = row_[i]; it != row_end[i]; it++) {
	mat[i + (*it)*_MatrixSize ] = true;
      }
    }
    // Lower part
    if (symmetric_){
      for( int j=0; j < _MatrixSize; j++) {
	for ( int* it = row_[j]; it != row_end[j]; it++) {
	  mat[(*it) + j*_MatrixSize ] = true;
	}
      }
    }
    return mat;
  }

  //! \brief Get a matrix element value.
  //! \param i row indice.
  //! \param j column indice.
  //! \return Return A(i,j) value.
  bool inline getMatrix(int i, int j) const { 
    if ( symmetric_ && (i > j) ) { int k = j; j = i; i = k; }
    bool found = false;
    for ( int* it = row_[i]; 
	  it != row_end[i] && !(found); it++ ) {
      found = ( *it == j ); 
    }
    return found;
  };

  //! \brief Get a dense (sub)adjacency matrix.
  //! \param k row, colomn indices.
  //! \return Return A(i,j) value i,j in set {index}.
  int *getMatrix(int *index, int n_index, bool transpose) const { 

    int *mat = new int[n_index*n_index];
    for (int *it=mat, *it_end=mat + n_index*n_index; it != it_end; ) {
      *it++=0;
    }

    int i, j , i_tmp;
    int i_index = 0;
    for ( int* it_i = index; i_index < n_index ; it_i++, i_index++ ) {
      int j_index = 0;
      for ( int* it_j = index; j_index < n_index; it_j++, j_index++ ) {
	i=*it_i; j=*it_j;
	if ( (symmetric_ && (i < j)) || (! symmetric_) ) {
	  bool found = false;
	  for ( int* it = row_[i]; 
		it != row_end[i] && !(found); it++ ) {
	    found = ( *it == j ); 
	  }
	  if( transpose) {
	    i_tmp = i;
	    i = j;
	    j = i_tmp;
	  }
	  if (found) {
	    mat[i_index + j_index*n_index] = 1; 
	    if( symmetric_ ) 
	      mat[j_index + i_index*n_index] = 1; 
	  } 
	}
      }
    }

    return mat;
  };


  //! \brief Get the sums of non zero elements for each row in the 
  //! Matrix (degress).
  //! Return the result of the upper-triangle values for a symmetric
  //! matrix.
  //! \return Return an NRow-size array. 
  const int *getRowSums( ) const { 
    return( row_size  );
  }

  //! \brief Get the sums of non zero elements for each row in the 
  //! Matrix (degress).
  //! Return the result of the whole matrix (upper and lower triangular) 
  //! values for a symmetric matrix.
  //! \return Return an NRow-size array. 
  const int *getMatrixRowSums( ) const {
    if (symmetric_) 
      return( all_row_size  );
    else 
      return( row_size  );
  }

  //! \brief Same as getRowSums.
  //! \return Return an NColumn-size array. 
  const int inline *getColSums( ) const { 
    return( col_size );
  }

  //! \brief Get the sum of non zero elements for the row i  
  //! in the Matrix (degree).
  //! Return the result of the upper-triangle values for a symmetric
  //! matrix.
  //! \return Return a single value. 
  int inline getRowSize( int i ) const { 

    return( row_size[i] );
  }

  //! \brief Get the sum of non zero elements for the whole row i  
  //! in the Matrix (degree).
  //! Return the result on the whole matrix (even for a symmetric one)
  //! 
  //! \return Return a single value. 
  int inline getAllRowSize( int i ) const { 

    return( all_row_size[i] );
  }

  //! \brief Get the sum of non zero elements for the column j  
  //! in the Matrix (degree).
  //! Return the result of the upper-triangle values for a symmetric
  //! matrix.
  //! \return Return a single value. 
  int inline getColSize( int j ) const { 

    return( col_size[j] );
  }


  //! \brief Get the i-row  in the Matrix.
  //! \return Return the list row=i of column indexes (const list<int> *)  .
  //! ??? TO TEST  
  const int inline *getRow( int i ) const { 
    return( row_[i] );
  }

  //! \brief Get the whole i-row  in the Matrix.
  //! \return Return the list row=i of column indexes + row indexes 
  //! (const list<int> *)  .
  //! ??? TO TEST  
  const int inline *getAllRow( int i ) const { 
    return( all_row_[i] );
  }

  //! \brief Get the j-col  in the Matrix.
  //! \return Return the list row=i of column indexes (const list<int> *)  .
  //! ??? TO TEST  
  const int inline *getCol( int j ) const { 
    return( col_[j] );
  }

  //! \brief Get the end pointer of  the i-row array in the Matrix.
  //! \return Return a const int*  .
  const int inline *getRowEnd( int i ) const { 
    return( row_end[i] );
  }

  //! \brief Get the end pointer of  the whole i-row array in the Matrix.
  //! \return Return a const int*  .
  const int inline *getAllRowEnd( int i ) const { 

    return( all_row_end[i] );
  }

  //! \brief Get the end pointer of  the j-col array in the Matrix.
  //! \return Return a const int*  .
  //! ??? TO TEST  
  const int inline *getColEnd( int j ) const { 

    return( col_end[j]  );
  }


  //! \brief Get node values.
  //! \return Return a reference on the first vector describing the
  //!   first node.  
  const int64_t *getNodeValues( ) const { return NodeValues_;};

  //! \brief Get the vector size describing one node .
  //! \return Return the vector size.  
  int getNodeValuesInc( )  const { return NodeValuesInc_;};

  bool inline isInAllRow( int i, int value ) {

    int *p_beg = all_row_[i] ;
    int *p_end = all_row_end[i] ;
    int *p_mid = p_beg + ( p_end - p_beg ) / 2 ;
    bool found = false;

    if (  p_beg != p_end ) {
      for ( ; p_beg != p_mid; ) {
	
	if( *p_mid < value ) {
	  p_beg = p_mid;
	  p_mid = p_beg + ( p_end - p_beg ) / 2 ;

	} else if ( *p_mid > value ) {
	  p_end = p_mid;
	  p_mid = p_beg + ( p_end - p_beg ) / 2 ;
	} else {
	  p_beg = p_mid;
	  p_end = p_mid;
	  found = true;
	}
      }

      if( !found ) {
	if (p_mid == all_row_[i])  {
	  found = (*p_mid == value);
	} else if ( p_mid == all_row_end[i])  {
	  found = (*p_mid == value);
	}
      } 
    }
    return found;
  } 

  bool inline isInRow( int i, int value ) {

    int *p_beg = row_[i] ;
    int *p_end = row_end[i] ;
    int *p_mid = p_beg + ( p_end - p_beg ) / 2 ;
    bool found = false;

    if (  p_beg != p_end ) {
      for ( ; p_beg != p_mid; ) {
	
	if( *p_mid < value ) {
	  p_beg = p_mid;
	  p_mid = p_beg + ( p_end - p_beg ) / 2 ;
	  
	} else if ( *p_mid > value ) {
	  p_end = p_mid;
	  p_mid = p_beg + ( p_end - p_beg ) / 2 ;
	} else {
	  p_beg = p_mid;
	  p_end = p_mid;
	  found = true;
	}
      }

      if( !found ) {
	if (p_mid == row_[i])  {
	  found = (*p_mid == value);
	} else if ( p_mid == row_end[i])  {
	  found = (*p_mid == value);
	}
      } 
    }
    return found;
  } 

  bool inline isInCol( int i, int value ) {

    int *p_beg = col_[i] ;
    int *p_end = col_end[i] ;
    int *p_mid = p_beg + ( p_end - p_beg ) / 2 ;
    bool found = false;

    if (  p_beg != p_end ) {
      for ( ; p_beg != p_mid; ) {

	if( *p_mid < value ) {
	  p_beg = p_mid;
	  p_mid = p_beg + ( p_end - p_beg ) / 2 ;
	  
	} else if ( *p_mid > value ) {
	  p_end = p_mid;
	p_mid = p_beg + ( p_end - p_beg ) / 2 ;
	} else {
	  p_beg = p_mid;
	  p_end = p_mid;
	  found = true;
	}
      }
      if( !found ) {
	if (p_mid == col_[i])  {
	  found = (*p_mid == value);
	} else if ( p_mid == col_end[i])  {
	  found = (*p_mid == value);
	}
      } 
    }
    return found;
  } 

  bool isConnex() {
    bool *ToDo      = new bool[_MatrixSize];
    bool *connected = new bool[_MatrixSize];

    for( int i=0; i < _MatrixSize; i++) {
      ToDo[i] = true;
      connected[i] = false;
    }
    nextNeighbours( 0, connected, ToDo );
    bool all = true;
    for( int i=0; (all) && (i < _MatrixSize); i++) {
      all &= connected[i];
    }

    delete [] ToDo;
    delete [] connected;

    return all;
  }

  // Undirected matrix are store in triangular matrix
  // format. Thus, for some cases, a node may not be 
  // connected with previous ones 
  bool optimizeConnexity( int *perm, int *array = 0) {

    // Degree array :
    //  >=0 : degree
    //  -1 : degree already ordered
    bool connex = true;

    // Degre of the node or -1 if already processed
    int *deg_ = new int[_MatrixSize]; 
    // Node order
    int *order = new int[_MatrixSize];

    for( int i=0; i < _MatrixSize; i++) {
      deg_[i] = 0;
    }

    // Number of nodes ordered
    int k=0;
    // Get the degre max
    int max_deg;
    if( symmetric_ ) {
      max_deg =  all_row_size[0];
    } else {
      max_deg =  row_size[0] + col_size[0];
    }
    int ind_max = 0;
    if( symmetric_ ) {
      for( int i=1; i < _MatrixSize; i++) {
	if ( max_deg < all_row_size[i] ){
	  max_deg = all_row_size[i];
	  ind_max = i;
	}
      }
    } else {
      // Outcoming + Incoming arcs
      for( int i=1; i < _MatrixSize; i++) {
	if ( max_deg < (row_size[i] + col_size[i])){
	  max_deg =  row_size[i] + col_size[i];
	  ind_max = i;
	}
      }
    }
    
    // Store it the rearragement array the 
    // max degre node
    order[k] = ind_max;
    perm[ind_max] = k;
    k++;
    deg_[ind_max] = -1;
    
    // Select the nodes which is the most connected
    // with the already selected nodes
    for( int m = 1; m < _MatrixSize; m++) {

      // To find the maximun degree
      max_deg = 0;
      ind_max = -1;
      for( int l = 0; l < _MatrixSize; l++ ) {
	// Not already ordered
	if( deg_[l] >= 0 ) {
	  deg_[l] = 0;
	  if ( symmetric_ ) {
	    for( int i=0; i < all_row_size[l]; i++) {
	      // Connected with selected nodes in "order"
	      if( deg_[ all_row_[l][i] ] == -1 ) {
		deg_[ l ] = deg_[ l ] + 1;
	      }
	    }
	  } else {
	    // Outcoming
	    for( int i=0; i < row_size[l]; i++) {
	      // Connected with selected nodes in "order"
	      if( deg_[ row_[l][i] ] == -1 ) {
		deg_[ l ] = deg_[ l ] + 1;
	      }
	    }	      
	    // Incomming
	    for( int i=0; i < col_size[l]; i++) {
	      // Connected with selected nodes in "order"
	      if( deg_[ col_[l][i] ] == -1 ) {
		deg_[ l ] = deg_[ l ] + 1;
	      }
	    }	      
	  }
	}
	if ( max_deg < deg_ [ l ] ){
	  max_deg = deg_[ l ];
	  ind_max =  l ;
	}
      }
      if ( max_deg > 0 ) {
	order[k] = ind_max;
	perm[ ind_max ] = k;
	k++;
	deg_[ ind_max ] = -1;
      } else {
	// Not a connex graph
	order[0] = -1;
	connex = false;
      }
    }
 
    if ( connex ) {
      permute( perm );
      
      if( array )
	_permute( perm, array, _MatrixSize );
    }

    delete [] deg_; 
    delete [] order; 

    return connex;

  }

  //! NodeNames is not token into account
  void permute( const int *sigma ) {
    
    int *ij_index = getIndexes( );

    for( int *it = ij_index; *it != NOT_INDEX;  it++) {
      *it = sigma[ *it ];
    }

    _deleteRowAndColumnContents();
    _copyRowAndColumnContents( ij_index );

    if ( NodeValues_ ) {
      _permute<int64_t>( sigma, NodeValues_, _MatrixSize, NodeValuesInc_);
    }

    delete [] ij_index;
  }

  // The permutation size is smaller then the matrix size 
  void permute( const int *sigma, int sigma_size ) {
    
    int *ij_index = getIndexes( );

    for( int *it = ij_index; *it != NOT_INDEX;  it++) {
      if(  *it <  sigma_size ) {
	*it = sigma[ *it ];
      }
    }

    _deleteRowAndColumnContents();
    _copyRowAndColumnContents( ij_index );

    if ( NodeValues_ ) {
      _permute<int64_t>( sigma, NodeValues_, _MatrixSize, NodeValuesInc_);
    }

    // ??? NodeNames
 
    delete [] ij_index;
  }

  //! \brief Shift all the matrix indexes.
  //!        Do a permutation.
  //! \param start index to start the shift.
  //! \param len   number of indexes to shift.  
  //! \param shift shifting value.
  //! NodeValues is not token into account
  void shiftNumbering( const int start, const int len, 
		       const int shift) {

    int max_val = len + start;
    int  *sigma = new int[ max_val];
    for( int i = 0; i < start; i++) {
      sigma[i] = i;
    } 
    for( int i = start; i < max_val; i++) {
      sigma[i] = i + shift;
    } 
    
    permute( sigma, max_val );
    
    // To Do
    delete [] sigma;
  }

  //! NodeNames and NodeValues are not token into account
  void setMatrix( const SparseAdjacency &A) {
    
    
    if( symmetric_ == A.getSymmetry() ){

      int *ij_index = getIndexes( );
      int *A_ij_index = A.getIndexes( );
      size_t nb_values = 0;
      size_t A_nb_values = 0;
      
      for( int *it = ij_index; *it != NOT_INDEX;  it++) {
	nb_values++;
      }
      for( int *it = A_ij_index; *it != NOT_INDEX;  it++) {
	A_nb_values++;
      }
      
      int *new_ij_index = new int[ nb_values + A_nb_values +2];
      
      int *ptr = new_ij_index;
      
      for( int *it = ij_index; *it != NOT_INDEX;  it++, ptr++) {
	*ptr = *it;
      }
      for( int *it = A_ij_index; *it != NOT_INDEX;  it++, ptr++) {
	*ptr = *it;
      }
      ptr[0] = NOT_INDEX;
      ptr[1] = NOT_INDEX;
      
      _deleteRowAndColumnContents();
      _copyRowAndColumnContents( new_ij_index );
      
      updateNbrOfValues_();   
      
      // NodesNames ???

 #   ifdef MSG     
      if( NodeValues_ != 0) 
	cerr << "setMatrix with values : not implemented" << endl;
#   endif

      delete [] A_ij_index;
      delete [] ij_index;
      delete [] new_ij_index;
      
    } else {
#   ifdef MSG
      LOGMSG( 1, std::cout, "Not the same matrix mapping (symetry)", "" );
#   endif
    }
  }

  //! \brief Remove the arc i -> j (i.e. the A[ dest, src ] 
  //!        element).  
  //! \param src source node
  //! \param dest destination node
  //! \return Return true if the arc is found
  //! ??? TO TEST  
  bool removeArc( int src, int dest ) {

    if( symmetric_ ) {

      // For row dest <= src
      if( src < dest ) {
	int tmp = src;
	src = dest;
	dest = tmp;
      }
    }

    bool found = _removeInRow( dest, src );

    if( found ) {
      _removeInCol( src, dest);
      nbr_values_--;
      if( src == dest )
	nbr_diagonal_values_--;
    }
    return found;
  } 

  //! \brief Get the triangular upper or lower part of the matrix.
  //! \param sup specifies the upper part otherwise the lower 
  //!   part is returned.
  //! \return Return a SparseAjacency
  //! ??? TO TEST  
  SparseAdjacency getTriangular( bool sup = true ) { 

    SparseAdjacency Tr( _MatrixSize, false );
    
    if( symmetric_ ) {   
      if ( sup ) {
	//
	// Triangular Sup
	//
	for( int i=0; i < _MatrixSize; i++) {
	  Tr.initRow( i, row_[i], row_end[i]);
	}
	for( int j=0; j < _MatrixSize; j++) {
	  Tr.initRow( j, col_[j], col_end[j]);
	}
	
      } else { 
	//
	// Triangular Inf
	//
	for( int i=0; i < _MatrixSize; i++) {
	  Tr.initRow( i, col_[i], col_end[i]);
	}
	for( int j=0; j < _MatrixSize; j++) {
	  Tr.initColumn( j, row_[j], row_end[j]);
	}
      }
    } else {

      int *start_it;
      int *end_it;

      if ( sup ) {
	//
	// Triangular Sup
	//

	for( int i=0; i < _MatrixSize; i++) {
	  if ( row_size[i] !=0 ) {
	    start_it = row_[i] ;
	    for (  int *it = row_[i]; 
		  it != row_end[i] && (*it < i); ) {
	      it++;
	      start_it = it;
	    }
	    Tr.initRow(i, start_it, row_end[i] );
	  }
	}
	for( int j=0; j < _MatrixSize; j++) {
	  if ( col_size[j] !=0 ) {
	    end_it =  col_[j];
	    for ( int *it = col_[j]; 
		  it != col_end[j] && (*it <= j); ) {
	      it++;
	      end_it = it;
	    }
	    Tr.initColumn( j, col_[j], end_it );
	  }
	}
      } else  {
	//
	// Triangular Inf
	//

	for( int i=0; i < _MatrixSize; i++) {
	  if ( row_size[i] !=0 ) {
	    end_it = row_[i] ;
	    for ( int *it = row_[i]; 
		  it != row_end[i] && (*it <= i); ) {
	      it++;
	      end_it = it;
	    }
	    Tr.initRow(i, row_[i], end_it );
	  }
	}
	for( int j=0; j < _MatrixSize; j++) {
	  if ( col_size[j] !=0 ) {
	    start_it =  col_[j];
	    for ( int *it = col_[j]; 
		  it != col_end[j] && (*it < j); ) {
	      it++;
	      start_it = it;
	    }
	    Tr.initColumn( j, start_it, col_end[j] );
	  }
	}
      }
    }
    updateNbrOfValues_();  

    Tr.setNodeValues( getNodeValues(), getNodeValuesInc() );

    return( Tr );
  }
  // -----------
  //  Operators
  // -----------
  //! \brief Operator which tests if 2 two matrix are equal
  //! WARNING : If two matrix have two different kind of storage (symmetric 
  //!           and not symmetric) then the operator return false.
  //! \return  a boolean is returned.
  bool inline operator==(const SparseAdjacency & A);

  //! \brief Operator which tests if the matrix elements are present in a
  //!        other one. Used if one motif is included in other one. 
  //! WARNING : If two matrix have two different kind of storage (symmetric 
  //!           and not symmetric) then the operator return false.
  //! \return  a boolean is returned.
  bool inline operator<=(const SparseAdjacency & A);

  // ???
  // Display all matrix
  // void print( const char *name, const char *mode="matrix");
  void print(  const char *name, const char *mode ="matrix") {
  
#   ifdef MSG
    std::cout << getContext() 
	      << std::endl;
    std::cout << getContext() 
	      <<  "Matrix " << name << "(" << _MatrixSize << "," << _MatrixSize << ")" 
	      << std::endl;
    std::cout << getContext() << std::endl;
    
    int *it;

    // Node Values
    if ( NodeValues_ != 0 ) { 
      std::cout << getContext() << "Node values " <<  " : " ;
      int64_t *jt_end = NodeValues_ + _MatrixSize * NodeValuesInc_;
      for (  int64_t *jt = NodeValues_ ; jt != jt_end; ) {
	cout << " (";
	for ( int i = 0; i < NodeValuesInc_; i++, jt++) {
	  cout << *jt ;
	  if( i != (NodeValuesInc_-1) ) 
	    cout << ", " ; 
	}
	  cout << ")";
      }
      std::cout << std::endl;
    }
    
    // By rows 
    for( int i=0; i < _MatrixSize; i++) {
      std::cout << getContext() << "Row " << i << " : " ;
      for ( it = row_[i] ; it != row_end[i]; it++ )
	cout << *it << ", " ;
      std::cout << std::endl;
    }
    
    // By columns 
    for( int i=0; i < _MatrixSize; i++) {
      std::cout << getContext() << "Columns " << i << " : " ;
      for ( it = col_[i] ; it != col_end[i]; it++ )
	cout << *it << ", " ;
      std::cout << std::endl;
    }
    cout << getContext() << "Symmetric       : " << symmetric_ << endl;
    
    string str( mode );
    if( str.compare("all") == 0 ) {
      if( symmetric_) {
	std::cout << getContext() << "All_Row_size " << " : " ;
	for( int i=0; i < _MatrixSize; i++) {
	  cout << all_row_size[i] << ", " ;
	}
	std::cout << std::endl;

	// By rows 
	for( int i=0; i < _MatrixSize; i++) {
	  std::cout << getContext() << "All_Row " << i << " : " ;
	  for ( it = all_row_[i] ; it != all_row_end[i]; it++ )
	    cout << *it << ", " ;
	  std::cout << std::endl;
	}
      }
      cout << getContext() << "Number of edges  : " << nbr_values_ << endl;
      cout << getContext() << "number of loops  : " << nbr_diagonal_values_ << endl;
      cout << endl;
    }
#   endif
  }

  //! \brief Write the matrix (dense format) in a stream 
  //! \param ootfile ofstream
  //! \return bool if OK

  bool write( ofstream &outfile ) {

    bool *mat = getDenseMatrix();
    int size  = getNbrRow();
    writeMatrix( outfile, "", mat, &mat[size * size], size, "", true);
    delete [] mat;
    return true;
  }

};
bool inline SparseAdjacency::operator==(const SparseAdjacency& A) {

  bool equal_ = false;
  int i;

  if( A._MatrixSize == _MatrixSize ) {
    equal_ = true;
    if ( symmetric_ == A.symmetric_ ) {
      for( i=0; i < _MatrixSize; i++) {
	if( row_size[i] == A.row_size[i] ) {
	  equal_ &= equal( row_[i], row_end[i], A.row_[i]);
	} else {
	  equal_ = false;
	}
      }
    }
    else {
      equal_ = false;
    }
  } 
  return equal_;
}

bool 
inline SparseAdjacency::operator<=(const SparseAdjacency& A) {

  bool in = true;
  int i;

  int Min_size = MIN( _MatrixSize,  A._MatrixSize );
  for ( i = Min_size; i < _MatrixSize; i++) {
    in &= ( row_size[i] == 0 );
  } 
  if ( symmetric_ == A.symmetric_ ) {
    for( i=0; i < Min_size; i++) {

      in &= includes( A.row_[i], A.row_end[i],
		      row_[i], row_end[i]
		      );
    }
  }
  else {
    in = false;
  } 
  return in;
}


template <typename T> 
void _permute( const int *sigma, T *array, int array_size, int inc) {

  T *old_values = new T[array_size*inc];

  // Copy
  T *it     = old_values;
  T *it_end = it + array_size * inc;
  T *it_src = array;
  for( ; it != it_end; it++, it_src++) {
    *it = *it_src ;
  }
  
  // Permute
  if( inc == 1) {
    for( int i = 0; i < array_size; i++) {
      array[ sigma[i] ] = old_values[ i ];
    }
  } else {
    for( int i = 0; i < array_size; i++) {
      it = &array[ sigma[i] * inc ];
      it_end = it + inc;
      it_src = &old_values[ i * inc ];
      
      for( ; it != it_end; it++, it_src++) {
	*it = *it_src ;
      }
    }
  }
  delete [] old_values;   
}


//! \brief Transforms a squared matrix into an index list of non-zeo values. 
//!        Return an array of indexes with the following structure:
//!        (i1, j1, i2, j2, ..., ik, jk, MAX_SIZE_T, MAX_SIZE_T)
//! \param matrix array of <T> values.
//! \param matrix_size number of rows of the matrix.
//! \param symetric specifies if the storage must take account 
//!        of the symetry symetric (if true) or not (if false).
//!        If symetric only the upper triangular part is read.  
//! \return array of pair (i,j) indexes. Last values contain 
//!         the pair (MAX_SIZE_T,MAX_SIZE_T). This array must
//!         freed by the consumer.            
template<typename T>
int *getEdgeList( const T *matrix, const size_t matrix_size, 
		    const bool symetric=false ) {

  // Directed take account only of the upper-triangular part
  int count = 0;
  int iend = matrix_size;

  // Count non zero values
  for(size_t j=0; j < matrix_size; j++) {

    if( symetric )
      iend = j + 1;
   
    for(int i=0; i < iend; i++) {
      if( matrix[ j*matrix_size + i ] != 0 ) {
	count++;
      }
    }
  }

  int *edge_list = new int[ 2*(count + 1) ];

  // Fill edge list
  count = 0;
  for(size_t j=0; j < matrix_size; j++) {

    if( symetric )
      iend = j + 1;
   
    for(int i=0; i < iend; i++) {
      if( matrix[ j*matrix_size + i ] != 0 ) {
	edge_list[ 2*count ]    = i ;
	edge_list[ 2*count + 1] = j ;
	count++;
      }
    }
  }

  // End the list
  edge_list[ 2*count ]    = NOT_INDEX ;
  edge_list[ 2*count + 1] = NOT_INDEX ;

  return edge_list;
}

#endif



