# include <vector>
# include "sparse-adjacency.h"

# ifndef _MOTIF_TREE_H_
# define _MOTIF_TREE_H_

using namespace std;

class MotifTree;

typedef vector<MotifTree*> siblings_t;

class MotifTree {
  
 public :
  int *_col; // used also as motif id pointer
  int *_col_end;
  int  _col_size;
  
  int *_row;
  int *_row_end;
  int  _row_size;
  vector<int*>  _constraint;
  vector<int>  _flux_out;
  vector<int>  _flux_in;
  int*  _flux_in_t;
  int*  _flux_in_end_t;
  int*  _flux_out_t;
  int*  _flux_out_end_t;
  int  _constr;
  int  _id;
  siblings_t *_child;

  MotifTree( 
	    int  *constraint,
	    siblings_t *child ) {
    _col          = new int[1]; // used also as motif id pointer
    _col_size     = 1;
    _col_end      = _col + _col_size;
    *_col = 0;
    _row      = 0;
    _row_end  = 0;
    _row_size = 0;
    _constraint.push_back( constraint );
    _id = 0;
    _child = child;
  };

  MotifTree(  const int *col, 
	      int  col_size,
	      int  *constraint,
	      int  cstr,

	      siblings_t *child = 0,
	      int id = 0 ) {

    new int[col_size];
    _col          = new int[col_size]; // used also as motif id pointer
    _col_size     = col_size;
    _col_end      = _col + _col_size;
    copy( col, col + col_size, _col );
    _constr = cstr;

    _row      = 0;
    _row_end  = 0;
    _row_size = 0;
    _constraint.push_back( constraint );
    _id = id;
    _child = child;

  };

  MotifTree(  const int *col, 
	      int  col_size,
	      const int *row,
	      int  row_size,
	      int constraint,
	      siblings_t *child) {
    _col          = new int[col_size]; // used also as motif id pointer
    _col_size     = col_size;
    _col_end      = _col + _col_size;
    copy( col, col + col_size, _col );
    _row          = new int[row_size]; // used also as motif id pointer
    _row_size     = row_size;
    _row_end      = _row + _row_size;
    copy( row, row + row_size, _row );
    _child = child;
    _id = 0;
  };

  ~MotifTree( ) { 
    if( _col ) delete [] _col;
    if( _row ) delete [] _row;
  };


  void inline addLink( siblings_t *child ) { _child = child; };  
};

void printTList( const char *ptr, 
		const int *begin, const int *end, 
		 bool endl_ );

void printTree(siblings_t *current_sibs, string &indent, int motif_size );


void addMotif( siblings_t *root, siblings_t *current_sibs, 
	       SparseAdjacency &motif, int k );

void addMotifConstraints( siblings_t *root, siblings_t *current_sibs, 
			  SparseAdjacency &motif, int *constraint, int k );
void optimizeConstraints( siblings_t *current_sibs, int k );

# endif
