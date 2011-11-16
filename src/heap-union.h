// 
# ifndef _HEAP_UNION_H_
# define _HEAP_UNION_H_

# include <iostream>
# include "util.h"
# include "sparse-adjacency.h"

using namespace std;

class HeapUnion {

 protected :

  const int **_set;
  const int **_set_end;
  int *_index;
  int _nb_set;
  int _nb_active_set;
  bool   _not_finished;

  // Location of last element picked up from the heap
  const int *found_in_set;
  int index_in_set; 

  void _sort( int left, int right);

 public:
  // No empty set
  HeapUnion(const int **set, const int **set_end, int nb_set ) {
    _set = set;
    _set_end = set_end;
    _nb_set = nb_set;
    _nb_active_set = nb_set;
    _index = new int[ _nb_set ];
    for( int *it=_index, *it_end=&_index[nb_set]; it != it_end; it++) {
      *it = 0;
    }
    if ( nb_set > 1) 
      _sort( 0, nb_set - 1 );
    _not_finished = true;
    
  }

  ~HeapUnion() { delete [] _index ; };

  void restart( void  ) {

    _nb_active_set = _nb_set;
    for( int *it=_index, *it_end=&_index[_nb_set]; it != it_end; it++) {
      *it = 0;
    }
    /// ???
    // if ( _nb_set > 1) 
	 //  _sort( 0, _nb_set - 1 );
    _not_finished = true;
    
  }

  void status( );

  int nextElement();

  const int inline *getElementLocation( ) { return found_in_set; }

  int inline getElementLocation(const int **initial) {
    int k = 0;
    const int **it     = initial;
    const int **it_end = initial + _nb_set;

    for ( ; (it != it_end) && (*it != found_in_set); it++) {
      k++;
    }
    if( it != it_end )
      return k;
    else
      // GG
      // exit( 0) 
      return( -1 );
  }

  int inline getElementIndex( ) { return index_in_set; }

};

void HeapUnion::_sort(int left_s, int right_s) {
  

  int left = left_s, right = right_s;
  int i = left, j = right;
  const int *tmp;
  int var;
  int pivot = _set[(left + right) / 2][0];
  
  /* partition */
  while (i <= j) {
    while (_set[i][0] < pivot)
      i++;
    while (_set[j][0] > pivot)
      j--;
    if (i <= j) {
      // swap
      tmp = _set[i];
      _set[i] = _set[j];
      _set[j] = tmp;
      tmp = _set_end[i];
      _set_end[i] = _set_end[j];
      _set_end[j] = tmp;
      var = _index[i];
      _index[i] = _index[j];
      _index[j] = var;
      i++;
      j--;
    }
  }
  
  /* recursion */
  if (left < j)
    _sort(left, j);
  if (i < right)
    _sort(i, right);
  // status();
}

int HeapUnion::nextElement() {

  int next_value;


  if ( _not_finished ) {
    next_value = _set[ 0 ] [ _index[0] ];
    found_in_set = _set[0];
    index_in_set = _index[0];

  } else {
    return NOT_INDEX;
  }

  // ??? To remove bool go_on = true;

  for( bool go_on = true ; go_on ; ) {

  // next element
  _index[ 0 ]++;

  int heap_bottom;
  int heap_top   ;
  int heap_mid   ;
  // int last_index_inf = 0;
  // int last_index_sup = 0;
  bool cont = true;
  bool val = false;

  // Is the last element 
  if( _set_end[ 0 ] == ( &_set[0] [ _index[ 0 ] ] ) ) {
    heap_top = _nb_active_set - 1;
    _nb_active_set--;

    
    if ( _nb_active_set == 0 ) {
      _not_finished = false;
    }
  } else {

    // Sort the new value in the heap 
    // Find the place of the new element
    heap_bottom = 0;
    heap_top    = _nb_active_set - 1;
    // last_index_inf = 0;
    // last_index_sup = 0;
    int new_value   = _set[ 0 ] [ _index[0] ];
    
    cont = (heap_bottom != heap_top);
    for( ; cont; ) {
      if ( _nb_active_set > 2 ) {
	// If no element in this set
	heap_mid = (heap_top + heap_bottom)/2;
	if( ( val = (new_value <= _set[ heap_mid ] [ _index[ heap_mid ] ] ) ) ) {
	  heap_top    = heap_mid;
	  // last_index_inf = _set[ heap_mid ] [ _index[ heap_mid ] ];
	} else {
	  heap_bottom = heap_mid+1;
	  // last_index_sup = _set[ heap_mid ] [ _index[ heap_mid ] ];
	}
	// cont = ( (heap_top != heap_mid ) &&  (heap_bottom != heap_mid ) );
	cont = ( (heap_top != (heap_bottom ) ) );
      } else {
	if( new_value <= _set[ 1 ] [ _index[ 1 ] ] ) {
	  heap_top = 0;
	  heap_bottom = 0;
	  val = false;
	} else {
	  heap_top    = 1;
	  heap_bottom = 1;
	  val = false;
	}
	cont = false;
      }
    } 
  }
  if(  (heap_top == (_nb_active_set -1) ) && val ) {
    heap_top = heap_top - 1;
    
  } 
  //  if ( (heap_top == (_nb_active_set -1)) && ( new_value < ) )
    
  //  Move heap
  //  save 
  if( heap_top > 0 ) {
  const int *save_set     = _set[0];
  const int *save_set_end = _set_end[0];
  int save_index    = _index[0];
  for( int i = 1; i < (heap_top + 1); i++){
    _set[i-1]     = _set[i];
    _set_end[i-1] = _set_end[i];
    _index[i-1]   = _index[i];
  } 
  _set[heap_top]     = save_set;
  _set_end[heap_top] = save_set_end;
  _index[heap_top]   = save_index;
  }
  // _index[0] must be a valid index and points to next_value
  if(  _index[0] < ( _set_end[0] - _set[0] ) ) { 
    go_on =  ( _set[ 0 ] [ _index[0] ] == next_value );
  } else 
    go_on = false;
  }
  // status();
  return next_value;
}


void HeapUnion::status() {

#   ifdef MSG
  cout << "Set content" << endl;
  for ( int k =0; k < _nb_set; k++) {
    cout << "  " ;
    for( const int *it=_set[k], *it_end=_set_end[k]; it != it_end; it++) {
      cout << *it << " " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << "Indexes" << endl;
  cout << "  " ;
  for ( int k =0; k < _nb_set; k++) {
    cout << _index[k] << " (" << (_set_end[k] - _set[k]) << ") " ;
  }
  cout << endl;
  cout << "Active set : " << _nb_active_set << " , finished : " << !_not_finished << endl;
#   endif
}

# endif
