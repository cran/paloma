//  combinatorics.h
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

//! \file    combinatorics.h
//! \author  Gilles Grasseau
//! \brief   Implement tools from combinatorics
//! \version 0.1
//! \date    September 2009

# ifndef _S_COMBINATORICS_H_
# define _S_COMBINATORICS_H_

# include <iostream>
# include <list>
# include <algorithm>
# include <iterator>
# include <string.h>
# include <stdint.h>


//not used # include "util.h"

using namespace std;

// TODO : use template for config

//! \brief Use to generate configurations greater than 2^64 
//   Provide tools to handle configurations of  _StateVectorLength
//   with _NumberOfStates states.
class ConfigurationIterator{
  int  _NumberOfStates;
  int  _StateVectorLength;
  int *_StateVector;
  bool _end;

 public:

  void initialize( int NumberOfStates, int StateVectorLength ) {
    _NumberOfStates    = NumberOfStates;
    _StateVectorLength = StateVectorLength;
    _StateVector       = new int [ _StateVectorLength ];
    _end               = false;

    restart();
  }

  ConfigurationIterator( );
  
  ConfigurationIterator( int NumberOfStates, int StateVectorLength ) {
    initialize( NumberOfStates, StateVectorLength );
  }

  ~ConfigurationIterator( ) {
    delete [] _StateVector;
    _StateVector = 0;
  }

  const int *restart() {

    for (int i=0; i < _StateVectorLength; i++) { _StateVector[i] = 0 ;}; 
    _end = false;
    return( _StateVector );
  }
 

  ConfigurationIterator& operator++() {

    bool hold = !(_end);
    for ( int i = 0; (i < _StateVectorLength) && (hold); i++) {
      _StateVector[i]++;
      if( _StateVector[i] == _NumberOfStates ) {
	_StateVector[i] = 0;
	hold = true;
      } else {
	hold = false;
      }
    }
    _end = hold;
    return *this;
  }

  bool isLastConfiguration() { return _end; }
  bool iterate() { operator++(); return !(_end); }
    
};

//! \brief use to generate combinaison to choose k elements in a set S 
//! with card(S) = n, generate (n k) combinaisons.  
class CombinationIterator {

 protected:
  int _IndexesLength;
  int *_Indexes;
  int _SetLength;

  bool _last;
  bool _end;

 public:

  const int *restart() {

    for (int i=0; i < _IndexesLength; i++) { _Indexes[i] = i ;}; 
    _end = false;
    return( _Indexes );
  }

  const int *initialize( int n, int k ) {
    _IndexesLength = k;
    _SetLength   = n;
    _Indexes = new int [_IndexesLength];
    _last = false;
    _end = false;
    restart();
    return _Indexes;

  }

  CombinationIterator() { _Indexes = 0; };

  CombinationIterator( int n, int k ) {
    initialize( n, k ); 
  }

  CombinationIterator& operator++(void) {
    bool hold = false;
    for( int k = 0; k < (_IndexesLength); k++) {
      if ( _Indexes[ _IndexesLength - k - 1 ] != ( _SetLength - k -1)) { 
	_Indexes[ _IndexesLength -k -1 ]++;
	if( hold ) {
	  for ( int j = 0; j < k; j++) {
	    _Indexes[ _IndexesLength -k + j ] = 
	      _Indexes[ _IndexesLength -k + j - 1  ] + 1;
	  }
	  hold = false;
	}
	return *this;
      } else {
	hold = true;
      }
    }

    if( hold ) 
      _end = true;
     
    return *this;
   }

  bool isEnd() { return _end;} 

  bool iterate() { operator++(); return !( isEnd() ); }

  ~CombinationIterator() {
    if( _Indexes != 0 ) 
      delete [] _Indexes ;
  }

};

 
// n,k must be positive
 
int64_t factorial( int64_t n, int64_t k = 1 );


double BinomCoef( int n , int k);

double MultinomCoef( const int64_t n,  const int64_t *k );

double MultinomCoef( const int64_t n,  const int64_t n1,  
		     const int64_t n2, const int64_t n3);

// Test : som( ni ) = n

double MultinomCoef( const int64_t n,  const int64_t n1,  
		     const int64_t n2, const int64_t n3, 
		     const int64_t n4 ) ;



#endif
