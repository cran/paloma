//  combinatorics.cpp
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

//! \file    combinatorics.cpp
//! \author  Gilles Grasseau
//! \brief   Implement tools from combinatorics
//! \version 0.1
//! \date    September 2009

# include <iostream>
# include <list>
# include <algorithm>
# include <iterator>
# include <string.h>
# include <stdint.h>


// not used # include "util.h"
# include "combinatorics.h"

using namespace std;

// TODO : Template
//! \brief Compute (n)!/(k)!
//! Factorial n is obtained with k=1
//! param n
//! param k
//! return return (n)!/(k)!

#ifdef USE_INTERNAL_C_FCT

static int64_t _binomcoef[11][11]   = { 
  1,  0,  0,   0,   0,   0,   0,   0,  0,  0, 0 ,
  1,  1,  0,   0,   0,   0,   0,   0,  0,  0, 0 ,
  1,  2,  1,   0,   0,   0,   0,   0,  0,  0, 0 , 
  1,  3,  3,   1,   0,   0,   0,   0,  0,  0, 0 , 
  1,  4,  6,   4,   1,   0,   0,   0,  0,  0, 0 , 
  1,  5, 10,  10,   5,   1,   0,   0,  0,  0, 0 ,
  1,  6, 15,  20,  15,   6,   1,   0,  0,  0, 0 ,
  1,  7, 21,  35,  35,  21,   7,   1,  0,  0, 0 ,
  1,  8, 28,  56,  70,  56,  28,   8,  1,  0, 0 ,
  1,  9, 36,  84, 126, 126,  84,  36,  9,  1, 0 ,
  1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1 
};
 

// n,k must be positive
int64_t BinomCoef_( int n , int k) {

  if ( k > n) {
    return 0;
  } else {
    if ( k == 0 )
      return 1L;
    else if ( k == 1 )
      return static_cast<int64_t> (n); 
    else if ( n < 11 ) {
      return _binomcoef[n][k];
    } else {
      // Vandermonde identity
      if ( k > (n/2)) 
	k = n - k;
      int m1 = n/2;
      int m2 = n - m1;
      int64_t res = 0;
      if ( m1 != m2 ) {
	// TODO : Suppress this form 
	/*
	for ( int i=0; i < (k+1) ; i++) {
	  res += (BinomCoef( m1, i ) * BinomCoef( m2, k-i ));
	}
	*/
	// Pascal triangle
	res = BinomCoef_( 2*m1, k) +  BinomCoef_( 2*m1, k-1);
      } else {
	int k_end ;
	if ( k % 2) {
	  // Odd
	  k_end = (k+1)/2; 
	} else {
	  // even
	  k_end = k/2;
	  res = BinomCoef_( m1, k_end );
	  res = res * res;
	}
	
	for ( int i=0; i < (k_end) ; i++) {
	  res += 2 * (BinomCoef_( m1, i ) * BinomCoef_( m1, k-i ));
	}
      }
      // cout << n << ", " << k << ": " << res << endl;
      
      return res;
    }
  }
}

 

double BinomCoef( int n , int k) {
  return static_cast<double> ( BinomCoef_( n , k) );
}

double MultinomCoef_( const int64_t n,  const int64_t *k ) {

  int64_t max = -1;
  int i_max = -1;
  for( int i=0; k[i] != -1; i++) {
    if( k[i] > max ) {
      i_max = i;
      max =  k[i];
    }
  }
  double res = static_cast<double> ( factorial( n, max ));
  int64_t coef = 1;
  for( int i=0; k[i] != -1; i++) {
    if( i != i_max ) {
      coef = coef *  factorial( k[i], 1); 
    }
  }
  res = res / coef;
  return res;
}

// ???
double BinomCoefOld( const int64_t n, 
		  const int64_t n1 ) {

  int64_t parameters[3];
  parameters[0] = n1;
  parameters[1] = n - n1;
  parameters[2] = -1;
  
  return MultinomCoef_( n, parameters );
  
}

double MultinomCoef( const int64_t n,  const int64_t n1,  
			    const int64_t n2, const int64_t n3) {

  int64_t parameters[4];
  parameters[0] = n1;
  parameters[1] = n2;
  parameters[2] = n3;
  parameters[3] = -1;
  
  return MultinomCoef_( n, parameters );
  
}

// Test : som( ni ) = n
double MultinomCoef( const int64_t n,  const int64_t n1,  
			    const int64_t n2, const int64_t n3, 
			    const int64_t n4 ) {

  int64_t parameters[5];
  parameters[0] = n1;
  parameters[1] = n2;
  parameters[2] = n3;
  parameters[3] = n4;
  parameters[4] = -1;
  
  return MultinomCoef_( n, parameters );
  
}

#else
#include <R.h>
#include <Rmath.h>

double BinomCoef( int n , int k) {
  double res = exp( lgammafn( n+1) -  lgammafn( k+1) -  lgammafn( n-k+1 ) );
  return res;
}


double MultinomCoef( const int64_t n,  const int64_t n1,  
		     const int64_t n2, const int64_t n3) {
  double res = exp( lgammafn( n+1) -  lgammafn( n1+1) -  lgammafn( n2+1) 
		    - lgammafn( n3+1) );
  return res;
}
// Test : som( ni ) = n

double MultinomCoef( const int64_t n,  const int64_t n1,  
		     const int64_t n2, const int64_t n3, 
		     const int64_t n4 ) {
  double res = exp( lgammafn( n+1) -  lgammafn( n1+1) -  lgammafn( n2+1) 
		    - lgammafn( n3+1)- lgammafn( n4+1) );
  return res;
}

#endif


int64_t factorial( int64_t n, int64_t k) {
  
  int64_t res = (n >= (k+1)) ? (k + 1) : 1 ; 
  for ( int64_t i = k+2; i <= n; i++) res *= i;
  return res;
} 
