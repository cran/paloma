//  util.h
// 
//  Copyright (C) 2008 Laboratoire Statistique & Genome
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

//! \file    util.h
//! \author  Gilles Grasseau
//! \brief   Basic tools for matrix
//! \version 0.1
//! \date    June 2008

# include <sstream>
# include <string>
# include <string.h>

# include "util.h"

using namespace std;

//! \cond

//
//    Log tools
//

static string str_ctxt = "";

void addLogContext( const char * str ) {
  str_ctxt += str;
}

void addLogTabulation() {
  str_ctxt += "  ";
}

const char* getContext() {
  return str_ctxt.data();
}

void removeLogContext( const char * str ) {
  string::iterator it = str_ctxt.end();

  str_ctxt.erase ( it-strlen(str), it );
}

void removeLogTabulation() {
  string::iterator it = str_ctxt.end();

  str_ctxt.erase ( it-2, it );
}

void printList( const char *ptr, 
		const int *begin, const int *end, 
		const char*sep, bool endl_ ) {

  cerr << ptr;
  for( const int *it = begin; (it != end ); it++) {
      cerr << sep << *it;
  }
  if( endl_ )
    cerr << endl;
}


void writeMatrix( ofstream &outfile, const char *ptr, 
		  const bool *begin, const bool *end, int stride, 
		const char*sep, bool endl_ ) {

  outfile << ptr;
  const bool *row_end = begin + stride;
  for( const bool *it_row = begin; (it_row != row_end ); it_row++) {
    const bool *it_end = it_row + stride * stride;
    for( const bool *it = it_row; (it != it_end ); it+=stride) {
      outfile << sep << *it ;
    }
    outfile << endl;
  }
  if( endl_ )
    outfile << endl;
}

void printMatrix( const char *ptr, 
		  const bool *begin, const bool *end, int stride, 
		  const char*sep, bool endl_ ) {

  cout << ptr;
  const bool *row_end = begin + stride;
  for( const bool *it_row = begin; (it_row != row_end ); it_row++) {
    const bool *it_end = it_row + stride * stride;
    for( const bool *it = it_row; (it != it_end ); it+=stride) {
      cout << sep << *it ;
    }
    cout << endl;
  }
  if( endl_ )
    cout << endl;
}

void printList( const char *ptr, 
		const int *begin, const int *end, 
	bool endl_ ) {

  cerr << ptr;
  for( const int *it = begin; (it != end ); it++) {
      cerr << " " << *it;
  }
  if( endl_ )
    cerr << endl;
}

void dummy() {
  LOGMSG( 10, cout , " ", " ");
}

const int n_count=10;
double counter[2*n_count];
struct timeval time_;


/*
#  define getTime( ID ) ( gettimeofday(&time_,0); counter[(ID)] += ( double( time_.tv_sec*1000 ) + double( time_.tv_usec)/1000.0 ); )
# else
#  define getTime( ID ) ( ID; )
*/


void getTime_( int ID ) {

  gettimeofday(&time_,0); 
  counter[(ID)] += ( double( time_.tv_sec*1000 ) 
		     + double( time_.tv_usec)/1000.0 ); 
}

// Do nothing
void getTime_null( int ID ) {
}


//! \endcond

