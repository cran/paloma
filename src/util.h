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

# ifndef _UTIL_H_
# define _UTIL_H_
# include <exception>
# include <iostream> 
# include <fstream>
# include <sstream>
# include <string>
#include <sys/time.h>

using namespace std;

//! \cond

//
//    Log tools
//

void addLogContext( const char * str );

void addLogTabulation( );

const char* getContext() ; 

void removeLogContext( const char * str );

void removeLogTabulation() ; 


//
//      Error Class
//

// Message Type
//  
//  0 - No message
//  1 - Error Message + exception
//  2 - 1 and Warning message 
//  3 - 2 and Info message
//  4 - 3 and Light Debug
//  5 - 4 and Debug

static int Msg_level = 1;

class ErrorClass: public std::exception
{
  const char* context;
  const char* message;
  const char* info;

 public :

  ErrorClass( const char* context_, const char* message_, const char* info_) {
    context = context_;
    message = message_;
    info    = info_;
  }  

  ErrorClass( const char* context_, const char* message_, long int info_) {
    context = context_;
    message = message_;

    std::stringstream str_io("");
    std::string str;
    str_io << info_;
    str_io >> str;
    info = str.data();

  }  
  const char *getContext() { return context; };
  const char *getMessage() { return message; };
  const char *getInfo()    { return info; };

  virtual const char* what() const throw()
  {
    std::string str(">>> Error in ");
    str = str + context + ": " + message + ", " + info;
    return str.data();
  }
};

// Utilities
# define MAX(a,b) ( ((a) > (b)) ? (a) : (b) ) 
# define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )

# define LOGMSG( a, s, b, c)  if ((a) <= Msg_level) {		\
    if ((a) != 1) { (s) << getContext() << (b) << (c) << std::endl;}	\
      else{ ErrorClass e(getContext(), (b), (c)) ; throw e; }		\
    } 

# define LOGEXCEPT(e) ( std::cout << std::endl << e.what() << std::endl )

# define LOGFLUSH     std::cout << std::flush

#ifdef DEBUG
#  define LOGInFct(s)  LOGMSG( 4, s, "Entering in : ",  __FUNCTION__ )
#  define LOGOutFct(s) LOGMSG( 4, s, "Living      : ",  __FUNCTION__ )
#else
#  define LOGInFct(s)
#  define LOGOutFct(s)
#endif

# define COUT std::cout << getContext()

//
//      Useful with size_t
//
#define SIZE_T_MAX ((size_t) -1)

template <class TypeAlloc_>
void setValue( TypeAlloc_ *ptr, TypeAlloc_ val, size_t size ) {
  for( TypeAlloc_ *it = ptr, *it_end = &ptr[size]; it != it_end; it++)
    *it = val;
}

template <class TypeAlloc_>
TypeAlloc_ *allocate( TypeAlloc_ val, size_t size) {
  TypeAlloc_ *ptr = new TypeAlloc_[ size ];
  setValue( ptr, val, size );
  return ptr;
}

void printList( const char *ptr, 
		const int *begin, const int *end, 
		const char*sep, bool endl_ );

void writeMatrix( ofstream &outfile, const char *ptr, 
		  const bool *begin, const bool *end, int stride, 
		  const char*sep, bool endl_ );

void printMatrix( const char *ptr, 
		  const bool *begin, const bool *end, int stride, 
		  const char*sep, bool endl_ );

void printList( const char *ptr, 
		const int *begin, const int *end, 
		bool endl_ );


void getTime_( int ID );
void getTime_null( int ID );
 
# ifdef TIMING
   extern const int n_count=10;
   extern double counter[2*n_count];
   extern struct timeval time_;
#  define getTime( ID ) (  getTime_( ID ) )

# else
#  define getTime( ID ) (  getTime_null( ID ) )

# endif

#endif


//! \endcond
