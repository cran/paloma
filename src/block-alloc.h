//  block-alloc.h
// 
//  Copyright (C) 2010 Laboratoire Statistique & Genome
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

//! \file    block-alloc.h
//! \author  Gilles Grasseau
//! \brief   Storage allocated by block 
//! \version 0.1
//! \date    March 2010

// TODO :
//   Copy constructor

# include <iostream>

# ifndef _BLOCK_ALLOC_H_
# define _BLOCK_ALLOC_H_

using namespace std;

//! \class ??? SparseAdjacency 
//! \brief ???
//! \brief Allocate array by memory block.
//! to avoid numerous small allocation.
//! Small arrays of the same size "chunk_size" 
//! are allocate by larger blocks of size 
//! "chunk_size" * "chunk_nbr" (See constructor).
//! 


template <typename T>
class BlockAllocation {

 protected:

  // Block allocation size
  int _block_size;

  // _block_size is a multiple of _chunk_size. See constructor
  int _chunk_size;

  // Block structure to handle data
  struct block_t;
  struct block_t { 
    block_t *_next_block;
    T* _data;
    T* _data_end;
  };
  
  // Block structure iterator - state of allocation 
  block_t  *_block_start;
  block_t *_it_block; 

  // Block structure iterator - to scan the data 
  block_t *_itr_block; 

  // Data iterator inside one block (the current block)
  T* _it;
  T* _it_end;
  
  block_t inline *_allocateBlock() {
    block_t *block = new block_t;
    block -> _data = new T[_block_size];
    block -> _next_block = 0;
    _it = block -> _data;
    _it_end = _it + _block_size;
    block -> _data_end = _it_end;
    return  block ;
  }

  void inline _deallocate(block_t *block) {
    if( block -> _next_block ) {
      _deallocate( block -> _next_block );
    }
    delete [] block ->_data;
    delete block;
  }

 public:
  //!  \brief Constructor of an array with 
  //!  chunk_size * chunk_nb elements
  //!  \param chunk_size is the constant array size to allocate.
  //!  \param chunk_nbr number of chunks to allocate per block.
  BlockAllocation( int chunk_size=1, int chunk_nbr=4096) {
    _chunk_size = chunk_size;
    _block_size  = chunk_size * chunk_nbr;
    _block_start = _allocateBlock();
    _it_block = _block_start;
  }

  ~BlockAllocation() {
    _deallocate( _block_start );
  }

  //! \brief Store data with a constant size
  //!        This size chunk_size have been specified in the constructor
  //!        WARNING : The size must be equal to with those specified to 
  //!        the constructor
  //! \param data array to store
  void inline storeData( const T *data ) {
    if( _it == _it_end ) {
      _it_block -> _next_block = _allocateBlock();
      // Current block become the new allocated
      _it_block = _it_block -> _next_block;
    }
    const T *it_src = data;
    T *it_end = _it + _chunk_size;
    for( ; _it < it_end; ) {
      *_it++ = *it_src++;
    }
  }

  //! \brief Prepare the block allocation to be iterate (read/scanned).
  //! New block allocations are not possible (future enhancement)
  bool inline setBlockBegin() {
    // WARNING : set the current iterator it_ to the last block (data) end
    _it_block -> _data_end = _it;
    _itr_block = _block_start;
    return ( _block_start != 0 ); 
  }

 //! \brief Iterate the block
  bool inline iterateBlock( ) {
    bool exist;
    if( _itr_block -> _next_block ) {
      _itr_block = _itr_block -> _next_block;
      exist = true;
    } else { 
      exist = false;
    }
    return exist; 
  }

  //! \brief Get starting data of the current block
  T inline *getBlockDataBegin() {
    return (_itr_block -> _data) ;
  }
  //! \brief Get last data data of the current block
  T inline *getBlockDataEnd() {
    return (_itr_block -> _data_end) ;
  }
  //! \brief Print
  void print() {
    for( bool cont = setBlockBegin(); cont ; cont = iterateBlock( )){
      for( int *it = getBlockDataBegin(), *it_end =  getBlockDataEnd(); 
	   it != it_end; ) {
	cout << *it++ << " " << endl;
      } 
      cout << endl;
    } 
  }
};

# endif
