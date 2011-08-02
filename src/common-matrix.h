//  common-matrix.h
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
//

//! \file    common-matrix.h
//! \author  Gilles Grasseau
//! \brief   Common definitions for matrix modules.
//! \version 0.1
//! \date    November 2009

#ifndef _COMMON_MATRIX_H_
#define _COMMON_MATRIX_H_

//! \cond
#define EQUAL_DIM( X ) ( (NRow_ == (X).NRow_) && (NCol_ == (X).NCol_) ) 

#define EQUAL_DIM_MAT( X, Y ) ( ((X).getRow() == (Y).getRow()) \
                             && ((X).getCol() == (Y).getCol()) )

#define EQUAL_NCOL_ROW( X ) ( (NCol_ == (X).NRow_) )  
#define EQUAL_COL( X ) ( (NCol_ == (X).NCol_) )  
#define EQUAL_ROW( X ) ( (NRow_ == (X).NRow_) )  

//! \endcond

static const char *not_same_dims   = "Matrix don't have the same dimensions";
static const char *not_dims_match  = "Matrix dimensions don't match";
static const char *not_scalar      = "Matrix can't be transformed to scalar";
static const char *out_of_dims     = "Index out of matrix dimensions";
static const char *not_square_matrix = "Not a squared Matrix";
static const char *bad_row_size = "Bad row number";
static const char *bad_column_size = "Bad column number";
static const char *not_implemented = "Not yet implemented";

#endif
