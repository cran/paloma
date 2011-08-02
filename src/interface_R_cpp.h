//  interface_R_cpp.h
//
//  Copyright (C) 2010 Gilles Grasseau
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

//! \file    interface_R_cpp.h
//! \author  Gilles Grasseau 
//! \brief   R/C interface
//! \version 0.1
//! \date    December 2010

# include <stdint.h>
# include "nauty_interface.h"

# include <R.h>
# include <Rdefines.h>


extern "C"  {
  SEXP getCanonic_C( SEXP S_M, SEXP S_M_nnodes, 
		     SEXP S_Directed);

  SEXP get_M_Mp_Occurrences_C( SEXP S_G, SEXP S_G_nnodes, 
			      SEXP S_M_nnodes,
			      SEXP S_Directed );

  SEXP get_M_Mp_PValues_C( SEXP S_G, SEXP S_G_nnodes, 
			   SEXP S_M_nnodes,
			   SEXP S_NodeToClass,
			   SEXP S_Pi,
			   SEXP S_NbrClasses,
			   SEXP S_PValue,
			   SEXP S_Directed );

  SEXP get_particular_M_Mp_Occurrences_C( SEXP S_G, SEXP S_G_nnodes, 
					 SEXP S_M_nnodes,
					 SEXP S_CanonicFilter,
					 SEXP S_DelClassFilter,
					 SEXP S_Directed ) ;

  SEXP get_Exceptional_M_Mp_C( SEXP S_G, SEXP S_G_nnodes, 
			   SEXP S_Min_nnodes, SEXP S_Max_nnodes,
			   SEXP S_NodeToClass,
			   SEXP S_Pi,
			   SEXP S_NbrClasses,
			   SEXP S_PValue,
			   SEXP S_Directed );

  SEXP countAllExactMotifs_C( SEXP S_G, SEXP S_G_nnodes, 
			      SEXP S_M_nnodes,
			      SEXP S_Directed );

  SEXP computePseudoMixNetMean_C( SEXP S_N, SEXP S_nclass,
				  SEXP S_Alpha, SEXP S_Pi, 
				  SEXP S_M, SEXP S_M_nedges,
				  SEXP S_Directed ) ;
}

