//  interface_R_cpp.cpp
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

//! \file    interface_R_cpp.cpp
//! \author  Gilles Grasseau 
//! \brief   R/C interface
//! \version 0.1
//! \date    December 2010

# include "sparse-adjacency.h"
# include "motif-search.h"
# include "graphical-models.h"

# include <R.h>
# include <Rdefines.h>

# include "interface_R_cpp.h"
# include "concentration.h"

SEXP getCanonic_C( SEXP S_M, SEXP S_M_nnodes, 
		   SEXP S_Directed) {

  // S_M adjacency matrix of the motif 
  // S_M_nodes node number in the motif
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  PROTECT( S_M        = AS_INTEGER(S_M)    ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  int *M          = INTEGER_POINTER( S_M );
  int M_nnodes     = INTEGER_POINTER( S_M_nnodes )[0];
  bool Directed   = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  SEXP S_Res;  
  PROTECT( S_Res = allocVector( INTSXP, M_nnodes*M_nnodes ) ); n_prot++; 

  int *canonic = INTEGER_POINTER( S_Res );

  int *perm = new int[ M_nnodes ];

  getCanonic( M , M_nnodes, Directed, canonic, perm );

  UNPROTECT( n_prot );

  delete [] perm;
  return S_Res;
}

SEXP get_Exceptional_M_Mp_C( SEXP S_G, SEXP S_G_nnodes, 
			 SEXP S_Min_nnodes, SEXP S_Max_nnodes,
			 SEXP S_NodeToClass,
			 SEXP S_Pi,
			 SEXP S_NbrClasses,
			 SEXP S_PValue,
			 SEXP S_Directed ) {
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M_nnodes number of nodes of M
  // S_NodeToClass : classe for each node
  // S_Pi  : Connectivity probabilities between classes.
  // S_NbrClasses  : number of classes
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_Min_nnodes  = AS_INTEGER(S_Min_nnodes) ); n_prot++;
  PROTECT( S_Max_nnodes  = AS_INTEGER(S_Max_nnodes) ); n_prot++;
  PROTECT( S_NodeToClass = AS_INTEGER(S_NodeToClass) ); n_prot++;
  PROTECT( S_Pi          = AS_NUMERIC(S_Pi)       ); n_prot++;
  PROTECT( S_NbrClasses  = AS_INTEGER(S_NbrClasses) ); n_prot++;
  PROTECT( S_PValue      = AS_NUMERIC(S_PValue)       ); n_prot++;
  PROTECT( S_Directed    = AS_INTEGER(S_Directed) ); n_prot++;

  // size_t n_    = LENGTH( S_G );
  
  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int M_Min_nnodes = INTEGER_POINTER( S_Min_nnodes )[0];
  int M_Max_nnodes = INTEGER_POINTER( S_Max_nnodes )[0];
  int *NodeToClass = INTEGER_POINTER(S_NodeToClass);
  double *Pi       = NUMERIC_POINTER( S_Pi );
  int NbrClasses   = INTEGER_POINTER(S_NbrClasses) [0]; 
  double PValue    = NUMERIC_POINTER( S_PValue )[0];
  bool Directed    = ( INTEGER_POINTER( S_Directed )[0] != 0 );



  SEXP S_ReturnList = getExceptional( G_edges, G_nnodes, 
				      M_Min_nnodes, M_Max_nnodes, 
				      NodeToClass, Pi, NbrClasses, 
				      PValue, Directed,
				      n_prot);

  UNPROTECT( n_prot );

  return S_ReturnList ;

}

SEXP countAllExactMotifs_C( SEXP S_G, SEXP S_G_nnodes, 
			    SEXP S_M_nnodes,
			    SEXP S_Directed ) 
{
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M_nnodes motif number of nodes 
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  // size_t n_    = LENGTH( S_G );

  
  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int M_nnodes = INTEGER_POINTER( S_M_nnodes )[0];
  bool Directed = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  int k = M_nnodes;

  //
  // Get all occurrences of size k = M_nnodes
  //

  FindMotif find( StoreMode, 0);

  vector<int*> *list;
  {
    SparseAdjacency G( G_edges, G_nnodes, true );
    // G.print( "G");

    list = find.findAllMotifs( G, k ); 
  }

  SEXP S_ReturnList;
  SparseAdjacency H( G_edges, G_nnodes, ! Directed  );

  // H.print( "H");

  S_ReturnList = count_m( H, k, list, n_prot );

  UNPROTECT( n_prot );

  return( S_ReturnList );

}

SEXP computePseudoMixNetMean_C( SEXP S_N, SEXP S_nclass,
			    SEXP S_Alpha, SEXP S_Pi, 
			    SEXP S_M, SEXP S_M_nedges,
			    SEXP S_Directed ) {
  // S_N Number of nodes. 
  // S_nclass Number of classes. 
  // S_Alpha number of nodes per classes
  // S_Pi  Connectivity probabilities between classes.
  // S_M adjacency matrix of the motif 
  // S_M_nedges edge number in the motif
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_N        = AS_INTEGER(S_N)        ); n_prot++;
  PROTECT( S_nclass   = AS_INTEGER(S_nclass)   ); n_prot++;
  PROTECT( S_Alpha    = AS_INTEGER(S_Alpha)    ); n_prot++;
  PROTECT( S_Pi       = AS_NUMERIC(S_Pi)       ); n_prot++;
  PROTECT( S_M        = AS_INTEGER(S_M)        ); n_prot++;
  PROTECT( S_M_nedges = AS_INTEGER(S_M_nedges) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  int64_t N     = INTEGER_POINTER( S_N )[0];
  int    nclass   = INTEGER_POINTER( S_nclass )[0];
  int *Alpha      = INTEGER_POINTER( S_Alpha );
  double *Pi      = NUMERIC_POINTER( S_Pi );
  int *M          = INTEGER_POINTER( S_M );
  size_t M_nedges = INTEGER_POINTER( S_M_nedges )[0];
  bool Directed   = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  // Allocate storage for results
  SEXP S_Res;
  PROTECT( S_Res = NEW_NUMERIC( 1 ) ); n_prot++;
  double *Res = NUMERIC_POINTER( S_Res );

  // Transform to size_t
  int max = 0; 
  int *M_t = new int[M_nedges*2+2];
  for ( size_t i=0; i < M_nedges*2; i++) {
    M_t [i] = M[i];
   max = MAX( max,  M_t[i] );
  }
  M_t [M_nedges*2]   = NOT_INDEX;
  M_t [M_nedges*2+1] = NOT_INDEX;
  int nnodes = max + 1;

  SparseAdjacency motif( M_t, nnodes, !Directed );

  // motif.print("motif :");
  NonRedondantPermutation NRPMotif( motif );

  PseudoMixnet pmixnet(N, nclass, Alpha, Pi, Directed);

  Res[0] = pmixnet.computeMean( motif, NRPMotif );
  // Res[1] = mixnet.computeMoment2( motif, NRPMotif );

  UNPROTECT(n_prot);

  delete [] M_t;

  return(S_Res);

}

SEXP get_particular_M_Mp_Occurrences_C( SEXP S_G, SEXP S_G_nnodes, 
				     SEXP S_M_nnodes,
				     SEXP S_CanonicFilter,
				     SEXP S_DelClassFilter,
				     SEXP S_Directed ) {
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M_nnodes number of nodes of M
  // S_CanonicFilter canonic form to filter
  // S_DelClassFilter deletion class to filter (one index)
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_CanonicFilter = AS_INTEGER(S_CanonicFilter) ); n_prot++;
  PROTECT( S_DelClassFilter = AS_INTEGER(S_DelClassFilter) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  // size_t n_    = LENGTH( S_G );

  
  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int M_nnodes = INTEGER_POINTER( S_M_nnodes )[0];
  int *CanonicFilter    = INTEGER_POINTER( S_CanonicFilter );
  int DelClassFilter = INTEGER_POINTER( S_DelClassFilter )[0];
  bool Directed = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  int k = M_nnodes;

  //
  //  Find all motif with size k
  //
  FindMotif find( StoreMode, 0);

  vector<int*> *list;
  {
    SparseAdjacency G( G_edges, G_nnodes, true );
    // G.print( "G");

    list = find.findAllMotifs( G, k ); 
  }

  SEXP S_ReturnList, S_Struct;
  SEXP S_Result;
  SparseAdjacency H( G_edges, G_nnodes, ! Directed  );

  // H.print( "H");

  // Compute canonical form
  int *perm = new int[ M_nnodes ];
  int *canonic = new int[ M_nnodes * M_nnodes ];

  getCanonic( CanonicFilter, M_nnodes, Directed, canonic, perm );

  S_Result = sort_m_mp_u( H,  G_nnodes, k, list, 
			  0, 0, 0, 
			  canonic, DelClassFilter, 
			  0.0,                     
			  n_prot, true );

  // k value
  SEXP S_k;
  PROTECT( S_k = allocVector( INTSXP, 1) ); n_prot++; 	  
  int *p_k = INTEGER_POINTER( S_k );
  p_k[0] = k;

  // Return structure
  PROTECT( S_Struct = allocVector( VECSXP, 2 ));  n_prot++;;
  SET_VECTOR_ELT( S_Struct, 0, S_k ); 
  SET_VECTOR_ELT( S_Struct, 1, S_Result ); 

  // Return list
  PROTECT( S_ReturnList = allocVector( VECSXP, 1 )); n_prot++;
  SET_VECTOR_ELT( S_ReturnList, 0, S_Struct );

  UNPROTECT( n_prot );

  delete [] perm ;
  delete [] canonic ;
  return( S_ReturnList );

}

// Get all occurrences
SEXP get_M_Mp_Occurrences_C( SEXP S_G, SEXP S_G_nnodes, 
			    SEXP S_M_nnodes,
			    SEXP S_Directed ) 
{
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M_nnodes number of nodes of M
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_Directed = AS_INTEGER(S_Directed) ); n_prot++;

  // size_t n_    = LENGTH( S_G );

  
  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int M_nnodes = INTEGER_POINTER( S_M_nnodes )[0];
  bool Directed = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  int k = M_nnodes;

  //
  // Get all occurrences of size k = M_nnodes
  //

  //
  //  Find all motif with size k
  //
  FindMotif find( StoreMode, 0);

  vector<int*> *list;
  {
    SparseAdjacency G( G_edges, G_nnodes, true );
    // G.print( "G");

    list = find.findAllMotifs( G, k ); 
  }

  SEXP S_ReturnList;
  SparseAdjacency H( G_edges, G_nnodes, ! Directed  );

  // H.print( "H");

  S_ReturnList = sort_m_mp_u( H,  G_nnodes, k, list, 
			      0, 0, 0, 
			      0, -1,
			      0.0,
			      n_prot, true );
  UNPROTECT( n_prot );

  return( S_ReturnList );

}

SEXP get_M_Mp_PValues_C( SEXP S_G, SEXP S_G_nnodes, 
			 SEXP S_M_nnodes,
			 SEXP S_NodeToClass,
			 SEXP S_Pi,
			 SEXP S_NbrClasses,
			 SEXP S_PValue,
			 SEXP S_Directed ) {
  // S_G edge list of the network G
  // S_G_nnodes number of nodes of G
  // S_M_nnodes number of nodes of M
  // S_NodeToClass : classe for each node
  // S_Pi  : Connectivity probabilities between classes.
  // S_NbrClasses  : number of classes
  // S_PValue : pvalue
  // S_Directed non-zero integer if the graph is directed

  int n_prot = 0;

  // Arguments
  PROTECT( S_G        = AS_INTEGER(S_G)    ); n_prot++;
  PROTECT( S_G_nnodes = AS_INTEGER(S_G_nnodes) ); n_prot++;
  PROTECT( S_M_nnodes = AS_INTEGER(S_M_nnodes) ); n_prot++;
  PROTECT( S_NodeToClass = AS_INTEGER(S_NodeToClass) ); n_prot++;
  PROTECT( S_Pi          = AS_NUMERIC(S_Pi)       ); n_prot++;
  PROTECT( S_NbrClasses  = AS_INTEGER(S_NbrClasses) ); n_prot++;
  PROTECT( S_PValue      = AS_NUMERIC(S_PValue)       ); n_prot++;
  PROTECT( S_Directed    = AS_INTEGER(S_Directed) ); n_prot++;

  // size_t n_    = LENGTH( S_G );
  
  // Passing argument by addresses
  int *G_edges = INTEGER_POINTER( S_G );
  int G_nnodes = INTEGER_POINTER( S_G_nnodes )[0];
  int M_nnodes = INTEGER_POINTER( S_M_nnodes )[0];
  int *NodeToClass = INTEGER_POINTER(S_NodeToClass);
  double *Pi       = NUMERIC_POINTER( S_Pi );
  int NbrClasses   = INTEGER_POINTER(S_NbrClasses) [0];
  double pvalue    = NUMERIC_POINTER( S_PValue) [0];
  bool Directed    = ( INTEGER_POINTER( S_Directed )[0] != 0 );

  int k = M_nnodes;

  //
  // Get all occurrences of size k = M_nnodes
  //

  //
  //  Find all motif with size k
  //
  FindMotif find( StoreMode, 0);

  vector<int*> *list;
  {
    SparseAdjacency G( G_edges, G_nnodes, true );
    // G.print( "G");

    list = find.findAllMotifs( G, k ); 
  }

  SEXP S_ReturnList;
  SparseAdjacency H( G_edges, G_nnodes, ! Directed  );

  // H.print( "H");

  S_ReturnList = sort_m_mp_u( H,  G_nnodes, k, list, Pi, 
			      NodeToClass, NbrClasses,
			      0, -1,
			      pvalue,
			      n_prot, false);

  UNPROTECT( n_prot );

  return( S_ReturnList );

}



