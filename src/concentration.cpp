//  concentration.cpp
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

//! \file    concentration.cpp
//! \author  Gilles Grasseau 
//! \brief   Implement all routines for motif concentration
//! \version 0.1
//! \date    December 2010

// notes : R --max-ppsize 100000

# include <iostream>
# include <fstream>
# include <list>
# include <ctime>
# include <stdlib.h>
# include <math.h>

# include "sparse-adjacency.h"
# include "motif-search.h"
# include "nauty_interface.h"
# include "graphical-models.h"

# include <R.h>
# include <Rinternals.h> 
# include <Rdefines.h>

# include "concentration.h"

# define MSG 0

static double Pi        = acos(-1.0);
static double TwoPi     =  2*Pi;
static double FourPi    =  4*Pi;
static double SixteenPi = 16*Pi;


// Struture to handle the sort on the canonic form
typedef struct{
  int *canonic;
  int *occurrence;
} canonic_list_t;

// Structure to handle the 2nd sort on the occurrences U 
typedef struct{
  int *occurrence;
  int remove;
  int flag;
} m_mp_list_t;

static int sort_canonic_nbr_elmnt=0;
static int sort_occurences_nbr_elmnt = 0;

// Comparaison function to sort canonical forms
int canonic_order (const void * a, const void * b) {

  int *it1 = ((canonic_list_t*) (a)) -> canonic;
  int *it2 = ((canonic_list_t*) (b)) -> canonic;
  int *it1_end = it1 + sort_canonic_nbr_elmnt; 
  for ( ; it1 != it1_end; it1++, it2++){
    if( *it1 != *it2 ) {
      return *it1 - *it2;
    }
  }
  return (*--it1 - *--it2);
}

// Comparaison function to sort motif occurences
// with one deleted node
inline int occurrences_order (const void * a, const void * b) {

  int *it1 = ((m_mp_list_t*) (a)) -> occurrence;
  int *it2 = ((m_mp_list_t*) (b)) -> occurrence;
  int remove_col1 =  ((m_mp_list_t*) (a)) -> remove;
  int remove_col2 =  ((m_mp_list_t*) (b)) -> remove;
  int *it1_end = it1 + sort_occurences_nbr_elmnt; 
  int *it2_end = it2 + sort_occurences_nbr_elmnt; 
  int i1 = 0, i2 = 0;
  for ( ; it1 != it1_end; it1++, it2++, i1++, i2++){
    if( i1 == remove_col1 ) {
      it1++; i1++;
      if( it1 == it1_end )
       	return 0 ;
    }
    if( i2 == remove_col2 ) {
      it2++; i2++;
      if( it2 == it2_end )
	return 0 ;
    }
    if( *it1 != *it2 ) {
	return *it1 - *it2;
    }
  }
  return (*--it1 - *--it2);
}

double *AllocateDoubleSEXP( size_t size, SEXP &S_exp, int &n_prot ) {
  PROTECT( S_exp = allocVector( REALSXP, size) ); n_prot++; 	  
  double *out = NUMERIC_POINTER( S_exp );
  return out;
}

int *AllocateIntSEXP( size_t size, SEXP &S_exp, int &n_prot ) {
  PROTECT( S_exp = allocVector( INTSXP, size) ); n_prot++; 	  
  int *out = INTEGER_POINTER( S_exp );
  return out;
}

void AllocateAndStoreInSEXP( int *ptr, size_t size, SEXP &S_exp, int &n_prot ) {
  PROTECT( S_exp = allocVector( INTSXP, size) ); n_prot++; 	  
  int *out = INTEGER_POINTER( S_exp );
  for( int *it=ptr, *it_end = ptr + size; 
       it != it_end; 
       it++, out++) {
    *out =  *it;
  } 
}

void AllocateAndStoreInSEXP( double *ptr, size_t size, SEXP &S_exp, int &n_prot ) {
  PROTECT( S_exp = allocVector( REALSXP, size) ); n_prot++; 	  
  double *out = NUMERIC_POINTER( S_exp );
  for( double *it=ptr, *it_end = ptr + size; 
       it != it_end; 
       it++, out++) {
    *out =  *it;
  } 
}

int *AllocateAndStoreOccurInSEXP( int *ptr, size_t size, size_t remove_id, 
				  SEXP &S_exp, int &n_prot ) {

  PROTECT( S_exp = allocVector( INTSXP, size - 1) ); n_prot++; 	  
  int *out = INTEGER_POINTER( S_exp );
  int *ret = out;
  for( size_t i=0 ; i < size; i++, ptr++){
    if( i != remove_id )
      *out++ =  *ptr;
  }
  return ret;
}

int *AllocateAndStoreOccurInArray( int *ptr, size_t size, size_t remove_id, 
				  int *out) {

  int *ret = out;
  for( size_t i=0 ; i < size; i++, ptr++){
    if( i != remove_id )
      *out++ =  *ptr;
  }
  return ret;
}

// Compute gu function see E.birmele article
double computegU( double LambdaU, double DeltaU ) {
  double tU, hU, gU;

  tU = 1.0 + 1.0/ (FourPi * LambdaU) * 
    ( 1.0 + sqrt(1 + SixteenPi * LambdaU) );
  if( DeltaU <= tU ) {
    hU = 1;
  } else {
    hU = sqrt( DeltaU+1) / 
      ( sqrt( TwoPi*LambdaU ) * (DeltaU-1) );
  }
  gU = LambdaU * ( (1.0 + DeltaU) * log(1 + DeltaU) - DeltaU ) 
    - log(hU);
  return gU;

}

// Is "delclass" index contained in the Topological Class
bool isInDelClass( int delclass,  vector< int > &CT, int *perm) {

  for ( vector<int>::iterator it_ct = CT.begin(),
	  it_ct_end = CT.end();
	it_ct != it_ct_end; it_ct++) {
    if( perm[*it_ct] == delclass )
      return true;
  }
  return false;
}

bool isEqualCanonic( int* canonic, int *filter, int k) {

  for (int *it_c = canonic, *it_c_end = canonic + (k*k), *it_filter = filter; 
       it_c != it_c_end; 
       it_c++, it_filter++) {
    if( *it_c != *it_filter ) 
      return false;
  }

  return true;
 
}

//  Return a sorted occurrence list 
//  First the occurences are sorted according to the motif m (canonic form)
//  Then the occurences are sorted according to the motif m' or mp 
//  obtained by removing one node from m.
//  Then a last sort is done according to the location U (of mp)
//  mode
//    a - if "NodeToClass" and "NbrOfClasses" are provided the p-value
//        is computed
//    b - if  "StoreOcc" is motif occurence positions are stored
//    c - "motif_filter" if true only this motif is considered
//    d - "del_class" if true provided an deletion class is applied after the 
//        "motif filter" 
SEXP sort_m_mp_u( SparseAdjacency &G, int G_N, int k, 
		  vector<int*> *occurrences, 
		  double *Pi, int *NodeToClass, int NbrOfClasses, 
		  int *MotifFilter, int DelClassFilter, 
		  double PValue,
		  int &n_prot, bool StoreOcc) {


  int directed; 
  if(G.getSymmetry()) 
    directed=0 ; 
  else 
    directed=1;

  // Active Filters
  bool useMotifFilter=false, useDelClassFilter=false;
  if( MotifFilter ) {
    useMotifFilter = true;
    if( DelClassFilter >= 0 )
      useDelClassFilter=true;
  }
    

  // If LambdaU must be computed
  bool ComputePValue = false;
  int  *Alpha = new int[NbrOfClasses];

  // Init Alpha
  if ( (NodeToClass != 0) && (Pi != 0) && (NbrOfClasses != 0) ) {
    ComputePValue = true;
    for( int *it=Alpha, *it_end=&Alpha[NbrOfClasses]; it != it_end; it++) {
      *it=0;
    }
    for( int *it = NodeToClass, *it_end = &NodeToClass[G_N]; 
	 it != it_end; it++) {
      Alpha[*it]++;
    }
  } else {
    NbrOfClasses = 0; 
  }

  // Instanciate PseudoMixnet
  PseudoMixnet pmixnet(G_N, NbrOfClasses, Alpha, Pi, directed);

  SEXP  S_M_List;

  int64_t TotalCount        = occurrences -> size();
  vector<int*>::const_iterator it = occurrences -> begin();
  int *ptr = *it;

  // Print occurences
  // cerr << "TotalCount " << TotalCount << endl;

  // for ( int64_t count=0; count < TotalCount; count++, it++) {
  //   printList( "Occurrence : ", *it, *it + k, " ", true);
  // }

  int *mat;
  int *gamma = new int[k];

  // Store 2 references : on canonic form and on an occurence
  canonic_list_t *canonic_list = new canonic_list_t[TotalCount];

  // To Store canonic form
  int *canonic = new int[TotalCount * k*k];

  // To sort canonic form (used by canonic_order function)
  sort_canonic_nbr_elmnt = k*k;

  // tmp array to perform permutation
  int *tempo = new int[k];

  //
  // Compute canonic form and reoder occurrence indexes
  //
  for ( int64_t count=0; count < TotalCount; count++, it++) {

    // Extract motif 
    // The matrix must be transposed (true) for nauty
    mat = G.getMatrix( *it, k, false );

    // Get canonic form
    // ToDo : Perf call getCanonic once
    int *cur_canonic = &canonic[k*k * count];
    getCanonic( mat , k, directed, cur_canonic, gamma );

    delete [] mat;

    canonic_list[count].canonic    = cur_canonic;
    canonic_list[count].occurrence = *it;

    //
    // Print canonic form
    // printList( "Occurences, canonic ", *it, *it+k, " ", false); 
    // cout << ", " ;

    //
    // Permute nodes
    //
    int *node=*it;
    for( int i=0; i < k; i++){
      tempo[i] = node[gamma[i]];
    } 
    // Permute node in the occurence list
    for( int i=0; i < k; i++){
      node[i] = tempo[i];
    } 
    
    // Print cannonic
    // printList("", cur_canonic, cur_canonic + k*k, "", true);
    
    // Print gamma
    // printList( "  gamma : ", gamma, gamma + k);

  }

  //
  //  Sort occurrences according to the canonical form of m
  //

  qsort( canonic_list, TotalCount,  sizeof(canonic_list_t), canonic_order);

  //
  //  Print occurrences after the canonic sort and reordering
  //
# if( MSG > 0)
  cout << " Occurrences after the canonic sort and reordering" << endl;
  cout << " -------------------------------------------------" << endl;
  it = occurrences -> begin( );
  for ( int64_t count=0; count < TotalCount; count++, it++) {

    // Print node list
    ptr = canonic_list[count].occurrence;
    printList( "  ",  ptr, ptr + k, false);
    cout << ", ";

    // Print cannonic
    ptr = canonic_list[count].canonic;
    printList("", ptr, ptr + k*k, "", true);
  }
  cout << endl; 
# endif
  //
  //  Get refences where start and end a motif( canonic form)
  //

  //
  //  Build list of reference (pointer of type : canonic_list_t*)
  //  pointing to the same canonical form of m.
  //  Stored in the MotifEnd, MotifStart (STL) list
  //

  vector<canonic_list_t*> MotifEnd, MotifStart;
  int nbr_motif=0;
  it = occurrences -> begin( );
  canonic_list_t *prev_canonic = &canonic_list[0];
  MotifStart.push_back( &canonic_list[ 0 ] );
  for ( int64_t count=0; count < TotalCount; count++, it++) {

    //  If the current canonc form is equal to the previous one
    if( canonic_order( &canonic_list[count], prev_canonic ) ) {

      // New canonic form
      prev_canonic = &canonic_list[count];
      
      MotifEnd.push_back( &canonic_list[count] );
      MotifStart.push_back( prev_canonic  );
      nbr_motif = 0;
    } else{
      nbr_motif++;
    }
  }
  MotifEnd.push_back( &canonic_list[ TotalCount ] );


  //
  //  For each motif m, the occurence list of m is sorted 
  //  by ignoring one node (belonging to the same TC - the 
  //  remove index is stored in "remove" field) 
  // 

  vector<canonic_list_t*>::const_iterator ite     =  MotifEnd.begin();
  // unused vector<canonic_list_t*>::const_iterator ite_end =  MotifEnd.end();
  FindMotif TopoClassesBuider(TCMode, 0);
  int *ij_index = new int[ 2*k*k+2];
  int *Color  = new int[k];
  int *perm   = new int[k];
  int *perm_1 = new int[k];

  // Init
  for( int i=0; i < k; i++) {
      Color[i] = 0;
  }

  // Configure the occurence sort
  sort_occurences_nbr_elmnt = k;
 

  SEXP S_M_CanonicList; 

  if ( useMotifFilter ) {
    PROTECT( S_M_CanonicList = allocVector( VECSXP, 1 )); n_prot++;
    PROTECT( S_M_List        = allocVector( VECSXP, 1 )); n_prot++;
  } else {
    PROTECT( S_M_CanonicList=
	     allocVector( VECSXP,  MotifStart.size() )); n_prot++;
    PROTECT( S_M_List=
	     allocVector(VECSXP,  MotifStart.size() )); n_prot++;
  }

  // Counter of the sorted motif m to manage the SEXP 
  int m_index = 0;
  // Counter of the deletion class  mp to mana the corresponding SEXP
  int mp_index = 0;

  for( vector<canonic_list_t*>::const_iterator its =  MotifStart.begin(), 
	 its_end=  MotifStart.end();
       its != its_end; its++, ite++ ) {

    //
    // Occurences having the same motif
    //

    //
    //   Motif filter
    //
    if( ! useMotifFilter 
	|| isEqualCanonic( (*its) -> canonic, MotifFilter, k) ) {

      // Print cannonic
#     if( MSG > 0)
      ptr = (*its) -> canonic; 
      printList("Motif m=", ptr, ptr + k*k, "", false);
      cout << ", " << endl;
#     endif      

      //
      // Get a SparseAdjacency on the current motif m
      //
      ptr = (*its) -> canonic; 
      int pos = 0;
      for( int i=0 ; i < k*k; i++, ptr++){
	if(*ptr) {
	  // row index
	  ij_index[pos+1]   = i % k;
	  // col index
	  ij_index[pos] = i / k;
	  pos += 2;
	}
      } 
      // End list
      ij_index[pos]   = NOT_INDEX;
      ij_index[pos+1] = NOT_INDEX;
      
      // Buid the SparseAdjacency 
      SparseAdjacency motif( ij_index, k, (directed==0) );
      // Get good property to compute TC, 
      // "perm" give the permutation done
      motif.optimizeConnexity( perm ); 
      
      // get Topological Classes (TC)
      vector< int > *CT = TopoClassesBuider.getEquivalentNodes( motif, Color );
    
      size_t nbr_motif_m = *ite - *its;
      
      // 
      // motif.print(" motif" );
      // cout << "Occurences number of motif m " << nbr_motif_m  << endl; 
      
      // Compute inverse permutation perm_1
      for( int i=0; i < k; i++) {
	perm_1[perm[i]] = i; 
      }
      
      // Count the number of m' (mp)
      int n_mp=0;
      for( int i=0; i < k; i++) {
	if( CT[i].size() != 0) {
	  n_mp++;
	}
      }
    
      // unused SEXP S_Fuse_MP_List;
      SEXP S_MP_List;
      if ( useDelClassFilter ) {
	PROTECT( S_MP_List=allocVector(VECSXP, 1  )); n_prot++;
      } else {
	PROTECT( S_MP_List=allocVector(VECSXP, n_mp  )); n_prot++;
      }

      //
      //  For each Topological Class of m
      //
      
      mp_index = 0;
      for( int i=0; i < k; i++) {
	
	// Deletion class Filter
	if( ! useDelClassFilter ||  
	    isInDelClass( DelClassFilter, CT[i], perm_1 ) ) { 
	  // Computed for class of deletion
	  double max_t = -1;
	  
	  // Print  Topological Class
#         if( MSG > 0)
	  cout << "  Couple (m, m'),  DelClass " << i << " : ";
	  for ( vector<int>::iterator it_ct = CT[i].begin(),
	        it_ct_end = CT[i].end();
	      it_ct != it_ct_end; it_ct++) {
	   cout << *it_ct << ", ";
	  }
	  cout << endl;
#         endif	  
	
	  // Build [m,_mp] occurences list
	  m_mp_list_t *m_mp_list   = new m_mp_list_t[ CT[i].size() * nbr_motif_m ];
	  m_mp_list_t *it_m_mp_end = &m_mp_list[ CT[i].size() * nbr_motif_m ];
	  canonic_list_t *it_m = *its;
	  for( m_mp_list_t *it_m_mp = m_mp_list; 
	       it_m_mp != it_m_mp_end; it_m++) {
	    
	    // Scan all possible node in the same TC
	    // and specify the node index to remove in the occ. 
	    // list
	    for ( vector<int>::iterator it_ct = CT[i].begin(),
		    it_ct_end = CT[i].end();
		  it_ct != it_ct_end; it_ct++) {
	      it_m_mp -> occurrence = (it_m) -> occurrence;
	      it_m_mp -> remove     = perm_1[*it_ct];
	      it_m_mp++; 
	    }
	  }
	  
	  // sort occurences with one column removing (same CT)
	  // Sort according to U
	  qsort( m_mp_list , CT[i].size() * nbr_motif_m, 
		 sizeof(m_mp_list_t), occurrences_order);
	  
	  // Print after the sort
	  m_mp_list_t *prev = m_mp_list;
	  int remove;
	
	  // Store temporary the extension
	  // stack pb 
	  // int ExtStack[ it_m_mp_end - m_mp_list];
	  // SEXP S_Exts[it_m_mp_end - m_mp_list];
	  // SEXP S_U[it_m_mp_end - m_mp_list];

	  int *ExtStack = new  int[ it_m_mp_end - m_mp_list ];
	  SEXP *S_Exts  = new SEXP[ it_m_mp_end - m_mp_list ];
	  SEXP *S_U     = new SEXP[ it_m_mp_end - m_mp_list ];

	  // Fuse the 4 lists
	  SEXP S_U_Fuse_List=0;      
	  size_t ExtStackSize = 0;
	  size_t U_index = 0;
	  SEXP S_Del_List;

	  // Store LambdaU and DeltaU
	  double LambdaU[2];
	  LambdaU[0] = 0.0;  LambdaU[1] = 0.0;
	  double DeltaU, gU;
	  
	  int *occ_tmp = new int[ k-1 ];

	  // New extension
#         if( MSG > 0)
	  cout << "    Position U" << endl;
#         endif
	  // Test empty list
	  if( m_mp_list != it_m_mp_end ) {
	    
	    for( m_mp_list_t *it_m_mp = m_mp_list; 
		 it_m_mp != it_m_mp_end; it_m_mp++, it_m++ ) {
	      ptr    = it_m_mp -> occurrence;
	      remove = it_m_mp -> remove;
	      
	      // New position
	      if( occurrences_order( it_m_mp, prev ) ) {
#               if( MSG > 0)
		cout << "    Position U" << endl;
#               endif 
		
		// To store occurences
		int *it_occ;
		
		//
		// Store previous U (position)
		//
		
		// Allocate and store previous extension 
		if( StoreOcc ) {
		  
		  // Store deletion class
		  /* Inv
		  int* ct_ptr = AllocateIntSEXP( CT[i].size(), 
						 S_Del_List, n_prot);
		  
		  for ( vector<int>::iterator it_ct = CT[i].begin(),
			  it_ct_end = CT[i].end();
			it_ct != it_ct_end; it_ct++, ct_ptr++ ) {
		    *ct_ptr = perm_1[ *it_ct ];
		  }
		  */
		  
		  AllocateAndStoreInSEXP( ExtStack, ExtStackSize, 
					  S_Exts[U_index], 
					  n_prot );
		  
		  // Allocate and store previous 'remove index' 
		  // AllocateAndStoreInSEXP( RemStack, ExtStackSize, 
		  // S_Remove[U_index], 
		  //			    n_prot );
		  
		  it_occ = AllocateAndStoreOccurInSEXP( prev->occurrence, k, 
							prev -> remove, 
							S_U[U_index], 
							n_prot );
		  
		} else if (ComputePValue)  {

		  // Compute LambdaU
		  
		  // Allocate and store previous occurence U
		  // Need to be stored even if StoreOcc == false
		  it_occ = AllocateAndStoreOccurInArray( prev->occurrence, k, 
							 prev -> remove, 
							 occ_tmp );
	       

		  
		  // The removed node must be permuted according to the 
		  // OptimizeGeometry permutation
		  LambdaU[ 0 ] = pmixnet.computeLambdaU( 
							G, 
							it_occ, 
							k-1,
							// First Ext
							ExtStack[0], 
							NodeToClass );
		  
		  DeltaU = (ExtStackSize -  LambdaU[ 0 ]) /  LambdaU[ 0 ];
		  
		  LambdaU[ 1 ] = DeltaU;
		  gU = computegU( LambdaU[ 0 ], DeltaU );
		  max_t = MAX( max_t, gU );

#                 if( MSG > 0)
		  cout << "    Lambda U, Delta U, gU : " 
		       << LambdaU[ 0 ] << " "
		       << DeltaU << " " 
		       << gU << " " 
		       << endl;
#                 endif
		}
		U_index++;
		ExtStackSize = 0;
		
		prev =  it_m_mp;
	      }

	      // Print occurrence
#             if( MSG > 0)
	      ptr    = it_m_mp -> occurrence;
	      cout << "      ";
	      for( int i_i=0 ; i_i < k; i_i++, ptr++){
		if( i_i != remove )
		  cout << *ptr << ", ";
	      } 
	      ptr    = it_m_mp -> occurrence;
	      cout << "[" << ptr[remove] << "]";
	      cout << endl;
#             endif
	      // Store extension and remove index
	      ExtStack[ExtStackSize ] = ptr[remove];
	      // RemStack[ExtStackSize ] = remove;
	      ExtStackSize++;
	      
	    }
	    
	    //
	    // Store last U (position)
	    //
	    
	    // To store Occurences
	    int *it_occ;
	    bool store_list=true;
	    
	    if ( StoreOcc ) {
	      
	      // Store deletion class
	      int* ct_ptr = AllocateIntSEXP( CT[i].size(), 
					   S_Del_List, n_prot);
	      
	      for ( vector<int>::iterator it_ct = CT[i].begin(),
		      it_ct_end = CT[i].end();
		    it_ct != it_ct_end; it_ct++, ct_ptr++ ) {
		*ct_ptr = perm_1[ *it_ct ];
	      }
	  
	      // Allocate and store previous extension
	      AllocateAndStoreInSEXP( ExtStack, ExtStackSize, S_Exts[U_index], 
				      n_prot );
	      
	      // AllocateAndStoreInSEXP( RemStack, ExtStackSize, 
	      // S_Remove[U_index], n_prot );
	      
	      it_occ = AllocateAndStoreOccurInSEXP( prev->occurrence, k, 
						    prev -> remove, 
						    S_U[U_index], 
						    n_prot );
	    } else if (ComputePValue)  {
	      
	      // Allocate and store previous occurence U
	      // Need to be stored even if StoreOcc == false
	      it_occ = AllocateAndStoreOccurInArray( prev->occurrence, k, 
						     prev -> remove, 
						     occ_tmp );


	      // Compute LambdaU
	      
	      LambdaU[ 0 ] = pmixnet.computeLambdaU( 
						    G, 
						    it_occ, 
						    k-1,
						    // First Ext
						    ExtStack[0], 
						    NodeToClass );
	      DeltaU = (ExtStackSize -  LambdaU[ 0 ]) /  LambdaU[ 0 ];
		  
	      LambdaU[ 1 ] = DeltaU;
	      gU = computegU( LambdaU[ 0 ], DeltaU );
	      max_t = MAX( max_t, gU );
#             if( MSG > 0)
	      cout << "    Lambda U, Delta U, gU : " 
		   << LambdaU[ 0 ] << ", "
		   << DeltaU << ", " 
		   << gU << ", " 
		   << endl;
#             endif
	    }
	    
	    // Number of occurences (m, m', U)
	    size_t  NbrOfOccur = U_index + 1;
	
	    // Store arrays in (NbrOfOccur) lists  
	    // unused SEXP S_RemList;
	    SEXP S_UList, S_ExtList;
	  
	    if (StoreOcc) {
	      PROTECT( S_UList  =allocVector(VECSXP, NbrOfOccur  )); n_prot++;
	      PROTECT( S_ExtList=allocVector(VECSXP, NbrOfOccur  )); n_prot++;
	      // PROTECT( S_RemList=allocVector(VECSXP, NbrOfOccur)); n_prot++;

	      for( size_t i_u = 0; i_u < NbrOfOccur; i_u++) { 
		SET_VECTOR_ELT( S_ExtList, i_u, S_Exts  [ i_u ] );
		SET_VECTOR_ELT( S_UList,   i_u, S_U     [ i_u ] );
		// SET_VECTOR_ELT( S_RemList, i_u, S_Remove[ i_u ] );
	      } 

	      // Store in the main fields
	      PROTECT( S_U_Fuse_List = allocVector(VECSXP, 3) ); n_prot++;
	      SET_VECTOR_ELT( S_U_Fuse_List, 0, S_Del_List);
	      SET_VECTOR_ELT( S_U_Fuse_List, 1, S_UList   );
	      SET_VECTOR_ELT( S_U_Fuse_List, 2, S_ExtList );

	    } else if ( ComputePValue) {

	      // Compute E( mp ) - mp submotif
	      // and p-value
	      
	      //
	      //   Get a SparseAdjacency on the current mp motif 
	      //
	      
	      if ( CT[i].begin() != CT[i].end() ) {
		
		// No empty list
		int *ij_mp_index = new int[ 2*(k-1)*(k-1)+2];
		ptr = (*its) -> canonic; 
		int rem = perm_1[ *( CT[i].begin() )];
		int pos = 0;
		for( int u=0 ; u < k*k; u++, ptr++){
		  if(*ptr) {
		    // row index - discard the removed node
		    int ii = u % k;
		    // col index - discard the removed node
		    int jj =  u / k;
		    if ( (ii != rem) && (jj != rem) ) {
		      ij_mp_index[pos+1]   = ii - ( ii > rem );
		      ij_mp_index[pos] = jj - ( jj > rem );
		      pos += 2;
		    }
		  }
		}
		// End list
		ij_mp_index[pos]   = NOT_INDEX;
		ij_mp_index[pos+1] = NOT_INDEX;
		
		SparseAdjacency mp( ij_mp_index, k-1, (directed==0) );
		// mp.print("mp :");
		
		// Not optimal : FindNotif can't be used because works 
		// on a connex network 
		NonRedondantPermutation NRPMotif( mp );
		
		SEXP S_PValue;
		SEXP S_Filter;
		
		double esp_mp;
		double proba = 0.0;
		esp_mp = pmixnet.computeMean( mp, NRPMotif );
		proba = exp( - max_t) * esp_mp;
		
#               if( MSG > 0)
		cout << "  max_t, esp_mp, proba : " << max_t << " "
		     << esp_mp << " " 
		     << proba << endl;
#               endif
		if ( proba <= PValue) {
		  // Store PValue
		  double *p_value = AllocateDoubleSEXP(1, S_PValue, n_prot);
		  *p_value = proba;

		  // allocate filter
		  int *p_filter = AllocateIntSEXP(4, S_Filter, n_prot);
		  for( int *p_=p_filter, *p_end=p_filter+4; p_ != p_end; p_++) {
		    *p_ = 1;
		  } 
		
		  // Store deletion class
		  int* ct_ptr = AllocateIntSEXP( CT[i].size(), 
						 S_Del_List, n_prot);
		  for ( vector<int>::iterator it_ct = CT[i].begin(),
			  it_ct_end = CT[i].end();
			it_ct != it_ct_end; it_ct++, ct_ptr++ ) {
		    *ct_ptr = perm_1[ *it_ct ];
		  }
	  
		  
		  PROTECT( S_U_Fuse_List = allocVector(VECSXP, 3) ); n_prot++;
		  SET_VECTOR_ELT( S_U_Fuse_List, 0, S_Del_List);
		  SET_VECTOR_ELT( S_U_Fuse_List, 1, S_Filter);
		  SET_VECTOR_ELT( S_U_Fuse_List, 2, S_PValue);
		} else {
		  // empty list
		  store_list = false;
		}
		delete [] ij_mp_index;
	      } else {
		// empty list
		store_list = false;
	      }
	      
	    }
	    
	    // Add an item to MP_List
	    if(store_list) { 
	      SET_VECTOR_ELT( S_MP_List, mp_index, S_U_Fuse_List);
	      mp_index++;
	    }
	    
	    
	  } // End new extension :  if( m_mp_list != it_m_mp_end )

	  delete [] m_mp_list;
	  delete [] ExtStack; 
	  delete [] S_Exts; 
	  delete [] S_U;
	  delete [] occ_tmp;

	} // DelClassFilter

      } // For each topologic class (m')
      
      if ( mp_index != 0 ) {
	if( useMotifFilter )
	  SET_VECTOR_ELT( S_M_List, 0, S_MP_List);
	else
	  SET_VECTOR_ELT( S_M_List, m_index, S_MP_List);

	// Store Canonic or Adjacency matrix
	
	SEXP S_Canonic;
	
	// Transpose canonic
	int *transp = new int[k*k];
	int *can =  (*its) -> canonic ;
	for ( int i=0; i < k; i++){
	  for ( int j=0; j < k; j++){
	    //	transp[ i + k * j] = can[  j + k * i];
	    transp[ i + k * j] = can[  i + k * j];
	  }
	}
	AllocateAndStoreInSEXP( transp, k*k, S_Canonic, n_prot );
	delete [] transp;

	// Add to canonic list
	if (useMotifFilter) {
	  SET_VECTOR_ELT( S_M_CanonicList, 0, S_Canonic );
	} else {
	  SET_VECTOR_ELT( S_M_CanonicList, m_index, S_Canonic );
	}
      }
      delete [] CT;

    } // Motif Filter
    if ( mp_index != 0) {
      m_index++;
    }
  } // Motif loop

  delete [] canonic;
  delete [] canonic_list;
  SEXP S_Result = R_NilValue;
  // Fuse the two fields
  if ( m_index !=0 ) {
    PROTECT( S_Result = allocVector(VECSXP, 2  )); n_prot++;
    SET_VECTOR_ELT( S_Result, 0, S_M_CanonicList );
    SET_VECTOR_ELT( S_Result, 1, S_M_List );
  } 
    
  // To avoid warnings
  LOGMSG(100, std::cout, "NULL message ", "" );

  delete [] Alpha;
  delete [] gamma;
  delete [] tempo;
  delete [] ij_index;
  delete [] Color;
  delete [] perm;
  delete [] perm_1;
  
  return S_Result;
}

//  Return the counts of motifs 

SEXP count_m( SparseAdjacency &G, int k, vector<int*> *occurrences, int &n_prot) {

  // unused SEXP  S_M_List;

  int64_t TotalCount        = occurrences -> size();
  vector<int*>::const_iterator it = occurrences -> begin();
  // unused int *ptr = *it;


  // for ( int64_t count=0; count < TotalCount; count++, it++) {
  //  printList( "Occurrence : ", *it, *it + k, " ", true);
  // }

  int *mat;
  int *gamma = new int[k];

  // Store 2 references : on canonic form and on an occurence
  canonic_list_t *canonic_list = new canonic_list_t[TotalCount];

  // To Store canonic form
  int *canonic = new int[TotalCount * k*k];

  // To sort canonic form (used by canonic_order function)
  sort_canonic_nbr_elmnt = k*k;

  int directed; 
  if(G.getSymmetry()) 
    directed=0 ; 
  else 
    directed=1;

  // tempory storage for permutations, transpositions
  int *tempo = new int[k];
  int *transp = new int[k*k];



  //
  // Compute canonic form and reoder occurrence indexes
  //

  for ( int64_t count=0; count < TotalCount; count++, it++) {

    // Extract motif
    // The matrix must be transposed (true) for nauty
    mat = G.getMatrix( *it, k, false );

    // Get canonic form
    // ??? Perf call getCanonic once
    int *cur_canonic = &canonic[k*k * count];
    getCanonic( mat , k, directed, cur_canonic, gamma );

    delete [] mat;

    canonic_list[count].canonic    = cur_canonic;
    canonic_list[count].occurrence = *it;

    //
    // Print canonic form
    // printList( "Occurences, canonic ", *it, *it+k, " ", false); 
    // cout << ", " ;

    //
    // Permute nodes
    //
    int *node=*it;
    for( int i=0; i < k; i++){
      tempo[i] = node[gamma[i]];
    } 
    // Permute node in the occurence list
    for( int i=0; i < k; i++){
      node[i] = tempo[i];
    } 
    
    // Print cannonic
    // printList("", cur_canonic, cur_canonic + k*k, "", true);
    
    // Print gamma
    // printList( "  gamma : ", gamma, gamma + k);

  }

  //
  //  Sort occurrences according to the canonical form of m
  //

  qsort( canonic_list, TotalCount,  sizeof(canonic_list_t), canonic_order);

  //
  //  Get refences where start and end a motif( canonic form)
  //

  //
  //  Build list of reference (pointer of type : canonic_list_t*)
  //  pointing to the same canonical form of m.
  //  Stored in the MotifEnd, MotifStart (STL) list
  //

  vector<canonic_list_t*> MotifEnd, MotifStart;
  int nbr_motif=0;
  it = occurrences -> begin( );
  canonic_list_t *prev_canonic = &canonic_list[0];
  MotifStart.push_back( &canonic_list[ 0 ] );
  for ( int64_t count=0; count < TotalCount; count++, it++) {

    //  If the current canonc form is equal to the previous one
    if( canonic_order( &canonic_list[count], prev_canonic ) ) {

      // New canonic form
      prev_canonic = &canonic_list[count];
      
      MotifEnd.push_back( &canonic_list[count] );
      MotifStart.push_back( prev_canonic  );
      nbr_motif = 0;
    } else{
      nbr_motif++;
    }
  }
  MotifEnd.push_back( &canonic_list[ TotalCount ] );


  //
  //  For each motif m, the occurence list of m is sorted 
  //  by ignoring one node (belonging to the same TC - the 
  //  remove index is stored in "remove" field) 
  // 

  vector<canonic_list_t*>::const_iterator ite     =  MotifEnd.begin();
  // unused vector<canonic_list_t*>::const_iterator ite_end =  MotifEnd.end();


  SEXP S_M_CanonicList; 
  SEXP S_M_CountList;

  PROTECT( S_M_CanonicList=allocVector( VECSXP,  MotifStart.size() )); n_prot++;
  PROTECT( S_M_CountList=allocVector( INTSXP,  MotifStart.size() )); n_prot++;
  int *M_CountList = INTEGER_POINTER( S_M_CountList );

  // ID of the sorted motif m
  int m_index = 0;

  for( vector<canonic_list_t*>::const_iterator its =  MotifStart.begin(), 
	 its_end=  MotifStart.end();
       its != its_end; its++, ite++, m_index++ ) {

#   if( MSG)
    // Print cannonic
    ptr = (*its) -> canonic; 
    printList("", ptr, ptr + k*k, "", false);
    cout << ", " << endl ;
#   endif

    // Store Canonic or Adjacency matrix
    SEXP S_Canonic;
    // Transpose canonic
    int *can =  (*its) -> canonic ;
    for ( int i=0; i < k; i++){
      for ( int j=0; j < k; j++){
	// transp[ i + k * j] = can[  j + k * i];
	transp[ i + k * j] = can[  i + k * j];
      }
    }

    AllocateAndStoreInSEXP( transp, k*k, S_Canonic, n_prot );

    // Add to canonic list
    SET_VECTOR_ELT( S_M_CanonicList, m_index, S_Canonic );

    // Store counts
    M_CountList[ m_index ] =  (int) ((*ite) -(*its)) ;

  }

  delete [] canonic;
  delete [] canonic_list;

  SEXP S_Result;

  // Fuse the two fields
  PROTECT( S_Result = allocVector( VECSXP, 2 )); n_prot++;
  SET_VECTOR_ELT( S_Result, 0,  S_M_CanonicList );
  SET_VECTOR_ELT( S_Result, 1,  S_M_CountList );

  delete [] gamma;
  delete [] tempo;
  delete [] transp;

  return S_Result;
}

// Fill the DelClass vector with the first node of the DelClass 
void getDelClassNode( SEXP  S_m_mp, int m_mp_len, int *del_class) {
  
  for( int ind = 0; ind < m_mp_len; ind++) {
    SEXP S_m_mp_i = VECTOR_ELT( S_m_mp, ind );
    if( S_m_mp_i != R_NilValue ) {
      del_class[ ind ] = INTEGER_POINTER( 
					 VECTOR_ELT( 
						    S_m_mp_i 
						    , 0 )
					  ) [0] ;
    } else {
      del_class[ ind ] = -1;
    }
  }
}

void setMotifFilter( SEXP  S_m_mp, int index, int selected, int parent, int motif_id, int del_id ) {
  
    int *filter = INTEGER_POINTER( 
			       VECTOR_ELT( 
			         VECTOR_ELT( S_m_mp, index) 
				 , 1 )
			      );

    if ( filter[0] == 1 ) {
      filter[0] = selected;
      filter[1] = parent;
      filter[2] = motif_id;
      filter[3] = del_id;
    } 
}

// S_m_mp[ index ] is NULL, return DBL_MAX
double getPValue( SEXP  S_m_mp, int index ) {

  SEXP S_m_mp_i = VECTOR_ELT( S_m_mp, index) ; 
  double pvalue = DBL_MAX;
  pvalue = NUMERIC_POINTER( 
			   VECTOR_ELT( 
				      S_m_mp_i 
				      , 2 )
			    ) [0];
  
  return pvalue;
}

SEXP getExceptional( int *G_edges, int G_nnodes, 
		int M_Min_nnodes, int M_Max_nnodes, 
		int *NodeToClass, double *Pi, int NbrClasses, 
		double PValue, int Directed, 
		int &n_prot) {
  //
  // Get all occurrences of size k = M_nnodes
  //

  //
  //  Find all motif with size k
  //
  FindMotif find( StoreMode, 0);
  FindMotif f( CountMode, 0);

  vector<int*> *list;
  SparseAdjacency G( G_edges, G_nnodes, true );
  // G.print( "G");
  SparseAdjacency H( G_edges, G_nnodes, ! Directed  );
  // H.print( "H");

  SEXP S_ReturnList; 

  // All components are store even if they are NULL
  PROTECT( S_ReturnList = allocVector( VECSXP, 
				   M_Max_nnodes - M_Min_nnodes + 1)); n_prot++;

  // Index in the Returned list
  int i_k = 0;

  /// Previous non NULL value of k index
  int prev_k_index = -1;

  for ( int k= M_Min_nnodes; k <= M_Max_nnodes; k++, i_k++) {

    // To optimize the motifs
    // int perm_0[k-1];
    int *perm_0 = new int[k];
    int *perm_1 = new int[k];

    SEXP S_Result; 

    // Occurence list 'list' must be deallocated"
    // Doesn't belong to 'find'
    // TODO To improve
    find.clearFoundList();
    list = find.findAllMotifs( G, k ); 
    
    S_Result = sort_m_mp_u( H, G_nnodes, k, list, 
			    Pi, NodeToClass, NbrClasses, 
			    0, -1,
			    PValue,
			    n_prot, false );

    if ( S_Result != R_NilValue ) {
      // Store k value
      SEXP S_k;
      PROTECT( S_k = allocVector( INTSXP, 1) ); n_prot++; 	  
      int *p_k = INTEGER_POINTER( S_k );
      p_k[0] = k;

      SEXP S_Struct; 
      PROTECT( S_Struct = allocVector( VECSXP, 2 ));  n_prot++;;
      SET_VECTOR_ELT( S_Struct, 0, S_k ); 
      SET_VECTOR_ELT( S_Struct, 1, S_Result ); 
      SET_VECTOR_ELT( S_ReturnList, i_k, S_Struct );

      SEXP S_Adj_1 = VECTOR_ELT( S_Result, 0 ); 
      int MotifNbr_1 = LENGTH( S_Adj_1 );

      // Skip First motif
      // if ( prev_k_index != -1 ) {

      for ( int k_index = 0; k_index <= (prev_k_index); k_index++ ) {

	SEXP S_k_value = VECTOR_ELT( 
				    VECTOR_ELT( S_ReturnList, 
						k_index ), 
				    0 );
	int prev_k_value = INTEGER_POINTER( S_k_value ) [0];
 
	// Get Adjacency matrix with i-1 nodes
	SEXP S_Adj_0 = VECTOR_ELT( 
				  VECTOR_ELT( 
					     VECTOR_ELT( S_ReturnList, 
							 k_index ), 
					     1 ), 
				  0);
	// Get M_list 
	SEXP S_m_0 = VECTOR_ELT( 
				VECTOR_ELT( 
					   VECTOR_ELT( S_ReturnList,
						       k_index ), 
					   1 ), 
				1);

	int MotifNbr_0 = LENGTH( S_Adj_0 );


	for(int m = 0; m <  MotifNbr_0; m++ ) {

	  // Test NULL SEXP
	  if( VECTOR_ELT(S_Adj_0, m) !=  R_NilValue ) {
	    // int debug = LENGTH( VECTOR_ELT(S_Adj_0, m) );
	    int *adj_0 =  INTEGER_POINTER( VECTOR_ELT(S_Adj_0, m) );
	    // printList( "canonic ", adj_0, adj_0+ (k-1)*(k-1), " ", true); 

	    SparseAdjacency motif0( adj_0, prev_k_value, "dense", (Directed == 0) );
	    SparseAdjacency motif0_canon( motif0 );
	    motif0.optimizeConnexity( perm_0 );

	    // motif0.print( "motif opt 0" );
	    // printList( " perm 0 ", perm_0, perm_0 + (k-1), true );

	    // Number of Del class
	    SEXP S_m_mp_0 = VECTOR_ELT( S_m_0, m );
	    int m_mp_len_0 = LENGTH( S_m_mp_0 );
	    int *del_class_0 = new int[ m_mp_len_0 ];

	    // get the first node ID of the del class
	    getDelClassNode( S_m_mp_0, m_mp_len_0, del_class_0 );
	       

	      // printList( "Del class 0 ", del_class_0, 
	      //	   del_class_0 + m_mp_len_0, " ", true);
 
	      // Get Degree of delete node
	      // int deg_in_0  = 0;
	      // int deg_out_0 = 0;

	      for(int p = 0; p <  MotifNbr_1; p++ ) {

		// Test NULL SEXP
		if( VECTOR_ELT(S_Adj_1, p) !=  R_NilValue ) {

		  //		  int debug_1 = LENGTH( VECTOR_ELT(S_Adj_1, p) );
		  int *adj_1  = INTEGER_POINTER( VECTOR_ELT(S_Adj_1, p) );


#if            ( MSG > 0 )
		  printList( "canonic 1 ", adj_1, adj_1 + (k)*(k), " ", true); 
#endif

		  SparseAdjacency motif1( adj_1, k, "dense", (Directed == 0) );
		  SparseAdjacency motif1_canonic( motif1 );

		  motif1.optimizeConnexity( perm_1 );

		  // get M List
		  SEXP S_m_1 = VECTOR_ELT( 
					  VECTOR_ELT( 
						     VECTOR_ELT( S_ReturnList, 
								 i_k ), 
						     1 ), 
					  1);

		  // Number of Del class
		  SEXP S_m_mp_1 = VECTOR_ELT( S_m_1, p );
		  int m_mp_len_1 = LENGTH( S_m_mp_1 ) ; 
		  int *del_class_1 = new int [ m_mp_len_1 ];

		  // get the first node ID of the del class
		  getDelClassNode( S_m_mp_1, m_mp_len_1, del_class_1 ); 
	  
		    // printList( "Del class 1 ", del_class_1, 
		    //	     del_class_1 + m_mp_len_1, " ", true); 
	    
		    // printList( "canonic ", adj_0, adj_0 + (k-1)*(k-1), " ", true); 
		    // printList( "canonic ", adj_1, adj_1 + (k)*(k), " ", true); 
	    

		    for( int del_0 = 0;  del_0 <  m_mp_len_0; del_0++) {

		      if( del_class_0[del_0] >= 0 ) {
			// Get Degree of delete node
			int deg_in_0  = 0;
			int deg_out_0 = 0;
			if ( Directed ) {
			  deg_in_0  = motif0_canon.getRowSize( del_class_0[ del_0 ] );
			  deg_out_0 = motif0_canon.getColSize( del_class_0[ del_0] );
			} else {
			  deg_in_0  = motif0_canon.getAllRowSize( del_class_0[ del_0 ] );
			}
	      
			for( int del_1 = 0;  del_1 <  m_mp_len_1; del_1++) {
			  if( del_class_1[ del_1 ] >= 0 ) {		    
			    if( getPValue( S_m_mp_1, del_1 ) <= PValue ) {
		      
			      if( getPValue( S_m_mp_0, del_0 ) <= PValue ) {
		    
				// Get Degree of delete node
				int deg_in_1 = 0;
				int deg_out_1 = 0;
				if ( Directed ) {
				  deg_in_1  = motif1_canonic.getRowSize( del_class_1[ del_1 ] );
				  deg_out_1 = motif1_canonic.getColSize( del_class_1[ del_1 ] );
				} else {
				  deg_in_1  = motif1_canonic.getAllRowSize( del_class_1[ del_1 ] );
				}
			
		    
				// cout << "del_0 :" << del_0 << ", del_1 :" << del_1 << endl;
			
				if ( ( deg_in_0 == deg_in_1 ) && (deg_out_0 == deg_out_1 ) ) {
		      
				  int *color_0 = new int [ prev_k_value ];
				  int *color_1 = new int [ k ];
		      
				  for( int *ptr=color_0, 
					 *ptr_end= color_0 + (prev_k_value); 
				       ptr != ptr_end; ptr++ ) {
				    *ptr=0;
				  }
				  for( int *ptr=color_1, *ptr_end= color_1 + (k); 
				       ptr != ptr_end; ptr++ ) {
				    *ptr=0;
				  }
		      
				  color_0[ perm_0[ del_class_0[ del_0 ] ]] = 1;
				  color_1[ perm_1[ del_class_1[ del_1 ] ]] = 1;
				  // printList( " color_0 ", color_0,color_0+ (k-1), true); 
				  // printList( " color_1 ", color_1,color_1+ (k), true); 
		      
		    
				  bool incl = f.isExactlyIncluded( motif0, motif1, color_0, color_1 );
				  if( incl ) {
				    setMotifFilter( S_m_mp_1, del_1, 0, 1, m, del_class_0[ del_0 ] );
				  } else {
				    setMotifFilter( S_m_mp_1, del_1, 1, 0, 0, 0 );
				  }
		      
				  // cout << "( k=" <<  k-1 <<  ", motif=" << m << ", del=" 
				  //        << del_class_0[ del_0 ]  << ") ";
				  // cout << "( k=" <<  k <<  ", motif=" << p  << ", del=" 
				  //     << del_class_1[ del_1 ] << ") : " <<  incl << endl ;
				  delete [] color_0;
				  delete [] color_1;
				    
				} else {
				  // Different degrees => over-represented
				  setMotifFilter( S_m_mp_1, del_1, 1, 0, 0, 0);
				}
		    
			      } else {
				// Parent motif not over-represented
				setMotifFilter( S_m_mp_1, del_1, 1, 0, 0, 0);
			      }
		  
			    } else {
			      // motif not over-represented
			      setMotifFilter( S_m_mp_1, del_1, 0, 0, 0, 0);	
			    }
			  }  // Valid DelClass - if( del_class_1[del_1] >= 0 ) 
			}    // for( int del_1 = 0; ...
		      }  // Valid DelClass - if( del_class_0[del_0] >= 0 ) 
		    }    // for( int del_0 = 0; ...

		    delete [] del_class_1;

		}   //  if( VECTOR_ELT(S_Adj_1, p) !=  R_NilValue ) {
	      }     //  for(int p = 0; p <  MotifNbr_1; p++ )

	      delete [] del_class_0;
	  }    // if( VECTOR_ELT(S_Adj_0, m) !=  R_NilValue ) {
	}      // for(int m = 0; m <  MotifNbr_0; m++ ) {
     
      }
      if (  prev_k_index == -1 ) {

	// First value of k

	SEXP S_m_0 = VECTOR_ELT( 
				VECTOR_ELT( 
					   VECTOR_ELT( S_ReturnList, i_k ), 
					   1 ), 
				1);

	int MotifNbr_0 = LENGTH( S_m_0 );

	for(int m = 0; m <  MotifNbr_0; m++ ) {
	  SEXP  S_m_mp_0 = VECTOR_ELT( S_m_0, m );
	  if (  S_m_mp_0 != R_NilValue ) {
	    int m_mp_len_0 = LENGTH( S_m_mp_0 );
	    for( int del_0 = 0;  del_0 <  m_mp_len_0; del_0++) {

	      if ( VECTOR_ELT( S_m_mp_0, del_0) != R_NilValue ) {
		if( getPValue( S_m_mp_0, del_0 ) < PValue ) {
		  setMotifFilter( S_m_mp_0, del_0, 1, 0, 0, 0);
		} else {
		  setMotifFilter( S_m_mp_0, del_0, 0, 0, 0, 0);
		}
	      }
	    }
	  }
	}  // for(int m = 0; m <  MotifNbr_0; m++ )

      }
      prev_k_index = i_k;
    } else {
      // S_Result NULL
    }
    // for ( vector<int*>::iterator it = list->begin(),
    //	    it_end=list->end();it != it_end; it++ ) {
    //  delete [] *it;
    // }
    FindMotif::deallocateFoundList( list );
    delete [] perm_0;
    delete [] perm_1;
  }
  
  return S_ReturnList;
}
