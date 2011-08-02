//  concentration.h
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

//! \file    concentration.h
//! \author  Gilles Grasseau 
//! \brief   Implement all routines for motif concentration
//! \version 0.1
//! \date    December 2010

SEXP sort_m_mp_u( SparseAdjacency &G, int G_N, int k, 
		  vector<int*> *occurrences, 
		  double *Pi, int *NodeToClass, int NbrOfClasses, 
		  int *MotifFilter, int DelClassFilter,
		  double PValue,
		  int &n_prot, bool StoreOcc);

SEXP count_m( SparseAdjacency &G, int k, 
	      vector<int*> *occurences, 
	      int &n_prot_ptr );

SEXP getExceptional( int *G_edges, int G_nnodes, 
		 int M_Min_nnodes, int M_Max_nnodes, 
		 int *NodeToClass, double *Pi, int NbrClasses, 
		 double PValue, int Directed, 
		 int &n_prot);
