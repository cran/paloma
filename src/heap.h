//  heap.h
//
//  Copyright (C) 2009 Michel Koskas
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

//! \file    sparse-adjacency.h
//! \author  Michel Koskas
//! \brief   Implement a sorted heap
//! \version 0.1
//! \date    December 2009

# include <list>
# include <stdint.h>
# include "sparse-adjacency.h"

# ifndef _S_HEAP_H_
# define _S_HEAP_H_

// TODO :
//   optimize inline, refactory, template, tests

//! \brief Make intersection of sets
//! The intersection operation use a sorted heap

class Heap {

public:
  const int **Lists;
  const int *ListsSizes;
  int *IndexInList;
  int NbLists;
  int *MyHeap;
  int HeapSize;
  bool ComputationFinished;
  Heap(const int **L,  const int *TL, int NbL);
  // Michel
  //  Heap(int **L,  int *TL, int NbL);
  ~Heap();
  void inline AddNode();
  int NextCommonElement();
  void Debug();
  void Union(int *Result, int &CardinalUnion);
  int *Union(int &CardinalUnion);

private:
  bool AllElementsEqual();
  void  RemoveHead();
};

bool IsIn(int *T, int TS, int x);

/*
// inline ???
void inline Swap(int &a, int &b)
{
  int x = a;
  a = b;
  b = x;
}

int Heap::NextCommonElement()
{
  if (ComputationFinished)
    return EMPTY_NODE;
  while ((!AllElementsEqual()) && (!ComputationFinished))
    RemoveHead();
  if (!ComputationFinished)
  {
    size_t Res = Lists[MyHeap[0]][IndexInList[MyHeap[0]]];
    for (size_t i = 0; i < NbLists; i++)
      RemoveHead();
    return Res;
  }
  return EMPTY_NODE;
}

void Heap::RemoveHead()
{
  if (ComputationFinished)
    return;
  IndexInList[MyHeap[0]]++;
  
  // ???
  // cout << " RemoveHeap : " << IndexInList[MyHeap[0]] << "==" << ListsSizes[MyHeap[0]] << endl;

  if (IndexInList[MyHeap[0]] == ListsSizes[MyHeap[0]] )
  {
    ComputationFinished = true;
    return;
  }
  size_t CurIndex = 0;
  if (2 * CurIndex + 1 < HeapSize)
  {
    while (2 * CurIndex + 2 < HeapSize)
    {
      size_t MinIndex = 2 * CurIndex + 1;
      if (Lists[MyHeap[MinIndex]][IndexInList[MyHeap[MinIndex]]] > Lists[MyHeap[MinIndex + 1]][IndexInList[MyHeap[MinIndex + 1]]])
        MinIndex++;
      if (Lists[MyHeap[CurIndex]][IndexInList[MyHeap[CurIndex]]] > Lists[MyHeap[MinIndex]][IndexInList[MyHeap[MinIndex]]])
      {
        Swap(MyHeap[CurIndex], MyHeap[MinIndex]);
        CurIndex = MinIndex;
      }
      else
        break;
    }
    if (2 * CurIndex + 1 < HeapSize)
      if (Lists[MyHeap[CurIndex]][IndexInList[MyHeap[CurIndex]]] > Lists[MyHeap[2 * CurIndex + 1]][IndexInList[MyHeap[2 * CurIndex + 1]]])
        Swap(MyHeap[CurIndex], MyHeap[2 * CurIndex + 1]);
  }
}
*/
#endif
