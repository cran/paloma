# ifndef _SHRINK_MOTIF_H_
# define _SHRINK_MOTIF_H_

# include <stdint.h>
# include "sparse-adjacency.h"

// ??? static int *ReverseOccurrence;

template <typename T>
void Swap(T &a, T &b)
{
  T x = a;
  a = b;
  b = x;
}

struct Neibourgh
{
int Target;
int *Sources;
  int NbSources;

  Neibourgh(){
    Sources = NULL;
    NbSources = 0;
  };

  ~Neibourgh() {};

  void DestroyMe() { 
    if (Sources != NULL)
      delete[] Sources; 
  } ;

  Neibourgh operator =(const Neibourgh &Original)
  {
    Target = Original.Target;
    Sources = Original.Sources;
    NbSources = Original.NbSources;
    return *this;
  } ;

  bool operator()(const Neibourgh &Left, const Neibourgh &Right)
  {
    if (Left.Target < Right.Target)
      return true;
    return false;
  };
};

/*
void Neibourgh::DestroyMe()
{
  if (Sources != NULL)
    delete[] Sources;
}
*/

 /*
bool operator<(const Neibourgh &Left, const Neibourgh &Right)
{
if (Left.Target < Right.Target)
return true;
return false;
}
 */

 unsigned long int ComputeEnumeration(const SparseAdjacency &G, Neibourgh *MultipleVertices, int FirstNeibourgh, int LastNeibourgh, int *OccurrenceVertices, int OccurenceSize, const SparseAdjacency &M, int *NbMultipleVerticesAsNeibourghs, int *ReverseOccurrence, int *NbAffectedVertices = NULL, int *NbNeibourghsInOccurrence = NULL);

unsigned long int ComputeNbOccurrences(const SparseAdjacency &G, int *OccurrenceVertices, int OccurenceSize, const SparseAdjacency *M, int *NbNeibourghsInOccurrence, Neibourgh *MultipleVerticesTwoSources, Neibourgh *MultipleVerticesManySources, int *ReverseOccurrence);

# endif 
