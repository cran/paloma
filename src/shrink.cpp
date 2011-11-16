# include "shrink.h"
# include "heap.h"

uint64_t Factorielle(int n)
{
  uint64_t Res = 1;
  for (int i = 1; i <= n; i++)
    Res *= i;
  return Res;
}

uint64_t Binomiale(int n, int p)
{
  if ((p > n) || (p < 0))
    return 0;
  int Greater = p > n - p ? p : n - p;
  int Lower = n - Greater;
  uint64_t Res = 1;
  for (int i = 0; i < Lower + 1; i++)
    Res = (Res * (Greater + i)) / (i + 1);
  return Res;
}

int *IntersectSortedArrays(const int *Left, int NbL, const int *Right, int NbR, int &ResultSize)
{
  int MinSize = NbL;
  if (NbR < NbL)
    MinSize = NbR;
int *Result = new int[MinSize];
  int Index = 0;
  int LeftIndex = 0, RightIndex = 0;
while((LeftIndex < NbL) && (RightIndex < NbR))
{
if (Left[LeftIndex] == Right[RightIndex])
{
Result[Index++] = Left[LeftIndex];
LeftIndex++;
RightIndex++;
continue;
}
if (Left[LeftIndex] < Right[RightIndex])
LeftIndex++;
else
RightIndex++;
}
  ResultSize = Index;
return Result;
}


/*
unsigned long int ComputeEnumeration(Neibourgh *MultipleVertices, int FirstNeibourgh, int LastNeibourgh, int *OccurrenceVertices, int OccurenceSize, const SparseAdjacency &motif, int *NbMultipleVerticesAsNeibourghs, int *NbAffectedVertices, int *NbNeibourghsInOccurrence)
*/
unsigned long int ComputeEnumeration(const SparseAdjacency &G, 
				     Neibourgh *MultipleVertices, int FirstNeibourgh, int LastNeibourgh, int *OccurrenceVertices, int OccurenceSize, const SparseAdjacency &motif, int *NbMultipleVerticesAsNeibourghs, int *ReverseOccurrence, int *NbAffectedVertices, int *NbNeibourghsInOccurrence)
{
  // GG
  // M.NbPoorNeibourghs - > NbPoorNeibourghs
  // OccurrenceVertices -> MapToG
  const int64_t *M_NbPoorNeibourghs =  motif.getNodeValues();

  if (FirstNeibourgh == LastNeibourgh)
  {
    unsigned long int PartialResult = 1;
    //    for (int i = 0; i < M.NbVertices; i++)
    for (int i = 0; i < motif.getNbrRow(); i++)
    {      
      for (int j = 0; j < M_NbPoorNeibourghs[i] - NbAffectedVertices[i]; j++)
      {
        // PartialResult *= (NbAdj[OccurrenceVertices[i]] - NbMultipleVerticesAsNeibourghs[i] - NbNeibourghsInOccurrence[i] - j);
        PartialResult *= ( G.getAllRowSize( OccurrenceVertices[i] ) - NbMultipleVerticesAsNeibourghs[i] - NbNeibourghsInOccurrence[i] - j);
        PartialResult /= (j + 1);
      }
      PartialResult *= Factorielle( M_NbPoorNeibourghs[i] );
    }
    if (PartialResult < 0)
    {
      std::cerr << "There's a bug here, folks. It's just about time to take a look." << std::endl;
      std::cerr << "Take courage. I'm leaving with errcode 272" << std::endl;
      // GG
      // exit(272);
      return( -1 );
    }
    // PrintSituation(OccurrenceVertices, OccurenceSize, NbMultipleVerticesAsNeibourghs, NbAffectedVertices, NbNeibourghsInOccurrence, PartialResult);
    // char c;
    // cin >> c;
    return PartialResult;
  }
  unsigned long int PartialResult = 0;
  PartialResult += ComputeEnumeration(G,MultipleVertices, FirstNeibourgh + 1, LastNeibourgh, OccurrenceVertices, OccurenceSize, motif, NbMultipleVerticesAsNeibourghs, ReverseOccurrence, NbAffectedVertices, NbNeibourghsInOccurrence);

  
  for (int i = 0; i < MultipleVertices[FirstNeibourgh].NbSources; i++)
  {
    int I = MultipleVertices[FirstNeibourgh].Sources[i];
    NbAffectedVertices[ReverseOccurrence[I]]++;
    if (M_NbPoorNeibourghs[ReverseOccurrence[I]] >= NbAffectedVertices[ReverseOccurrence[I]])
      PartialResult += ComputeEnumeration(G, MultipleVertices, FirstNeibourgh + 1, LastNeibourgh, OccurrenceVertices, OccurenceSize, motif, NbMultipleVerticesAsNeibourghs, ReverseOccurrence, NbAffectedVertices, NbNeibourghsInOccurrence);
    NbAffectedVertices[ReverseOccurrence[I]]--;
  }
  return PartialResult;
}


// Dans cette fonction, NbNeibourghsInOccurrence est, pour chaque sommet de l'occurrence du motif réduit, le nombre de ses voisins appartenant au motif réduit.
// Modification 10 mai 2011 : on passe en argument deux familles de MultipleVertices : ceux à deux sources et ceux à nombre quelconque de sources. 
// Dansles deux cas, le cardinal en est majpré par n(n-1)/2 avec n le nombre de sommets du motif.
unsigned long int ComputeNbOccurrences(const SparseAdjacency &G, int *OccurrenceVertices, int OccurenceSize, const  SparseAdjacency *motif, int *NbNeibourghsInOccurrence, Neibourgh *MultipleVerticesTwoSources, Neibourgh *MultipleVerticesManySources, int *ReverseOccurrence)
{
  // GG
  const int64_t *M_NbPoorNeibourghs =  motif->getNodeValues();


  int NbMultipleVertices = 0;
  
  //TODO: ce qui suit n'est à faire quepour les sommets ayant au moins un PoorNeibourgh.
  Neibourgh *MultipleVertices = MultipleVerticesTwoSources;
  // On commence par calculer la liste des sommets du graphe ayant au moins deux voisins dans l'occurrence du motif réduit.

  
for (int i = 0; i < OccurenceSize; i++)
    if (M_NbPoorNeibourghs[i] > 0)
      for (int j = i + 1; j < OccurenceSize; j++)
        if (M_NbPoorNeibourghs[j] > 0)
        {
          // Source1 et Source2 sont deux sommets de l'occurrence du motif réduit
          int Source1 = OccurrenceVertices[i];
          int Source2 = OccurrenceVertices[j];

          
          // On calcule la liste de leurs voisins communs
          int IntersectionListCardinal;
	  // GG
	  /*
          int *IntersectionList = IntersectSortedArrays(MyAdj[Source1], NbAdj[Source1], MyAdj[Source2], NbAdj[Source2], IntersectionListCardinal);
	  */
          int *IntersectionList = IntersectSortedArrays(
		   G.getAllRow(Source1), G.getAllRowSize( Source1 ), 
		   G.getAllRow(Source2), G.getAllRowSize( Source2 ),
		   IntersectionListCardinal);
          // On en retire le nombre de voisins communs qui sont dans l'occurrence du motif réduit
          int NbToRemove = 0;
          for (int i = 0; i < IntersectionListCardinal; i++)
            if (ReverseOccurrence[IntersectionList[i]] != -1)
              NbToRemove++;
          // On en déduit le nombre de sommets n'appartenant pas à l'occurrence du motif réduit et ayant au moins deux voisins dans l'occurrence du motif réduit.
          NbMultipleVertices += IntersectionListCardinal - NbToRemove;
          delete[] IntersectionList;
        }
  if (NbMultipleVertices > 0)
  {
    // Maintenant on va stocker les sommets du graphe n'appartenant pas à l'occurrence du motif réduit et ayant au moins deux voisins dedans.
    // POur chaque sommet de l'occurrence du motif réduit on veut calculer le nombre de ses voisins n'appartenant pas à l'occurrence mais étant aussi voisins d'au moins un autre sommet de l'occurrence du motif réduit.
    int NumMultipleVertex = 0;
    // On alloue la place nécessaire
    // TODO: free this mem (next line)
    int *NbMultipleVerticesAsNeibourghs = new int[OccurenceSize];
    // On initialise.
    for (int i = 0; i < OccurenceSize; i++)
      NbMultipleVerticesAsNeibourghs[i] = 0;
    // On y va.
    for (int i = 0; i < OccurenceSize; i++)
      if (M_NbPoorNeibourghs[i] > 0)
        for (int j = i + 1; j < OccurenceSize; j++)
          if (M_NbPoorNeibourghs[j] > 0)
          {
            // Source1 et Source2 sont à nouveau deux sommets de l'occurrence du motif réduit.
            int Source1 = OccurrenceVertices[i];
            int Source2 = OccurrenceVertices[j];
            // On calcule la liste de leurs voisins communs.
            int IntersectionListCardinal; 
	    // GG
	    /*
            int *IntersectionList = IntersectSortedArrays(MyAdj[Source1], NbAdj[Source1], MyAdj[Source2], NbAdj[Source2], IntersectionListCardinal);
	    */

	    int *IntersectionList = IntersectSortedArrays(
		   G.getAllRow(Source1), G.getAllRowSize( Source1 ), 
		   G.getAllRow(Source2), G.getAllRowSize( Source2 ),
		   IntersectionListCardinal);

            // On parcourt cette liste
            if (IntersectionListCardinal > 0)
              for (int i = 0; i < IntersectionListCardinal; i++)
              {
                int I = IntersectionList[i];
                // Si l'élément courant n'est pas dans l'occurrence du motif réduit
                if (ReverseOccurrence[I] == -1)
                {
                  // On l'ajoute à la liste des MultipleVertices
                  MultipleVertices[NumMultipleVertex].Target = I;
                  if (Source1 > Source2)
                    Swap<int>(Source1, Source2);
                  MultipleVertices[NumMultipleVertex].Sources = new int[2];
                  MultipleVertices[NumMultipleVertex].Sources[0] = Source1;
                  MultipleVertices[NumMultipleVertex].Sources[1] = Source2;
                  NumMultipleVertex++;
                }
              }
            delete[] IntersectionList;
          }
    // Il faut maintenant régler le cas des sommets extérieurs à l'occurrence et ayant au moins trois voisins dans l'occurrence.
    // Ces cas se reconnaissent au fait qu'on a des Neibourghs ayant la même valeur de Target. On trie donc MultipleVertices en comparant les Target.
    sort(MultipleVertices, MultipleVertices + NbMultipleVertices, MultipleVertices[0]);
    int NbNewMultipleVertices = NbMultipleVertices;
    // On calcule le nombre de fois qu'un Target apparait plusieurs fois dans les Neibourghs et on en déduit le nombre de Neibourghs qu'on obtiendra après simplification. 
    for (int i = 1; i < NbMultipleVertices; i++)
      if (MultipleVertices[i - 1].Target == MultipleVertices[i].Target)
        NbNewMultipleVertices--;
    // S'il y a simplification possible on la fait.
    if (NbNewMultipleVertices < NbMultipleVertices)
    {
      int NbMultipleVerticesInVector = NbMultipleVertices;
      int Begin, End = 0;
      int CurrentNewMultipleVertexIndex = 0;
      while (End < NbMultipleVerticesInVector)
      {
        Begin = End;
        End++;
        while ((End < NbMultipleVerticesInVector) && (MultipleVertices[Begin].Target == MultipleVertices[End].Target))
          End++;
        int **TheSources = new int*[End - Begin], *TheSizes = new int[End-Begin];
        for (int i = 0; i < End - Begin; i++)
        {
          TheSizes[i] = 2;
          TheSources[i] = MultipleVertices[Begin + i].Sources;
        }
        MultipleVerticesManySources[CurrentNewMultipleVertexIndex].Target = MultipleVertices[Begin].Target;
        Heap H(
	       (const int** ) TheSources, 
	       (const int*) TheSizes, 
	       End - Begin);
        H.Union(MultipleVerticesManySources[CurrentNewMultipleVertexIndex].Sources, MultipleVerticesManySources[CurrentNewMultipleVertexIndex].NbSources);
        delete[] TheSources;
        delete[] TheSizes;
        CurrentNewMultipleVertexIndex++;
      }
      MultipleVertices = MultipleVerticesManySources;
      NbMultipleVertices = NbNewMultipleVertices;
    }

    
    // On calcule maintenant pour chaque sommet de l'occurrence du motif simplifie le nombre de ses voisins qu'il partage avec un autre sommet de l'occurrence du motif simplifie.
    for (int NumMultipleVertex = 0; NumMultipleVertex < NbMultipleVertices; NumMultipleVertex++)
      for (int i = 0; i < MultipleVertices[NumMultipleVertex].NbSources; i++)
        NbMultipleVerticesAsNeibourghs[ReverseOccurrence[MultipleVertices[NumMultipleVertex].Sources[i]]]++;
    uint64_t Result;

    
    // On calcule enfin le nombre d'occurrences du motif non simplifie.
    int *NbAffectedVertices = new int[OccurenceSize];
    for (int i = 0; i < OccurenceSize; i++)
      NbAffectedVertices[i] = 0;

    Result = ComputeEnumeration( 
		     G, MultipleVertices, 0, NbMultipleVertices, 
		     OccurrenceVertices, OccurenceSize, *motif, 
		     NbMultipleVerticesAsNeibourghs, 
		     ReverseOccurrence,
		     NbAffectedVertices, 
		     NbNeibourghsInOccurrence);
    delete[] NbAffectedVertices;
    delete[] NbMultipleVerticesAsNeibourghs;
    return Result;
  }
else
  {
    uint64_t PartialResult = 1;
    for (int i = 0; i < motif->getNbrRow(); i++)
      for (int j = 0; j < M_NbPoorNeibourghs[i]; j++)
      {// Modified on May 4th 2011 (M->MyDegre[i] transformed in M->NbAdj[i])
	/* GG
        PartialResult *= (NbAdj[OccurrenceVertices[i]] - M->MyDegre[i] - j);
	*/
        PartialResult *= (  G.getAllRowSize( OccurrenceVertices[i] )  - motif->getAllRowSize(i) - j);
        // PartialResult /= (j + 1);
      }
    return PartialResult;
  }
}

