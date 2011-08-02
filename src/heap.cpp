# include "util.h"
# include "heap.h"


void inline Swap(int &a, int &b)
{
  int x = a;
  a = b;
  b = x;
}

bool IsIn(int *T, int TS, int x)
{
  for (int i = 0; i < TS; i++)
    if (T[i] == x)
      return true;
  return false;
}

Heap::Heap(const int **L, const int *TL, int Nb)
{
  Lists = L;
  ListsSizes = TL;
  NbLists = Nb;
  IndexInList = new int[NbLists];
  for (int i = 0; i < NbLists; i++)
    IndexInList[i] = 0;
  MyHeap = new int[NbLists];
  for (int i = 0; i < NbLists; i++)
    MyHeap[i] = i;
  HeapSize = 0;
  ComputationFinished = false;
  for (int i = 0; i < NbLists; i++)
    AddNode();
}

/* Michel
Heap::Heap(int **L, int *TL, int Nb)
{
  Lists = L;
  ListsSizes = TL;
  NbLists = Nb;
  IndexInList = new int[NbLists];
  for (int i = 0; i < NbLists; i++)
    IndexInList[i] = 0;
  MyHeap = new int[NbLists];
  for (int i = 0; i < NbLists; i++)
    MyHeap[i] = i;
  HeapSize = 0;
  ComputationFinished = false;
  for (int i = 0; i < NbLists; i++)
    AddNode();
}
*/


Heap::~Heap()
{
  delete[] IndexInList;
  delete[] MyHeap;
}

void Heap::AddNode()
{
  int CurIndex = HeapSize;
  while ((CurIndex > 0) && (Lists[MyHeap[CurIndex]][0] < Lists[MyHeap[(CurIndex - 1) / 2]][0]))
  {
    Swap(MyHeap[CurIndex], MyHeap[(CurIndex - 1) / 2]);
    CurIndex = (CurIndex - 1) / 2;
  }
  HeapSize++;
}

bool Heap::AllElementsEqual()
{
  // Debug();
  for (int64_t i = HeapSize - 1; i >= HeapSize / 2; i--) {
    // Trace
    // cout << "AllElement i =" << i << endl;
    // cout << "  MyHeap[i]="<< MyHeap[i] << endl;
    // cout << "  MyHeap[0]="<< MyHeap[0] << endl;
    // cout << "  IndexInList[MyHeap[i]]="<< IndexInList[MyHeap[i]] << endl;
    // cout << "  IndexInList[MyHeap[0]]="<< IndexInList[MyHeap[0]] << endl;
    // cout << "  Lists[MyHeap[i]][IndexInList[MyHeap[i]]]="<<  Lists[MyHeap[i]][IndexInList[MyHeap[i]]] << endl;
    // cout << "  Lists[MyHeap[0]][IndexInList[MyHeap[0]]]="<<  Lists[MyHeap[0]][IndexInList[MyHeap[0]]] << endl;

    //    if ( ( IndexInList[MyHeap[i]] < ListsSizes[MyHeap[i]] ) && ( IndexInList[MyHeap[0]] < ListsSizes[MyHeap[0]] ) ) { 
    if ( ( IndexInList[MyHeap[i]] < ListsSizes[MyHeap[i]] ) ) { 
      if ( Lists[MyHeap[i]][IndexInList[MyHeap[i]]] 
	   != Lists[MyHeap[0]][IndexInList[MyHeap[0]]]) {
	return false;
      }
    } else {
      return false;
    }
  }
  return true;
}

void Heap::RemoveHead()
{
  if (ComputationFinished)
    return;
  IndexInList[MyHeap[0]]++;
  
  // Trace
  // cout << " RemoveHeap : " << IndexInList[MyHeap[0]] << "==" << ListsSizes[MyHeap[0]] << endl;

  if (IndexInList[MyHeap[0]] == ListsSizes[MyHeap[0]] )
  {
    ComputationFinished = true;
    return;
  }
  int CurIndex = 0;
  if (2 * CurIndex + 1 < HeapSize)
  {
    while (2 * CurIndex + 2 < HeapSize)
    {
      int MinIndex = 2 * CurIndex + 1;
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


int Heap::NextCommonElement()
{
  if (ComputationFinished)
    return NOT_INDEX;
  while ((!AllElementsEqual()) && (!ComputationFinished))
    RemoveHead();
  if (!ComputationFinished)
  {
    int Res = Lists[MyHeap[0]][IndexInList[MyHeap[0]]];
    for (int i = 0; i < NbLists; i++)
      RemoveHead();
    return Res;
  }
  return NOT_INDEX;
}


void Heap::Debug()
{
  cout << "Heap status : " << NbLists << " lists" << endl;
  cout << "  Lists sizes: ";
  for (int i = 0; i < NbLists; i++)
    cout  << ListsSizes[i] << " ";
  cout << endl;
  cout << "  Contents of lists: " << endl;
  for (int i = 0; i < NbLists; i++) {
    cout << "    list " << i << " :"; 
      for (int j = 0; j <  ListsSizes[i] ; j++)
	cout  << (Lists[i])[j] << " ";
      cout << endl;
  }
  cout << "  Index in lists :";
  for (int i = 0; i < NbLists; i++)
    cout << IndexInList[i] << " ";
  cout << endl;
  cout << "  Heap size :" << HeapSize << endl;

}
int *Heap::Union(int &CardinalUnion)
{
  CardinalUnion = 0;
  for (int i = 0; i < NbLists; i++)
    CardinalUnion += ListsSizes[i];
  int *Result = new int[CardinalUnion];
  int Index = 0;
  while (HeapSize > 0)
  {
    int CurrentElement = Lists[MyHeap[0]][IndexInList[MyHeap[0]]];
    do RemoveHead(); while ((HeapSize > 0) && (CurrentElement == Lists[MyHeap[0]][IndexInList[MyHeap[0]]]));
    Result[Index++] = CurrentElement;
  }
  CardinalUnion = Index;
  return Result;
}

void Heap::Union(int *Result, int &CardinalUnion)
{
  CardinalUnion = 0;
  while (HeapSize > 0)
  {
    int CurrentElement = Lists[MyHeap[0]][IndexInList[MyHeap[0]]];
    do RemoveHead(); while ((HeapSize > 0) && (CurrentElement == Lists[MyHeap[0]][IndexInList[MyHeap[0]]]));
    Result[CardinalUnion++] = CurrentElement;
  }
}



