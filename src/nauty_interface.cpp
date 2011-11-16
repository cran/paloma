# include "nauty_interface.h"

static DEFAULTOPTIONS(options);
static statsblk stats;


// Get the canonic form of the adjacency matrix m (n vertices)
// Return the the canonic form and the permutation to obtain it.
// The storage of the reurn values must be provided.  
void getCanonic(const int* m, const size_t n, const int directed, int *canon, int *reorder) {

  int *ptn    = new int[n];
  int *orbits = new int[n];
  int worksize = 50*n;
  setword *workspace = new setword[worksize];
  
  options.writeautoms = FALSE;
  options.getcanon    = TRUE;

  // Directed graph
  options.digraph = (directed != 0);

  // Defaut partition
  options.defaultptn = TRUE;

  int N = n; 
  int M = n/WORDSIZE;
  M = M + ((n % WORDSIZE) != 0);
  graph *g          = new setword[N*M];
  graph *canon_word = new setword[N*M];
  setword w    = 0;

  // Put the adjacency matrix in bit mode storage.
  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      w = 0;
      for(int l=0; l<MIN(WORDSIZE, N); l++) {
	// Transpose
	w = w | (m[i+(j*WORDSIZE+l)*n] << (WORDSIZE-1-l));
	// w = w | (m[j+(i*WORDSIZE+l)*n] << (WORDSIZE-1-l));
      }
      g[i+N*j] = w;
    } 
  } 

  nauty(g, reorder, ptn, NULL, orbits, &options, &stats, 
			  workspace, worksize, M, N, canon_word);
  int remain;
  setword bit_ = 0;

  // Inverse transform : bit storage to integer storage
  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      w = 0;
      bit_ = 0;

      w = canon_word[i+j*M];
      remain = N - j*WORDSIZE;
      for(int l=0; l<MIN(WORDSIZE, remain); l++) {
	bit_ =   ( (w & (1 << (31-l))) != 0) ? 1 : 0;
	canon[i + ((j*WORDSIZE)+l)*N] = bit_;
      } 
    } 
  }
  // ???  if (reorder != NULL ) {
  // 
  //  for(int i=0; i<N; i++) {
  //    reorder[i] = lab[i];
  //  }
  // }
  delete [] g;
  delete [] canon_word;
  delete [] ptn;
  delete [] orbits;
  delete [] workspace;

  // To avoid warnings

# ifdef MSG 
  LOGMSG(100, std::cout, "NULL message ", "" );
# endif

}
