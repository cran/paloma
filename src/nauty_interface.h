
# include "util.h"

extern "C" {
   #include "nauty.h"
   #include "nautinv.h"
 
}

void getCanonic( const int* m, const size_t n, const int directed, 
		   int *canon, int *reorder);


