/*
 Defines usefull constructs (memory constructs) and
  free that supports several memory organization approaches
  using proper usage of defines

  MKL:  uses MKL malloc and free. 
    Define ALIGNMENT32, ALIGNMENT64, ALIGNMENT128, as well
    to define the alignment. 

  else: uses standard C++ allocation via new and delete of arrays.

*/

#include "vector.hpp"

FTYPE * vector(int N){
  FTYPE * x;
  
  #ifdef MKL
  x = (FTYPE*)mkl_malloc(sizeof(FTYPE)*N, ALIGNMENT);
  #else
  x = new FTYPE [N];
  #endif

  if( x == NULL ){
    printf("Vector.cpp: Could not allocate %d bytes of memory", 
           (int)sizeof(FTYPE)*N);
    exit(EXIT_FAILURE);
  }

  return x;
}


void del_vector(FTYPE * x){

  #ifdef MKL
  mkl_free(x);
  #else
  delete [] x;
  #endif
  
}
