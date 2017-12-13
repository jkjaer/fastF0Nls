#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <sys/time.h>

#include "../src/single_pitch.hpp"

int main(int argc, char **argv){

  int L = 20; //controlled model order
  int ell = 4;
  int N = 1000;
  int F = 5*N*L;
  FTYPE pitch_bounds[] = {0.001, 0.45};
  int i, j;
  struct timeval tim;
  double l1, l2;
  int rpt = 10;

  /* construct solver object */
  gettimeofday(&tim, NULL);
  l1 = tim.tv_sec + (tim.tv_usec/1000000.0);
  single_pitch * sp = new single_pitch(L, F, N, pitch_bounds);
  gettimeofday(&tim, NULL);
  l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
 
  /* construct signal */
  FTYPE * x = new FTYPE [N];
  FTYPE omega = 0.1;

  for( j = 0 ; j < N ; ++j ){
    x[j] = 0;
    for( i = 1 ; i <= ell ; ++i ){
      x[j] += sin(M_PI*i*omega*j);
    }
  }

  /* repeat */
  gettimeofday(&tim, NULL);
  l1 = tim.tv_sec + (tim.tv_usec/1000000.0);

  for( j = 0 ; j < rpt ; ++j){
    sp->compute(x);
  }

  gettimeofday(&tim, NULL);
  l2 = tim.tv_sec + (tim.tv_usec/1000000.0);


  /* repeat with refinement */
  FTYPE * omega_0h;
  
  gettimeofday(&tim, NULL);
  l1 = tim.tv_sec + (tim.tv_usec/1000000.0);

  for( j = 0 ; j < rpt ; ++j){
    sp->compute(x);
    omega_0h = sp->refine(x);
  }

  gettimeofday(&tim, NULL);
  l2 = tim.tv_sec + (tim.tv_usec/1000000.0);


  delete sp;
  delete [] x;
}


