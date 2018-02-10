#include <stdlib.h>
#include <cstdio>
#include <cmath>
#include <sys/time.h>

#include "../src/single_pitch.hpp"

int main(int argc, char **argv){

  int L = 20; //controlled model order
  int ell = 4;
  int N = 160;
  int F = 5*N*L;
  FTYPE pitch_bounds[] = {0.001, 0.45};
  int i, j;
  struct timeval tim;
  double l1, l2;
  int rpt = 100;
  double eps = 1e-3;
  
  printf("Repetitions %d, Signal length N = %d, Max model order L = %d, eps (accuracy) = %1.1e\n", rpt, N, L, eps);
  /* construct solver object */
  gettimeofday(&tim, NULL);
  l1 = tim.tv_sec + tim.tv_usec*1e-6;
  single_pitch * sp = new single_pitch(L, F, N, pitch_bounds);
  gettimeofday(&tim, NULL);
  l2 = tim.tv_sec + tim.tv_usec*1e-6;
  printf("Solver constuction time : %2.5f [s]\n", (l2-l1));
 
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
  l1 = tim.tv_sec + tim.tv_usec*1e-6;

  for( j = 0 ; j < rpt ; ++j){
    sp->est(x, 1.0, eps);
  }

  gettimeofday(&tim, NULL);
  l2 = tim.tv_sec + tim.tv_usec*1e-6;

  printf("Time total with standard ordering     : %2.5f [s]\n", (l2-l1));
  printf("Time per solve with standard ordering : %2.5f [s]\n", (l2-l1)/rpt);


  /* repeat with refinement */
  FTYPE * omega_0h;
  
  gettimeofday(&tim, NULL);
  l1 = tim.tv_sec + tim.tv_usec*1e-6;

  for( j = 0 ; j < rpt ; ++j){
    sp->est_fast(x, 1.0, eps);
  }

  gettimeofday(&tim, NULL);
  l2 = tim.tv_sec + tim.tv_usec*1e-6;

  printf("Time total with permuted ordering (fast)     : %2.5f [s]\n", (l2-l1));
  printf("Time per solve with permuted ordering (fast) : %2.5f [s]\n", (l2-l1)/rpt);

  delete sp;
  delete [] x;
}


