/*
Methods for computing the solution to a linear system with coefficient matrix

   R = T + H

where T is symmetric Toeplitz and H is Hankel

Based on:

I. Gohberg and I. Koltracht
Efficient Algorithms for Toeplitz Plus Hankel Matrices
Integral Equations and Operator Theory
Vol. 12
1989

Tobias L. Jensen
March 2015
*/
#ifndef __TH__H_
#define __TH__H_

#include <stdlib.h>

#include "tools.hpp"
#include "vector.hpp"

/* Solves the system with f as right hand side. 
   t: length N+1
   h: lenght 2*N+1
   f: lenght N+1
   R is N+1 x N+1 
 
   returns the solution in x. 
*/
void solve(int N, FTYPE *t, FTYPE *h, FTYPE *f, FTYPE *gamma, FTYPE *x);
 
/* Computes the gamma values for use in solve
   t: length N+1
   h: lenght 2*N+1
   gamma: length (N+2)*(N+1)/2
          gamma_0: lenght 1
          gamma_1: length 2
          gamma_2: length 3
          
          In general,
          gamma_k: starts at (k+1)*k/2 and has length k+1
*/
void th(int N, FTYPE *t, FTYPE *h, FTYPE *gamma);


/* Computes the gamma values for use in solve with coefficient matrix
 
   R = T - H

   t: length N+1
   h: lenght 2*N+1 (should be called with -h)
   gamma: length (N+2)*(N+1)/2
          gamma_0: lenght 1
          gamma_1: length 2
          gamma_2: length 3
          
          In general,
          gamma_k: starts at (k+1)*k/2 and has length k+1

  Special case with t_2 = h_0, t_3 = h_1, etc
*/
void thp(int N, FTYPE *t, FTYPE *h, FTYPE *gamma);

#endif
