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
#include "th.hpp"


void solve(int N, FTYPE *t, FTYPE *h, FTYPE *f, FTYPE *gam, FTYPE *x){
  FTYPE *tp, *t1, lambda_k, zero = 0.0, *gamma;
  int K = ((N+2)*(N+1))/2;

  tp = vector(N+1);
  t1 = vector(N+1);

  if(gam==NULL){
    gamma = vector(K);
    th(N, t, h, gamma);
  }
  else
    gamma = gam;

  
  CvmdReverse(N+1, t, tp);
  CvmdInit(N+1, zero, x);
  
  x[0] = f[0]/(t[0] + h[0]);
  
  for(int k = 1 ; k < N+1 ; ++k){
    CvmdAdd(k, tp+N-k, h+k, t1);
    lambda_k = f[k] - CvmdDot(k, t1, x);

    /*  x[:k+1] += lambda_k*gamma[k] */
    CvmdAxpy(k+1, lambda_k, gamma+((k+1)*k)/2, x);
  }

  del_vector(tp);
  del_vector(t1);

  if(gam == NULL)
    del_vector(gamma);
}


void th(int N, FTYPE *t, FTYPE *h, FTYPE *gamma){

  FTYPE *a, *tp, *alpha, *phi_k, *psi_k, *t1, *r_kp1, *b_k;
  FTYPE zero = 0.0, R00, R01, R11, s, lambda_k, mu_k, nu_kp1, gammap;
  int K, Kp1;

  // test two special cases 
  if( N == 0 ){
    R00 = t[0] + h[0];
    gamma[0] = 1/R00;
    return;
  }
  if( N == 1 ){    
    R00 = t[0] + h[0];
    R01 = t[1] + h[1];
    R11 = t[0] + h[2];
  
    gamma[0] = 1/R00;

    s = 1/(R00*R11 - R01*R01);
    gamma[1] = -s*R01;
    gamma[2] = s*R00;
    return;
  }

  /* allocate */
  a= vector(N); 
  tp = vector(N+1);
  alpha = vector(N+1);
  phi_k = vector(N+1);
  psi_k = vector(N+1);
  t1 = vector(N+1);
  r_kp1 = vector(N+1);
  b_k = vector(N+1);

  CvmdShift(N, t+1, h, a);
  CvmdReverse(N+1, t, tp);

  CvmdInit(N+1, zero, phi_k);
  CvmdInit(N+1, zero, psi_k);
  CvmdInit(N+1, zero, alpha);

  R00 = t[0] + h[0];
  R01 = t[1] + h[1];
  R11 = t[0] + h[2];
  psi_k[0] = 1/R00;
  
  gamma[0] = 1/R00;
  phi_k[0] = a[0]/R00;
  alpha[1] = R01/R00;

  s = 1/(R00*R11 - R01*R01);
  gamma[1] = -s*R01;
  gamma[2] = s*R00;

  Kp1 = 1;

  for(int k = 1 ; k < N ; ++k){

    K = Kp1;
    Kp1 = (k+2)*(k+1)/2;

    /* lambda_k = a[k] - np.dot(tp[N-k:N]+h[k:2*k], phi_k[:k]) */
    CvmdAdd(k, tp+N-k, h+k, t1);
    lambda_k = a[k] - CvmdDot(k, t1, phi_k);
    
    /* mu_k = - np.dot(tp[N-k:N]+h[k:2*k], psi_k[:k]) */
    mu_k = -CvmdDot(k, t1, psi_k);

    /* phi_k[:k+1] += lambda_k*gamma[k] */
    CvmdAxpy(k+1, lambda_k, gamma+K, phi_k);

    /* psi_k[:k+1] += mu_k*gamma[k] */
    CvmdAxpy(k+1, mu_k, gamma+K, psi_k);

    /*r_kp1 = tp[N-k-1:N] + h[k+1:2*k+2]*/
    CvmdAdd(k+1, tp + N - k - 1, h + k + 1, r_kp1);
    
    /* alpha[k+1] = np.dot(tp[N-k-1:N]+h[k+1:2*(k+1)], gamma[k]) */
    alpha[k+1] = CvmdDot(k+1, r_kp1, gamma + K);

    /* b_k = ((alpha[k] - alpha[k+1])*gamma[k] + np.hstack([0, gamma[k][:k]])
               + np.hstack([gamma[k][1:k+1], 0])
               -np.hstack([gamma[k-1], 0])
               + psi_k[k]*phi_k[:k+1] - phi_k[k]*psi_k[:k+1]) */
    CvmdInit(k+1, 0.0, b_k);
    CvmdAxpy(k+1, alpha[k]-alpha[k+1], gamma + K, b_k);
    CvmdAxpy(k, 1.0, gamma + K, b_k + 1);
    CvmdAxpy(k, 1.0, gamma + K + 1, b_k);
    CvmdAxpy(k, -1.0, gamma + ((k)*(k-1))/2, b_k);
    CvmdAxpy(k, psi_k[k], phi_k, b_k);
    CvmdAxpy(k, -phi_k[k], psi_k, b_k);

    /* nu_kp1 = (1/gamma[k][k]) * np.dot(tp[N-k-1:N]+h[k+1:2*(k+1)], b_k) */
    nu_kp1 = (1/gamma[K + k]) 
      * CvmdDot(k+1, r_kp1, b_k);

    gammap = 1.0/(nu_kp1 + tp[N] + h[2*k+2]);

    /* gamma.append(np.hstack([(gammap/gamma[k][k])*b_k, gammap])) */
    CvmdScal(k+1, gammap/gamma[K + k], b_k);
    CvmdCopy(k+1, b_k, gamma + Kp1);
    gamma[Kp1 + k + 1] = gammap;
  }

  /* de-allocate */
  del_vector(a);
  del_vector(tp);
  del_vector(alpha);
  del_vector(phi_k);
  del_vector(psi_k);
  del_vector(t1);
  del_vector(r_kp1);
  del_vector(b_k);
}


void thp(int N, FTYPE *t, FTYPE *h, FTYPE *gamma){

  FTYPE *tp, *alpha, *r_kp1, *b_k;
  FTYPE zero = 0.0, R00, R01, R11, s, nu_kp1, gammap;
  int K, Kp1;

  // test two special cases 
  if( N == 0 ){
    R00 = t[0] + h[0];
    gamma[0] = 1/R00;
    return;
  }
  if( N == 1 ){    
    R00 = t[0] + h[0];
    R01 = t[1] + h[1];
    R11 = t[0] + h[2];
  
    gamma[0] = 1/R00;

    s = 1/(R00*R11 - R01*R01);
    gamma[1] = -s*R01;
    gamma[2] = s*R00;
    return;
  }

    
  /* allocate */
  tp = vector(N+1);
  alpha = vector(N+1);
  r_kp1 = vector(N+1);
  b_k = vector(N+1);

  CvmdReverse(N+1, t, tp);

  CvmdInit(N+1, zero, alpha);

  R00 = t[0] + h[0];
  R01 = t[1] + h[1];
  R11 = t[0] + h[2];
  
  gamma[0] = 1/R00;
  alpha[1] = R01/R00;

  s = 1/(R00*R11 - R01*R01);
  gamma[1] = -s*R01;
  gamma[2] = s*R00;

  Kp1 = 1;

  for(int k = 1 ; k < N ; ++k){

    K = Kp1;
    Kp1 = (k+2)*(k+1)/2;

    /*r_kp1 = tp[N-k-1:N] + h[k+1:2*k+2]*/
    CvmdAdd(k+1, tp + N - k - 1, h + k + 1, r_kp1);
    
    /* alpha[k+1] = np.dot(tp[N-k-1:N]+h[k+1:2*(k+1)], gamma[k]) */
    alpha[k+1] = CvmdDot(k+1, r_kp1, gamma + ((k+1)*k)/2);

    /* b_k = ((alpha[k] - alpha[k+1])*gamma[k] + np.hstack([0, gamma[k][:k]])
               + np.hstack([gamma[k][1:k+1], 0])
               -np.hstack([gamma[k-1], 0])
                */
    CvmdInit(k+1, 0.0, b_k);
    CvmdAxpy(k+1, alpha[k]-alpha[k+1], gamma + ((k+1)*k)/2, b_k);
    CvmdAxpy(k, 1.0, gamma + K, b_k + 1);
    CvmdAxpy(k, 1.0, gamma + K + 1, b_k);
    CvmdAxpy(k, -1.0, gamma + ((k)*(k-1))/2, b_k);

    /* nu_kp1 = (1/gamma[k][k]) * np.dot(tp[N-k-1:N]+h[k+1:2*(k+1)], b_k) */
    nu_kp1 = (1/gamma[K + k]) 
      * CvmdDot(k+1, r_kp1, b_k);

    gammap = 1.0/(nu_kp1 + tp[N] + h[2*k+2]);

    /* gamma.append(np.hstack([(gammap/gamma[k][k])*b_k, gammap])) */
    CvmdScal(k+1, gammap/gamma[K + k], b_k);
    CvmdCopy(k+1, b_k, gamma + Kp1);
    gamma[Kp1 + k + 1] = gammap;
  }

  /* de-allocate */
  del_vector(tp);
  del_vector(alpha);
  del_vector(r_kp1);
  del_vector(b_k);
}
