/* 

  Various elementwise tools. Developed specifically for the
  use in single pitch ML estimation

  Tobias L Jensen, 
  March 2015
*/


#include "tools.hpp"

// Copy y <= x
void CvmdCopy(int n, FTYPE *x, FTYPE *y){

  #ifdef MKL
  COPY(n, x, 1, y, 1);
  #else
  for(int i = 0 ; i < n ; ++i)
    y[i] = x[i];
  #endif

}

// Scal x <= alpha*x
void CvmdScal(int n, FTYPE alpha, FTYPE *x){
  
  #ifdef MKL
  SCAL(n, alpha, x, 1);
  #else
  for(int i = 0 ; i < n ; ++i)
    x[i] = alpha*x[i];
  #endif
}

// returns  y <= alpha x + y
void CvmdAxpy(int n, FTYPE alpha, FTYPE *x, FTYPE *y){

  #ifdef MKL
  AXPY(n, alpha, x, 1, y, 1);
  #else
  for(int i = 0 ; i < n ; ++i)
    y[i] += alpha*x[i];
  #endif
}

// returns inner dot product a^T b
FTYPE CvmdDot(int n, FTYPE *a, FTYPE *b){

  #ifdef MKL
  return DOT(n, a, 1, b, 1);
  #else
  FTYPE c = 0.0;
  for(int i = 0 ; i < n ; ++i)
    c += a[i]*b[i];
  return c;
  #endif


}

// Adds c_i <- a + b_i
void CvmdAddConstant(int n, FTYPE a, FTYPE *b, FTYPE *c){

  for(int i = 0 ; i < n ; ++i)
    c[i] = a + b[i];
}

// Adds c <- a + b
void CvmdAdd(int n, FTYPE *a, FTYPE *b, FTYPE *c){

  #ifdef MKL
  VADD(n, a, b, c);
  #else
  for(int i = 0 ; i < n ; ++i)
    c[i] = a[i] + b[i];
  #endif
}

// Adds c <- a - b
void CvmdSub(int n, FTYPE *a, FTYPE *b, FTYPE *c){

  #ifdef MKL
  VSUB(n, a, b, c);
  #else  
  for(int i = 0 ; i < n ; ++i)
    c[i] = a[i] - b[i];
  #endif
}

// Adds a <- a + b
void CvmdAddInplace(int n, FTYPE *a, FTYPE *b){

  #ifdef MKL
  VADD(n, a, b, a);
  #else  
  for(int i = 0 ; i < n ; ++i)
    a[i] += b[i];
  #endif
}

// Adds a <- a - b
void CvmdSubInplace(int n, FTYPE *a, FTYPE *b){

  #ifdef MKL
  VSUB(n, a, b, a);
  #else    
  for(int i = 0 ; i < n ; ++i)
    a[i] -= b[i];
  #endif
}

//forms reverse vector b <- a_[n:1:-1]
void CvmdReverse(int n, FTYPE *a, FTYPE *b){ 

  for(int i=0 ; i < n ; ++i)
    b[i] = a[n-i-1];

}

//forms symmetric vector b <- [a_[1:n]; a_[n-1:1:-1]]
void CvmdSymmetric(int n, FTYPE *a, FTYPE *b){ 

  for(int i=0 ; i < n ; ++i)
    b[i] = a[i];

  for(int i=0 ; i < n-1 ; ++i)
    b[n+i] = a[n-i-2];
}

//computer c <- a_[1:n] + [0; b_[1:n-1]]
void CvmdShift(int n, FTYPE *a, FTYPE *b, FTYPE *c){ 
  c[0] = a[0];
  for(int i=1 ; i < n; ++i)
    c[i] = a[i] + b[i-1];
}

//Computes elementwise y_i <- a - x_i^2
void CvmdScalSquare(int n, FTYPE a, FTYPE * x, FTYPE * y){
  
  for(int i = 0 ; i < n ; ++i){
    y[i] = a - x[i]*x[i];
  }
}



//Initialize a vector with value a
void CvmdInit(int n, FTYPE a, FTYPE * x){
  
  for(int i = 0 ; i < n ; ++i){
    x[i] = a;
  }
}


//Computes elementwise sin y_i <- sin(a * x_i), where x_i is from min to max
void CvmdSinRange(int min, int max, FTYPE a, FTYPE * y){
  
  FTYPE * yp = y;
  for(int i = min ; i <= max ; ++i){
    *yp++ = sin(a*i);
  }
}

//Computes elementwise sin y_i <- sin(a * x_i), 
 void CvmdSinRange(int n, FTYPE * x, FTYPE a, FTYPE * t, FTYPE * y){

  #ifdef MKL
  CvmdInit(n, 0.0, t);
  AXPY(n, a, x, 1, t, 1);
  VSIN(n, t, y);
  #else
  for(int i = 0 ; i < n ; ++i){
    y[i] = sin(a*x[i]);
  }
  #endif
}


//Computes elementwise complex exponential y_i <- exp(1i*alpha * x_i)
// where x_i is from min to max
// and alpha and xi are real
// y_i is in interleaved complex format
void CvmdExpRange(int min, int max, FTYPE a, FTYPE * y){
  FTYPE v;
  FTYPE *yp = y;

  for(int i = min ; i <= max ; ++i){
    v = a*i;
    *yp++ = cos(v);
    *yp++ = sin(v);
  }
}

//Computes elementwise complex exponential y_i <- exp(1i*alpha * x_i)
// where alpha and xi are real
// y_i is in interleaved complex format
// t1 and t2 are temp vectors of length n
 void CvmdExpRange(int n, FTYPE * x, FTYPE a, 
    FTYPE * t1, FTYPE * t2, FTYPE * y){

  #ifdef MKL
  CvmdInit(n, 0.0, y);
  AXPY(n, a, x, 1, y, 1);
  VSINCOS(n, y, t1, t2);
  VUNPACKI(n, t2, y, 2);
  VUNPACKI(n, t1, y+1, 2);
  #else
  FTYPE v;
  FTYPE *yp = y;

  for(int i = 0 ; i < n ; ++i){
    v = a*x[i];
    *yp++ = cos(v);
    *yp++ = sin(v);
  }
  #endif
}

//Computes elementwise 
//   y_i = cos(alpha * x_i)
// and
//   z_i = sin(alpha * x_i)
// t is a temp vectors of length n
void CvmdCosSinRange(int n, FTYPE * x, FTYPE a, FTYPE *t,
                     FTYPE * y, FTYPE * z){

  #ifdef MKL
  CvmdInit(n, 0.0, t);
  AXPY(n, a, x, 1, t, 1);
  VSINCOS(n, t, z, y);
  #else
  FTYPE v;

  for(int i = 0 ; i < n ; ++i){
    v = a*x[i];
    y[i] = cos(v);
    z[i] = sin(v);
  }
  #endif
}


//Computes elementwise multiplication y_i <- x_i * z_i
void CvmdMul(int n, FTYPE * x, FTYPE * z, FTYPE * y){
  
  #ifdef MKL
  VMUL(n, x, z, y);
  #else
  for(int i = 0 ; i < n ; ++i){
    y[i] = x[i]*z[i];
  }
  #endif
}


//Computes elementwise multiplication v + iw <- x_i * z_(si) 
// with x, z complex and step s
// Real part goes into v, imaginary part into w
// t is a temp vector with n complex number
void CvmdMulz(int n, FTYPE * x, FTYPE * z, int s, FTYPE *t,
              FTYPE * v, FTYPE * w){

  #ifdef MKL
  vzPackI(n, (MKL_COMPLEX*)z, s, (MKL_COMPLEX*)t);
  vzMul(n, (MKL_COMPLEX*)x, (MKL_COMPLEX*)t, (MKL_COMPLEX*)t);
  vdPackI(n, t, 2, v);
  vdPackI(n, t+1, 2, w);
  #else
  int ii = 0;
  for(int i = 0 ; i < n ; ++i){
    v[i] = x[2*i]*z[ii] - x[2*i+1]*z[ii+1];
    w[i] = x[2*i]*z[ii+1] + x[2*i+1]*z[ii];
    ii += 2*s;
  }
  #endif
}

//Computes elementwise multiplication v + iw <- x_i * z_i 
// with x, z complex
// Real part goes into v, imaginary part into w
// t is a temp vector with n complex number
void CvmdMulz(int n, FTYPE * x, FTYPE * z, FTYPE *t,
             FTYPE * v, FTYPE * w){

  #ifdef MKL
  vzMul(n, (MKL_COMPLEX*)x, (MKL_COMPLEX*)z, (MKL_COMPLEX*)t);
  vdPackI(n, t, 2, v);
  vdPackI(n, t+1, 2, w);
  #else
  int ii = 0;
  for(int i = 0 ; i < n ; ++i){
    v[i] = x[2*i]*z[ii] - x[2*i+1]*z[ii+1];
    w[i] = x[2*i]*z[ii+1] + x[2*i+1]*z[ii];
    ii += 2;
  }
  #endif
}



//Computes elementwise division y_i <- x_i / z_i
void CvmdDiv(int n, FTYPE * x, FTYPE * z, FTYPE * y){
  
  #ifdef MKL
  VDIV(n, x, z, y);
  #else
  for(int i = 0 ; i < n ; ++i){
    y[i] = x[i]/z[i];
  }
  #endif
}

//Computes elementwise inverse y_i <- 1.0 / x_i
void CvmdInverse(int n, FTYPE * x, FTYPE * y){
  
  #ifdef MKL
  VINV(n, x, y);
  #else
  for(int i = 0 ; i < n ; ++i){
    y[i] = 1/x[i];
  }
  #endif
}

// packs x from step s to unit stride in y
void CvmdPack(int n, FTYPE * x, int s, FTYPE * y){
  int ii = 0;
  for( int i = 0 ; i < n ; ++i){
    y[i] = x[ii];
    ii += s;
  }
}


int CvmdArgmax(int n, FTYPE * x){

  FTYPE max;
  int arg;

  max = x[0];
  arg = 0;

  for( int k = 1 ; k < n ; ++k){
    if( x[k] > max){
      max = x[k];
      arg =  k;
    }
  }

  return arg;

}
