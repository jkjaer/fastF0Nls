#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <stdio.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "gtest/gtest.h"

#include "../src/single_pitch.hpp"

#ifdef DOUBLE
#define EPS 1e-10 // pure floating
#define EPS2 1e-8 //iterative algs
#else
#define EPS 1e-1
#define EPS2 1e-1
#endif

class RealDataTest : public ::testing::Test{
protected:

  RealDataTest(){
  }

  virtual ~RealDataTest(){
  }

  virtual void SetUp(){
  }

  virtual void TearDown(){
  }
};

TEST_F(RealDataTest, ComputeGamma1) {

  hid_t file;
  herr_t status;


  file = H5Fopen("data_files/unittest1.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  int Mp = F/2+1;
  char label[50];
  
  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif

  double * cost_function = new double [L*Mp];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  /* construct object and that will compute gamma1 and gamma2  */
  single_pitch * sp = new single_pitch(L, F, N, pb);


  /* Tjeck that gamma is computed correctly */
  double * gammah;

  for(int ell = 1 ; ell <= L ; ell++){

    gammah = new double[ell*(int)dim[ell-1]];

    /* test for gamma_1 */
    sprintf(label, "/gamma_1_%d", ell);
    H5LTread_dataset(file, label, 
                     H5T_NATIVE_DOUBLE, gammah);

    for(int k = 0 ; k < ell*dim[ell-1] ; ++k ){
      ASSERT_NEAR(gammah[k], sp->Gamma1[k + Mp*(ell-1)], EPS) 
        << "Gamma1: Vectors differ at index " << k << ", gamma level " << ell;
    }

    /* test for gamma_2 */
    sprintf(label, "/gamma_2_%d", ell);
    H5LTread_dataset(file, label, 
                     H5T_NATIVE_DOUBLE, gammah);

    for(int k = 0 ; k < ell*dim[ell-1] ; ++k ){
      ASSERT_NEAR(gammah[k], sp->Gamma2[k + Mp*(ell-1)], EPS) 
        << "Gamma2: Vectors differ at index " << k << ", gamma level " << ell;
    }

    delete [] gammah;

  }
  
  /* close file */
  status = H5Fclose(file);

  del_vector(x);
  delete [] cost_function;
  delete [] dim;
  delete sp;

}

TEST_F(RealDataTest, ComputeCost1) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest1.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  char label[50];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif


  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  int M = (int) dim[0];

  double * cost_function = new double [L*M];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  status = H5Fclose(file);


  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  del_vector(xs);
  #else
  sp->compute(x);
  #endif

  /* Tjeck that the cost function is correctly computed */
  for(int ell = 1 ; ell <= L ; ell++){

    for(int k = 0 ; k < dim[ell-1] ; ++k ){
      ASSERT_NEAR(cost_function[M*(ell-1)+k], 
                  sp->costFunctionMatrix[(ell-1)][k], EPS) 
        << "Cost function differ at index " << k << ", level " << ell;
    }
  }


  
  delete sp;
  del_vector(x);
  delete [] cost_function;
  delete [] dim;
}


TEST_F(RealDataTest, ComputeGamma2) {

  hid_t file;
  herr_t status;


  file = H5Fopen("data_files/unittest2.h5", H5F_ACC_RDONLY, H5P_DEFAULT);

  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  int Mp = F/2+1;
  char label[50];
  
  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif

  double * cost_function = new double [L*Mp];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  /* construct object and that will compute gamma1 and gamma2  */
  single_pitch * sp = new single_pitch(L, F, N, pb);


  /* Tjeck that gamma is computed correctly */
  double * gammah;

  for(int ell = 1 ; ell <= L ; ell++){

    gammah = new double[ell*(int)dim[ell-1]];

    /* test for gamma_1 */
    sprintf(label, "/gamma_1_%d", ell);
    H5LTread_dataset(file, label, 
                     H5T_NATIVE_DOUBLE, gammah);

    for(int k = 0 ; k < ell*dim[ell-1] ; ++k ){
      ASSERT_NEAR(gammah[k], sp->Gamma1[k + Mp*(ell-1)], EPS) 
        << "Gamma1: Vectors differ at index " << k << ", gamma level " << ell;
    }

    /* test for gamma_2 */
    sprintf(label, "/gamma_2_%d", ell);
    H5LTread_dataset(file, label, 
                     H5T_NATIVE_DOUBLE, gammah);

    for(int k = 0 ; k < ell*dim[ell-1] ; ++k ){
      ASSERT_NEAR(gammah[k], sp->Gamma2[k + Mp*(ell-1)], EPS) 
        << "Gamma2: Vectors differ at index " << k << ", gamma level " << ell;
    }

    delete [] gammah;

  }
  
  /* close file */
  status = H5Fclose(file);

  del_vector(x);
  delete [] cost_function;
  delete [] dim;
  delete sp;

}

TEST_F(RealDataTest, ComputeCost2) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest2.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  char label[50];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif


  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  int M = (int) dim[0];

  double * cost_function = new double [L*M];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  status = H5Fclose(file);


  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  del_vector(xs);
  #else
  sp->compute(x);
  #endif

  /* Tjeck that the cost function is correctly computed */
  for(int ell = 1 ; ell <= L ; ell++){

    for(int k = 0 ; k < dim[ell-1] ; ++k ){
      ASSERT_NEAR(cost_function[M*(ell-1)+k], 
                  sp->costFunctionMatrix[(ell-1)][k], EPS) 
        << "Cost function differ at index " << k << ", level " << ell;
    }
  }


  
  delete sp;
  del_vector(x);
  delete [] cost_function;
  delete [] dim;
}


TEST_F(RealDataTest, ComputeSingleCost3) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest3.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  char label[50];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif


  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  int M = (int) dim[0];

  double * cost_function = new double [L*M];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  status = H5Fclose(file);


  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  del_vector(xs);
  #else
  sp->compute(x);
  #endif

  /* Tjeck that the cost function is correctly computed */
  for(int ell = 1 ; ell <= L ; ell++){

    for(int k = 0 ; k < dim[ell-1] ; ++k ){
      ASSERT_NEAR(cost_function[M*(ell-1)+k], 
                  sp->costFunctionMatrix[(ell-1)][k], EPS) 
        << "Cost function differ at index " << k << ", level " << ell;
    }
  }

  int ell = 2;

  FTYPE max = *std::max_element(cost_function + M*(ell-1), 
                               cost_function + M*(ell-1) + sp->nPitches(ell)-1);
  FTYPE maxh = sp->max_obj(ell);
  
  ASSERT_NEAR(max, maxh, EPS)  << "max differs";
  

  FTYPE omega = sp->argmax_obj(ell);
  FTYPE maxhh = sp->compute_obj(omega, x, ell);

  ASSERT_NEAR(maxhh, max, EPS)  << "max differs from computed objective";

  delete sp;
  del_vector(x);
  delete [] cost_function;
  delete [] dim;
}


TEST_F(RealDataTest, ComputeSingleCost4) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest4.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  char label[50];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif


  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  int M = (int) dim[0];

  double * cost_function = new double [L*M];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  status = H5Fclose(file);


  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  del_vector(xs)
  #else
  sp->compute(x);
  #endif

  /* Tjeck that the cost function is correctly computed */
  for(int ell = 1 ; ell <= L ; ell++){

    for(int k = 0 ; k < dim[ell-1] ; ++k ){
      ASSERT_NEAR(cost_function[M*(ell-1)+k], 
                  sp->costFunctionMatrix[(ell-1)][k], EPS) 
        << "Cost function differ at index " << k << ", level " << ell;
    }
  }

  int ell = 4;

  FTYPE max = *std::max_element(cost_function + M*(ell-1), 
                               cost_function + M*(ell-1) + sp->nPitches(ell)-1);
  FTYPE maxh = sp->max_obj(ell);
  
  ASSERT_NEAR(max, maxh, EPS)  << "max differs";
  

  FTYPE omega = sp->argmax_obj(ell);
  FTYPE maxhh = sp->compute_obj(omega, x, ell);

  ASSERT_NEAR(maxhh, max, EPS)  << "max differs from computed objective";

  delete sp;
  del_vector(x);
  delete [] cost_function;
  delete [] dim;
}


TEST_F(RealDataTest, Refine4) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest5.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  //H5LTread_dataset(file, "/omega_0", H5T_NATIVE_DOUBLE, d);
  //double omega_0 =  d[0];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);

  double * omega_0h = new double [L];
  H5LTread_dataset(file, "/omega_0h", H5T_NATIVE_DOUBLE, omega_0h);

  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif


  double * dim = new double [L];
  H5LTread_dataset(file, "/dim", H5T_NATIVE_DOUBLE, dim);

  int M = (int) dim[0];

  double * cost_function = new double [L*M];
  H5LTread_dataset(file, "/costFunction", H5T_NATIVE_DOUBLE, cost_function);

  status = H5Fclose(file);


  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  del_vector(xs);
  FTYPE * omega_0hh = sp->refine(xs, 1e-8);
  #else
  sp->compute(x);
  FTYPE * omega_0hh = sp->refine(x, 1e-8);
  #endif

  
  for( int ell = 1 ; ell < L ; ++ell)
    
    ASSERT_NEAR(omega_0h[ell-1], omega_0hh[ell-1], EPS2)  
      << "refinement differs at iOrder = " << ell ;

  delete sp;
  del_vector(x);
  delete [] cost_function;
  delete [] dim;
  delete [] omega_0h;

}


TEST_F(RealDataTest, ModelOrderSelectionNoNoise) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest6.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  H5LTread_dataset(file, "/modelOrder", H5T_NATIVE_DOUBLE, d);
  int modelOrder = (int) d[0];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);


  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif

  status = H5Fclose(file);


  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  del_vector(xs)
  FTYPE * omega_0hh = sp->refine(xs, 1e-8);
  #else
  sp->compute(x);
  FTYPE * omega_0hh = sp->refine(x, 1e-8);
  #endif

  int estModelOrder = sp->model_order_selection(x);


  // Known values
  ASSERT_NEAR(22.971993868687033, sp->lnBF[1], EPS2);
  ASSERT_NEAR(22.4532605908562, sp->lnBF[2], EPS2);

  // Without noise the method should estimate the model order correctly
  ASSERT_EQ(modelOrder, estModelOrder);

  delete sp;
  del_vector(x);

}


TEST_F(RealDataTest, ModelOrderSelectionNoise) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest7.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  H5LTread_dataset(file, "/modelOrder", H5T_NATIVE_DOUBLE, d);
  int modelOrder = (int) d[0];

  
  H5LTread_dataset(file, "/pitch", H5T_NATIVE_DOUBLE, d);
  double pitch = d[0];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);


  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif

  status = H5Fclose(file);

  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  sp->refine(xs, 1e-8);
  del_vector(xs);
  #else
  sp->compute(x);
  sp->refine(x, 1e-8);
  #endif

  int estModelOrder = sp->model_order_selection(x);

  // Control the Bayes Factors
  
  FTYPE lnBF[] = {22.976505179621622, 22.444499748782523, 423.7940152029966, 
                  412.11543962457455, 400.92575920553287};

  for( int k = 1 ; k <= L ; ++k )
    ASSERT_NEAR(lnBF[k-1], sp->lnBF[k], EPS2)  
      << "Bayes Factor differs at order = " << k ;

  // With this realization the model order estimation should work
  ASSERT_EQ(modelOrder, estModelOrder);

  FTYPE omega_0h = sp->est(x);
  ASSERT_NEAR(omega_0h/(2*M_PI), pitch, 1e-6);  

  delete sp;
  del_vector(x);
}

TEST_F(RealDataTest, ModelOrderSelectionNoiseBayes) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest7.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  H5LTread_dataset(file, "/modelOrder", H5T_NATIVE_DOUBLE, d);
  int modelOrder = (int) d[0];

  H5LTread_dataset(file, "/pitch", H5T_NATIVE_DOUBLE, d);
  double pitch = d[0];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);


  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif

  status = H5Fclose(file);

  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  sp->refine(xs, 1e-8);
  del_vector(xs);
  #else
  sp->compute(x);
  sp->refine(x, 1e-8);
  #endif

  int estModelOrder = sp->model_order_selection(x);

  // Control the Bayes Factors
  
  FTYPE lnBF[] = {22.976505179621622, 22.444499748782523, 423.7940152029966, 
                  412.11543962457455, 400.92575920553287};

  for( int k = 1 ; k <= L ; ++k )
    ASSERT_NEAR(lnBF[k-1], sp->lnBF[k], EPS2)  
      << "Bayes Factor differs at order = " << k ;

  // With this realization the model order estimation should work
  ASSERT_EQ(modelOrder, estModelOrder);

  FTYPE omega_0h = sp->est(x, 20.0, 1e-6);
  estModelOrder = sp->modelOrder();
  ASSERT_NEAR(omega_0h/(2*M_PI), pitch, 1e-6);
  ASSERT_EQ(modelOrder, estModelOrder);

  delete sp;
  del_vector(x);
}

TEST_F(RealDataTest, ModelOrderSelectionNoiseBayesRefinementNewOrder) {

  hid_t file;
  herr_t status;

  /* extract data and output */
  file = H5Fopen("data_files/unittest7.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
  
  double d[1];
  H5LTread_dataset(file, "/nData", H5T_NATIVE_DOUBLE, d);
  int N = (int) d[0];

  H5LTread_dataset(file, "/L", H5T_NATIVE_DOUBLE, d);
  int L = (int) d[0];

  H5LTread_dataset(file, "/F", H5T_NATIVE_DOUBLE, d);
  int F = (int) d[0];

  H5LTread_dataset(file, "/modelOrder", H5T_NATIVE_DOUBLE, d);
  int modelOrder = (int) d[0];

  H5LTread_dataset(file, "/pitch", H5T_NATIVE_DOUBLE, d);
  double pitch = d[0];

  double * x = vector(N);
  H5LTread_dataset(file, "/data", H5T_NATIVE_DOUBLE, x);


  double pitch_bounds[2];
  H5LTread_dataset(file, "/pitchBounds", H5T_NATIVE_DOUBLE, pitch_bounds);
  #ifdef DOUBLE
  double * pb = pitch_bounds;
  #else
  float pb[2];
  pb[0] = (float)pitch_bounds[0]; pb[1] = (float)pitch_bounds[1];
  #endif

  status = H5Fclose(file);

  /* construct object and compute  */
  single_pitch * sp = new single_pitch(L, F, N, pb);

  #ifdef SINGLE
  float * xs = vector(N);
  for(int l = 0 ; l < N ; ++l)
    xs[l] = (float)x[l];
  sp->compute(xs);
  sp->refine(xs, 1e-8);
  del_vector(xs);
  #else
  sp->compute(x);
  sp->refine(x, 1e-8);
  #endif
  
  FTYPE omega_0h = sp->est_fast(x, 0.0, 1e-3);
  int estModelOrder = sp->modelOrder();
  ASSERT_NEAR(omega_0h/(2*M_PI), pitch, 1e-3); 
  ASSERT_EQ(modelOrder, estModelOrder);

  delete sp;
  del_vector(x);
}

int main(int argc, char **argv){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


