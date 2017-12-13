#include "../src/th.hpp"
#include "gtest/gtest.h"
#include <stdlib.h>
#include <math.h>


#ifdef DOUBLE
#define FTYPE double 
#else
#define FTYPE float
#endif


class RealTest : public ::testing::Test{
protected:

  RealTest(){
  }

  virtual ~RealTest(){
  }

  virtual void SetUp(){
  }

  virtual void TearDown(){
  }
};

TEST_F(RealTest, th) {


  int N = 3;

  FTYPE t[] = { 0.2548401 , -1.19102012,  0.93416293, -0.37818981};
  FTYPE h[] = { -0.18860475, -0.33219973,  0.35968834, 
                -2.22949682, -0.3830574, 1.13525328,  0.37186914};

  FTYPE *gammah;
  int K = ((N+2)*(N+1))/2;
  gammah = vector(K);

  FTYPE gamma[] = { 15.0976782, -0.66822683, -0.02905703, 0.36881507, 
                    -0.14570701, -0.19041787, -0.38019016, -0.15277086, 
                    0.16812068,  0.16299686};

  th(N, t, h, gammah);


  for(int i = 0 ; i < K ; ++i){
    ASSERT_NEAR(gamma[i], gammah[i], 1e-8) << "Vectors differ at index " << i;
  }

  del_vector(gammah);
}


TEST_F(RealTest, thp) {


  int N = 3;

  FTYPE t[] = { 0.2548401 , -1.19102012,  0.93416293, -0.37818981};
  FTYPE h[] = { 0.93416293, -0.37818981,  0.35968834, -2.22949682, 
                -0.3830574, 1.13525328,  0.37186914};

  FTYPE *gammah;
  int K = ((N+2)*(N+1))/2;
  gammah = vector(K);

  FTYPE gamma[] = { -1.47205416, -1.37892347,  1.15243512, 1.35466476, 
                    -0.41218666,  1.01869934, 0.52769438,  0.22693268, 
                    0.10547072,  0.26054153};

  for(int k = 0; k < 2*N+1 ; ++k)
    h[k] = -h[k];

  thp(N, t, h, gammah);


  for(int i = 0 ; i < K ; ++i){
    ASSERT_NEAR(gamma[i], gammah[i], 1e-8) << "Vectors differ at index " << i;
  }

  del_vector(gammah);
}

TEST_F(RealTest, solve) {


  int N = 3;

  FTYPE t[] = { 0.2548401 , -1.19102012,  0.93416293, -0.37818981};
  FTYPE h[] = { -0.18860475, -0.33219973,  0.35968834, 
                -2.22949682, -0.3830574, 1.13525328,  0.37186914};
  FTYPE f[] = {-0.13876892,  0.1205414, -1.00767956,  0.05949563};
  FTYPE x[] = { 0.0107896, 0.3009554, -0.00617545, -0.12537081};
  FTYPE * xh = vector(N+1);
    
  solve(N, t, h, f, NULL, xh);

  for(int i = 0 ; i < N + 1 ; ++i){
    ASSERT_NEAR(x[i], xh[i], 1e-8) << "Vectors differ at index " << i;
  }

  del_vector(xh);
}


int main(int argc, char **argv){
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
