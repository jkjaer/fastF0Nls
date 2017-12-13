#ifndef __TOOLS_H__
#define __TOOLS_H__

#include <cmath>

#ifdef DOUBLE
#define FTYPE double 

#define ABS fabs

#define FFTW_MALLOC fftw_malloc
#define FFTW_PLAN fftw_plan
#define FFTW_PLAN_DFT_R2C_1D fftw_plan_dft_r2c_1d
#define FFTW_PLAN_DFT_C2R_1D fftw_plan_dft_c2r_1d
#define FFTW_EXECUTE fftw_execute
#define FFTW_DESTROY_PLAN fftw_destroy_plan
#define FFTW_FREE fftw_free
#define FFTW_COMPLEX fftw_complex

#define COPY cblas_dcopy
#define SCAL cblas_dscal
#define AXPY cblas_daxpy
#define DOT cblas_ddot 
#define VADD vdAdd
#define VSUB vdSub
#define VSIN vdSin
#define VSINCOS vdSinCos
#define VUNPACKI vdUnpackI
#define VMUL vdMul
#define VDIV vdDiv
#define VINV vdInv

#define MKL_COMPLEX MKL_Complex16

#else
#define FTYPE float

#define ABS fabsf

#define FFTW_MALLOC fftwf_malloc
#define FFTW_PLAN fftwf_plan
#define FFTW_PLAN_DFT_R2C_1D fftwf_plan_dft_r2c_1d
#define FFTW_PLAN_DFT_C2R_1D fftwf_plan_dft_c2r_1d
#define FFTW_EXECUTE fftwf_execute
#define FFTW_DESTROY_PLAN fftwf_destroy_plan
#define FFTW_FREE fftwf_free
#define FFTW_COMPLEX fftwf_complex

#define COPY cblas_scopy
#define SCAL cblas_sscal
#define AXPY cblas_saxpy
#define DOT cblas_sdot 
#define VADD vsAdd
#define VSUB vsSub
#define VSIN vsSin
#define VSINCOS vsSinCos
#define VUNPACKI vsUnpackI
#define VMUL vsMul
#define VDIV vsDiv
#define VINV vsInv

#define MKL_COMPLEX MKL_Complex8
#endif

#define MAX(A,B) (A > B ? A : B)
#define MIN(A,B) (A < B ? A : B)

#ifdef MKL
#include <mkl.h>
#endif

void CvmdSinRange(int min, int max, FTYPE a, FTYPE * y);

void CvmdSinRange(int n, FTYPE * x, FTYPE a, FTYPE * t, FTYPE * y);

void CvmdExpRange(int min, int max, FTYPE a, FTYPE * y);

void CvmdExpRange(int n, FTYPE * x, FTYPE a, FTYPE *t1, FTYPE *t2, FTYPE * y);

void CvmdCosSinRange(int n, FTYPE * x, FTYPE a, FTYPE *t,
                     FTYPE * y, FTYPE * z);

void CvmdScalSquare(int n, FTYPE a, FTYPE * x, FTYPE * y);

void CvmdInit(int n, FTYPE a, FTYPE * x);

void CvmdDiv(int n, FTYPE * x, FTYPE * z, FTYPE * y);

void CvmdInverse(int n, FTYPE * x, FTYPE * y);

void CvmdMul(int n, FTYPE * x, FTYPE * z, FTYPE * y);

void CvmdMulz(int n, FTYPE * x, FTYPE * z, int s, FTYPE *t, FTYPE * v, FTYPE * w);

void CvmdMulz(int n, FTYPE * x, FTYPE * z, FTYPE *t, FTYPE * v, FTYPE * w);

void CvmdShift(int n, FTYPE *a, FTYPE *b, FTYPE *c);

void CvmdSymmetric(int n, FTYPE *a, FTYPE *b);

void CvmdReverse(int n, FTYPE *a, FTYPE *b);

FTYPE CvmdDot(int n, FTYPE *a, FTYPE *b);

void CvmdAddConstant(int n, FTYPE a, FTYPE *b, FTYPE *c);

void CvmdAdd(int n, FTYPE *a, FTYPE *b, FTYPE *c);

void CvmdSub(int n, FTYPE *a, FTYPE *b, FTYPE *c);

void CvmdAddInplace(int n, FTYPE *a, FTYPE *b);

void CvmdSubInplace(int n, FTYPE *a, FTYPE *b);

void CvmdAxpy(int n, FTYPE alpha, FTYPE *x, FTYPE *y);

void CvmdScal(int n, FTYPE alpha, FTYPE *x);

void CvmdCopy(int n, FTYPE *x, FTYPE *y);

void CvmdPack(int n, FTYPE * x, int s, FTYPE * y);

int CvmdArgmax(int n, FTYPE * x);

#endif
