
# Mandatory 
GTEST = /usr/lib/
HDF5 = /usr/lib/x86_64-linux-gnu/hdf5/serial

# Optional - but remove link option -DMKL
MKLROOT = /opt/intel

EXTLIB = -L$(HDF5)
INCLUDE = -I$(HDF5)/include
LIBSHDF = $(EXTLIB) -lhdf5_hl -lhdf5 -lz 

#Dynamic linking
MKL = -m64 -I$(MKLROOT)/mkl/include -I$(MKLROOT)/mkl/include/fftw -L$(MKLROOT)/mkl/lib/intel64 -L$(MKLROOT)/compilers_and_libraries/linux/lib/intel64_lin -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -DALIGNMENT32

FTYPE = -DDOUBLE

#LINK =  -DMKL $(MKL) -DCHEBYSHEV
#LINK = -DMKL $(MKL)
LINK = -lfftw3 -DCHEBYSHEV

all:
	g++ -g -fno-inline -Wno-write-strings $(FTYPE) $(INCLUDE) -isystem $(GTEST)/include ../src/single_pitch.cpp ../src/tools.cpp ../src/th.cpp ../src/vector.cpp test.cpp $(GTEST)/libgtest.a $(INCLUDE) $(LIBSHDF) $(LINK) -pthread -o test

	g++ -g -Wno-write-strings -DDOUBLE -isystem $(GTEST)/include ../src/th.cpp ../src/tools.cpp ../src/vector.cpp test_th.cpp $(GTEST)/libgtest.a $(LINK) -pthread -o test_th
