# Optional - but remove link option -DMKL
MKLROOT = /opt/intel

#Dynamic linking
MKL = -m64 -I$(MKLROOT)/mkl/include -I$(MKLROOT)/mkl/include/fftw -L$(MKLROOT)/mkl/lib/intel64 -L$(MKLROOT)/compilers_and_libraries/linux/lib/intel64_lin -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -DALIGNMENT32

#FTYPE = -DDOUBLE -DREGULARIZE
FTYPE = -DDOUBLE 
#LINK =  -DMKL $(MKL) -DALIGNMENT32 -DCHEBYSHEV
LINK = -lfftw3 -DCHEBYSHEV


all:
	g++ -O3 -march=native -fPIC --shared -Wl,-soname,single_pitch.so -Wno-write-strings -DLIB $(FTYPE) ../src/single_pitch.cpp ../src/tools.cpp ../src/th.cpp ../src/vector.cpp $(LINK) -pthread  -o single_pitch.so
