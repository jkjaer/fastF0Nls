#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdio.h>
#include <stdlib.h>

#include "tools.hpp"

#ifdef ALIGNMENT32
#define ALIGNMENT 32
#endif

#ifdef ALIGNMENT64
#define ALIGNMENT 64
#endif

#ifdef ALIGNMENT128
#define ALIGNMENT 128
#endif

FTYPE * vector(int N);

void del_vector(FTYPE * x);


#endif 
