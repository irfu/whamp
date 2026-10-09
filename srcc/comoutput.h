/* comoutput.h -- translation of module comoutput (src/comoutput.f90)
 *
 * The Fortran allocatable arrays are represented as malloc'ed buffers
 * with 1-based indexing: a Fortran array A(kperpSize, kparSize) becomes
 * a C array A[1..kperpSize][1..kparSize] stored column-major (Fortran
 * order) in a flat buffer indexed by AOUT(i, j).
 */
#ifndef COMOUTPUT_H
#define COMOUTPUT_H

#include "whamp.h"

extern int kperpSize, kparSize;

extern double *kperpOUT; /* kperpOUT[1..kperpSize] */
extern double *kparOUT;  /* kparOUT[1..kparSize] */

extern cd *fOUT;
extern cd *ExOUT;
extern cd *EyOUT;
extern cd *EzOUT;
extern cd *BxOUT;
extern cd *ByOUT;
extern cd *BzOUT;
extern cd *SxOUT;
extern cd *SyOUT;
extern cd *SzOUT;

extern double *EBOUT;
extern double *VGPOUT;
extern double *VGZOUT;
extern double *SGPOUT;
extern double *SGZOUT;
extern double *uOUT;

extern int *flagSolutionFoundOUT;
extern int *flagTooHeavilyDampedOUT;
extern int *flagNoConvergenceOUT;

/* 1-based element accessors for the 2D output arrays (Fortran order) */
#define F2D(arr, i, j) ((arr)[(size_t)((i) - 1) + (size_t)((j) - 1) * (size_t)kperpSize])

/* Allocate 1-based array of n elements (element 0 unused) */
#define ALLOC1(arr, n, type)                                                     \
    do {                                                                         \
        if (arr) free(arr);                                                      \
        arr = (type *)calloc((size_t)(n) + 1, sizeof(type));                     \
    } while (0)

/* Allocate 1-based 2D array (kperpSize x kparSize) in Fortran order */
#define ALLOC2(arr, type)                                                        \
    do {                                                                         \
        if (arr) free(arr);                                                      \
        arr = (type *)calloc((size_t)kperpSize * (size_t)kparSize,               \
                             sizeof(type));                                      \
    } while (0)

#endif /* COMOUTPUT_H */
