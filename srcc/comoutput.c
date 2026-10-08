/* comoutput.c -- storage for module comoutput (src/comoutput.f90) */
#include "comoutput.h"

int kperpSize, kparSize;

double *kperpOUT;
double *kparOUT;

cd *fOUT;
cd *ExOUT;
cd *EyOUT;
cd *EzOUT;
cd *BxOUT;
cd *ByOUT;
cd *BzOUT;
cd *SxOUT;
cd *SyOUT;
cd *SzOUT;

double *EBOUT;
double *VGPOUT;
double *VGZOUT;
double *SGPOUT;
double *SGZOUT;
double *uOUT;

int *flagSolutionFoundOUT;
int *flagTooHeavilyDampedOUT;
int *flagNoConvergenceOUT;
