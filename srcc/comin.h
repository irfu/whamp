/* comin.h -- translation of module comin (src/comin.f90) */
#ifndef COMIN_H
#define COMIN_H

#include "whamp.h"

extern double ASS[11];   /* particle mass in unit of mp, if 0 then use e- mass */
extern double AA[11][3]; /* alpha parameters in distribution function, AA(I,K) K=1,2 */
extern double DD[11];    /* Delta_j in distribution function */
extern double VD[11];    /* drift velocity normalized against parallel thermal velocity */
extern double DN[11];    /* number density of each plasma component cm^-3 */
extern double TA[11];    /* temperature in keV */
extern double PM[4];     /* array (pmin,pmax,pstep), indices 1..3 */
extern double ZM[4];     /* array (zmin,zmax,zstep), indices 1..3 */
extern double XOI;       /* initialization frequency (real part) */
extern double XC;        /* XC - gyrofrequency (kHZ) */
extern double PZL;       /* option L, =0 linear input, =1 logarithmic input of P and Z */

extern int printDebugInfo;

extern int cycleZFirst; /* which direction vary first (KFS in WHAMP) */

/* ------ below are variables calculated by WHAMP or WHAMP_ENGINE ------- */

extern int JMA;        /* the number of last plasma species with non-zero density */
extern double PP[11];  /* ?? */
extern double ZZ[11];  /* ?? */
extern double CV;      /* speed of light / thermal velocity of 1st species */
extern cd XX[11];
extern cd XP[11];      /* wpj^2/w^2 */
extern double PX;      /* plasma frequency */
extern double DEN;     /* total electron density */

/* ------ root finding parameters */

extern int maxIterations;

#endif /* COMIN_H */
