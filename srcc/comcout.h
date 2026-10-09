/* comcout.h -- translation of module comcout (src/comcout.f90) */
#ifndef COMCOUT_H
#define COMCOUT_H

#include "whamp.h"

extern cd X;        /* frequency */
extern cd EFL[4];   /* E field, indices 1..3 */
extern cd BFL[4];   /* B field, indices 1..3 */
extern cd D;        /* dispersion function */
extern cd DX;       /* D derivative wrt frequency */
extern cd DZ;       /* D derivative wrt Z */
extern cd DP;       /* D derivative wrt P */
extern cd E[7][5];  /* epsilon and its derivatives x eps_x, z eps_z, p eps_p; E(I,K) I=1..6 K=1..4 */
extern cd RI;       /* complex refractive index */

extern double P;
extern double Z;
extern double VG[3];
extern double SG[3];
extern double ENE;

#endif /* COMCOUT_H */
