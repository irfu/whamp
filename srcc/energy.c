/*
 * energy.c -- translation of SUBROUTINE ENERGY (src/energy.f90)
 *
 * The polarization vectors are here calculated from the
 * matrix of cofactors (e.g., D. B. Melrose, Plasma Astrophysics,
 * Vol I, p. 43, Gordon and Breach, 1980).
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"

void ENERGY(cd U1, cd U3, cd U2, cd U12, cd U32)
{
    int i, j;
    double Q, V;
    cd A, B, C;
    cd ERG, CM[4][4], DDER[4][4], ECON, EE[7][5];

    for (i = 1; i <= 6; i++)
        for (j = 1; j <= 4; j++)
            EE[i][j] = E[i][j];

    DDER[1][1] = EE[1][1] + EE[1][2] + U32;
    DDER[1][2] = EE[2][1] + EE[2][2];
    DDER[1][3] = EE[3][1] + EE[3][2] - U1 * U3;
    DDER[2][2] = EE[4][1] + EE[4][2] + U2;
    DDER[2][3] = EE[5][1] + EE[5][2];
    DDER[3][3] = EE[6][1] + EE[6][2] + U12;
    DDER[2][1] = -DDER[1][2];
    DDER[3][1] = DDER[1][3];
    DDER[3][2] = -DDER[2][3];

    EE[1][1] = EE[1][1] - U32;
    EE[3][1] = EE[3][1] + U1 * U3;
    EE[4][1] = EE[4][1] - U2;
    EE[6][1] = EE[6][1] - U12;
    /*  Form the elements of the cofactor matrix */
    CM[1][1] = EE[4][1] * EE[6][1] + EE[5][1] * EE[5][1];
    /***** Correction 1994-04-13 so that Im E_y < 0 now corresponds
     *     to L-mode waves again.   Kjell R.
     *****  CM(2,1)=EE(2,1)*EE(6,1)
     */
    CM[1][2] = EE[2][1] * EE[6][1] + EE[5][1] * EE[3][1];
    CM[1][3] = EE[2][1] * EE[5][1] - EE[4][1] * EE[3][1];
    CM[2][2] = EE[1][1] * EE[6][1] - EE[3][1] * EE[3][1];
    /***** Correction 1994-04-13 so that Im E_y < 0 now corresponds
     *     to L-mode waves again.   Kjell R.
     *  ***  CM(3,2)=EE(5,1)*EE(1,1)
     */
    CM[2][3] = EE[5][1] * EE[1][1] + EE[3][1] * EE[2][1];
    CM[3][3] = EE[1][1] * EE[4][1] + EE[2][1] * EE[2][1];
    /***** Correction 1994-04-13 so that Im E_y < 0 now corresponds
     *     to L-mode waves again.   Kjell R.
     *  ***   CM(1,2)=-CM(2,1)
     */
    CM[2][1] = -CM[1][2];
    CM[3][1] = CM[1][3];
    /*  ***   CM(2,3)=-CM(3,2) */
    CM[3][2] = -CM[2][3];

    V = 1. / 299.792458;
    for (i = 1; i <= 3; i++) {
        EFL[i] = 0.;
        for (j = 1; j <= 3; j++) {
            EFL[i] = EFL[i] + CM[j][i];
        }
    }

    if (cabs(EFL[1]) < 1000. * cabs(EFL[3])) {
        ECON = conjg(EFL[3]);
    } else {
        ECON = conjg(EFL[1]);
    }

    A = EFL[1] * ECON;
    B = EFL[2] * ECON;
    C = EFL[3] * ECON;

    Q = creal(A * conjg(A) + B * conjg(B) + C * conjg(C));
    Q = sqrt(Q);

    EFL[1] = A / Q;
    EFL[2] = B / Q;
    EFL[3] = C / Q;
    BFL[1] = -V * U3 * EFL[2];
    BFL[2] = V * (U3 * EFL[1] - U1 * EFL[3]);
    BFL[3] = V * U1 * EFL[2];

    ERG = 0.;
    for (i = 1; i <= 3; i++) {
        for (j = 1; j <= 3; j++) {
            ERG = ERG + conjg(EFL[i]) * DDER[i][j] * EFL[j];
        }
    }

    ENE = creal(ERG * conjg(ERG) / (ERG + conjg(ERG)));
}
