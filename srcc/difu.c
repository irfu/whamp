/*
 * difu.c -- translation of SUBROUTINE DIFU (src/difu.f90)
 *
 *        ARGUMENTS: KOL     =2 D AND ITS X-DERIVATIVE COMPUTED.
 *                           =4 ALL DERIVATIVES AND WAVE FIELDS.
 *                   JMAX    =1 - 6, NUMBER OF COMPONENTS.
 *                   IERR    ERROR FLAG
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"

void DIFU(int KOL, int JMAX, int *IERR)
{
    int IB, J, KOLLA;
    double Q;
    cd XSI[24];
    cd DF, U1, U2, U3, U12, U13, U32;
    cd Av, Bv, Cv, DA, DB, DC;

#define XS(I, K) XSI[(size_t)((I) - 1) + (size_t)((K) - 1) * 6]

    /*              *********** FORM DIELECTRIC TENSOR ************ */

    for (int i = 1; i <= 6; i++)
        for (int k = 1; k <= 4; k++)
            E[i][k] = 0.;
    E[1][1] = 1.;
    E[4][1] = 1.;
    E[6][1] = 1.;

    JMA = JMAX;
    KOLLA = KOL < 4 ? KOL : 4;
    for (J = 1; J <= JMA; J++) {
        E[1][1] = E[1][1] - XP[J];
        E[4][1] = E[4][1] - XP[J];
        E[6][1] = E[6][1] - XP[J];
        if (AA[J][1] == AA[J][2]) {
            AA[J][2] = 0.;
            DD[J] = 1.;
        }
        IB = 1;
        DF = XP[J] / (AA[J][1] * (AA[J][1] - AA[J][2]));
        Q = AA[J][1] - DD[J] * AA[J][2];
        for (;;) {
            CHI(XSI, J, IB, KOL, IERR);
            if (*IERR != 0)
                return;
            for (int k = 1; k <= KOLLA; k++)
                for (int i = 1; i <= 6; i++)
                    E[i][k] = E[i][k] + DF * Q * XS(i, k);

            if (IB == 2)
                break; /* cycle species_loop */
            Q = (DD[J] - 1.) * AA[J][1];
            if (Q == 0.)
                break; /* cycle species_loop */
            IB = 2;
        }
    }
    /*                       *** DIELECTRIC TENSOR COMPUTED *** */

    /*       ******* FORM REFRACTIVE INDEX, CV=SPEED OF LIGHT/THERM. SPEED. ** */
    U1 = PP[1] * CV / XX[1];
    U3 = ZZ[1] * CV / XX[1];
    U12 = U1 * U1;
    U32 = U3 * U3;
    U2 = U12 + U32;
    U13 = 2. * U1 * U3;

    /*                      ******** FORM DISPERSION FUNCTION ******** */
    Av = U12 * E[1][1] + U13 * E[3][1] + U32 * E[6][1];
    /*************  Sign error corrected in Feb. 1989.   Kjell R.  *********/
    /* ***  B=U2*(E(1,1)*E(6,1)-E(3,1)**2)+(U3*E(5,1)+U1*E(2,1))**2 */
    Bv = U2 * (E[1][1] * E[6][1] - E[3][1] * E[3][1]) +
         (U3 * E[5][1] - U1 * E[2][1]) * (U3 * E[5][1] - U1 * E[2][1]);
    Cv = (E[1][1] * E[6][1] - E[3][1] * E[3][1]) * E[4][1] + E[6][1] * E[2][1] * E[2][1];
    Cv = Cv + (E[1][1] * E[5][1] + 2. * E[2][1] * E[3][1]) * E[5][1];

    D = (U2 - E[4][1]) * Av - Bv + Cv;
    if (KOL <= 1)
        return;
    /*        ****** COMPLETE X-DERIVATIVE OF DIELECTRIC TENSOR ****** */
    E[1][2] = E[1][2] - 2. * (E[1][1] - 1.);
    E[2][2] = E[2][2] - 2. * E[2][1];
    E[3][2] = E[3][2] - 2. * E[3][1];
    E[4][2] = E[4][2] - 2. * (E[4][1] - 1.);
    E[5][2] = E[5][2] - 2. * E[5][1];
    E[6][2] = E[6][2] - 2. * (E[6][1] - 1.);
    /*          ****** X-DERIVATIVE OF DISPERSION FUNCTION ******* */
    DA = (E[1][2] - 2. * E[1][1]) * U12 + U13 * (E[3][2] - 2. * E[3][1]) +
         (E[6][2] - 2. * E[6][1]) * U32;
    DB = 2. * (U3 * E[5][1] - U1 * E[2][1]) *
             (U3 * (E[5][2] - E[5][1]) - U1 * (E[2][2] - E[2][1])) +
         U2 * ((E[1][2] - E[1][1]) * E[6][1] - 2. * (E[3][2] - E[3][1]) * E[3][1] +
               (E[6][2] - E[6][1]) * E[1][1]);
    DC = (E[1][2] * E[6][1] + E[1][1] * E[6][2] - 2. * E[3][1] * E[3][2]) * E[4][1];
    DC = DC + (E[1][1] * E[6][1] - E[3][1] * E[3][1]) * E[4][2];
    DC = DC + E[5][1] * (E[1][2] * E[5][1] + 2. * E[1][1] * E[5][2] +
                         E[2][1] * E[3][2] + 2. * E[2][2] * E[3][1]);
    DC = DC + E[2][1] * (E[6][2] * E[2][1] + 2. * E[6][1] * E[2][2] +
                         E[5][1] * E[3][2] + 2. * E[5][2] * E[3][1]);

    DX = ((U2 - E[4][1]) * DA - (2. * U2 + E[4][2]) * Av - DB + DC) / XX[1];
    if (KOL <= 2)
        return;
    DZ = 0.;
    if (ZZ[1] != 0.) {
        /*         ****** Z-DERIVATIVE OF DISPERSION FUNCTION ****** */
        DA = U12 * E[1][3] + U13 * (E[3][3] + E[3][1]) + U32 * (E[6][3] + 2. * E[6][1]);
        /*      DB=2.*(U3*E(5,1)-U1*E(3,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,1))+ */
        /* CCCCCCorrection 1991- 03 - 11, Kjell R. */
        /*      DB=2.*(U3*E(5,1)-U1*E(2,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,1))+ */
        /* CCCCCCorrection of error found by Scott Boardsen, NASA/MSFC. 1992-03-12, Kjell R. */
        DB = 2. * (U3 * E[5][1] - U1 * E[2][1]) * (U3 * (E[5][3] + E[5][1]) - U1 * E[2][3]) +
             2. * U32 * (E[1][1] * E[6][1] - E[3][1] * E[3][1]) +
             U2 * (E[1][3] * E[6][1] + E[1][1] * E[6][3] - 2. * E[3][1] * E[3][3]);
        DC = (E[1][3] * E[6][1] + E[1][1] * E[6][3] - 2. * E[3][1] * E[3][3]) * E[4][1];
        DC = DC + (E[1][1] * E[6][1] - E[3][1] * E[3][1]) * E[4][3];
        DC = DC + E[5][1] * (E[1][3] * E[5][1] + 2. * E[1][1] * E[5][3] +
                             E[2][1] * E[3][3] + 2. * E[2][3] * E[3][1]);
        DC = DC + E[2][1] * (E[6][3] * E[2][1] + 2. * E[6][1] * E[2][3] +
                             E[5][1] * E[3][3] + 2. * E[5][3] * E[3][1]);

        DZ = ((U2 - E[4][1]) * DA + (2. * U32 - E[4][3]) * Av - DB + DC) / ZZ[1];
    }
    if (KOL <= 3)
        return;
    /*        ****** P-DERIVATIVE OF DISPERSION FUNCTION ****** */
    DP = 0.;
    if (PP[1] != 0.) {
        DA = U12 * (E[1][4] + 2. * E[1][1]) + U13 * (E[3][4] + E[3][1]) + U32 * E[6][4];
        DB = 2. * (U3 * E[5][1] - U1 * E[2][1]) * (U3 * E[5][4] - U1 * (E[2][4] + E[2][1])) +
             2. * U12 * (E[1][1] * E[6][1] - E[3][1] * E[3][1]) +
             U2 * (E[1][4] * E[6][1] + E[1][1] * E[6][4] - 2. * E[3][1] * E[3][4]);
        DC = (E[1][4] * E[6][1] + E[1][1] * E[6][4] - 2. * E[3][1] * E[3][4]) * E[4][1];
        DC = DC + (E[1][1] * E[6][1] - E[3][1] * E[3][1]) * E[4][4];
        DC = DC + E[5][1] * (E[1][4] * E[5][1] + 2. * E[1][1] * E[5][4] +
                             E[2][1] * E[3][4] + 2. * E[2][4] * E[3][1]);
        DC = DC + E[2][1] * (E[6][4] * E[2][1] + 2. * E[6][1] * E[2][4] +
                             E[5][1] * E[3][4] + 2. * E[5][4] * E[3][1]);

        DP = ((U2 - E[4][1]) * DA + (2. * U12 - E[4][4]) * Av - DB + DC) / PP[1];
    }
    /*                     ******** COMPUTE ELECTRIC FIELD ******** */
    ENERGY(U1, U3, U2, U12, U32);
    /*       THE ELECTRIC FIELD IS 1 MV/M. */
    /*       THE MAGNETIC FIELD WILL BE IN GAMMA */
}
