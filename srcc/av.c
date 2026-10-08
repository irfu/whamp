/*
 * av.c -- translation of SUBROUTINE AV (src/av.f90)
 *
 * This subroutine calculates the individual contribution
 * of plasma species to the energy density and the energy
 * flux.
 * Energy density units - electric field energy density is 1
 * Energy flux units - part of the total energy flux
 *
 * Loop through all species
 * For each species
 *
 * Variables
 *           EE(6,4) - dielectric tensor and its derivatives
 *             ENDEN - energy density in particular species
 *           ENFLUXZ - energy flux z-direction        in part. species
 *           ENFLUXP - energy flux p-direction in part. species
 *               JMA - number of species
 *           ENFIELD - energy in fields
 *              POYN - Poynting flux
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"

void AV(void)
{
    int i, IB, IERR, J, k;
    double coef, ENFLUXZ, ENFLUXP, ENDEN, ENFIELD, Q;
    double A, B;
    cd XSI[24];
    cd DF, U1, U2, U3, U12, U13, U32;
    cd ERG, DDER[4][4], EE[7][5];
    cd POYN[4];
    cd ENFLZ, ENFLP, DEP[4][4], DEZ[4][4];

    static const double PI = 3.1415926535897;

#define XS(I, K) XSI[(size_t)((I) - 1) + (size_t)((K) - 1) * 6]

    /* *** estimate the number of species JMA */
    for (J = 1; J <= 10; J++) {
        if (DN[J] == 0.)
            break;
        JMA = J;
    }
    /* *** calculate Poynting flux POYN in whamp units */
    coef = 10.0 / 4.0 / PI / 2.0;
    coef = coef * 4 / F32(8.8542e-12) * CV / F32(2.9979e8);
    POYN[1] = creal(EFL[2] * conjg(BFL[3]) - EFL[3] * conjg(BFL[2])) * coef;
    POYN[2] = creal(EFL[3] * conjg(BFL[1]) - EFL[1] * conjg(BFL[3])) * coef;
    POYN[3] = creal(EFL[1] * conjg(BFL[2]) - EFL[2] * conjg(BFL[1])) * coef;
    (void)POYN; /* computed but unused, as in the original Fortran */

    /* *** calculate energy in fields */
    A = creal(EFL[1] * conjg(EFL[1]) + EFL[2] * conjg(EFL[2]) + EFL[3] * conjg(EFL[3]));
    B = creal(BFL[1] * conjg(BFL[1]) + BFL[2] * conjg(BFL[2]) + BFL[3] * conjg(BFL[3]));
    ENFIELD = 1.0 + B / A * (double)(float)299.79 * (double)(float)299.79;

    /* *** Loop through all species */
    for (J = 1; J <= JMA; J++) {
        /*              *********** FORM DIELECTRIC TENSOR ************ */

        for (int ii = 1; ii <= 6; ii++)
            for (int kk = 1; kk <= 4; kk++)
                E[ii][kk] = 0.;
        E[1][1] = 1.;
        E[4][1] = 1.;
        E[6][1] = 1.;

        /* TODO: double check that this does not include vacuum energy */
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
        IERR = 0;
        for (;;) {
            CHI(XSI, J, IB, 4, &IERR);
            if (IERR != 0)
                return;
            for (k = 1; k <= 4; k++) {
                for (i = 1; i <= 6; i++) {
                    E[i][k] = E[i][k] + DF * Q * XS(i, k);
                }
            }

            if (IB == 2)
                break;
            Q = (DD[J] - 1.) * AA[J][1];
            if (Q == 0.)
                break;
            IB = 2;
        }
        /*                       *** DIELECTRIC TENSOR COMPUTED *** */

        /*       ******* FORM REFRACTIVE INDEX, CV=SPEED OF LIGHT/THERM. SPEED. ** */
        U1 = PP[1] * CV / XX[1];
        U3 = ZZ[1] * CV / XX[1];
        U12 = U1 * U1;
        U32 = U3 * U3;
        U2 = U12 + U32;
        U13 = 2. * U1 * U3;
        (void)U13; /* computed but unused, as in the original Fortran */

        /*        ****** COMPLETE X-DERIVATIVE OF DIELECTRIC TENSOR ****** */
        E[1][2] = E[1][2] - 2. * (E[1][1] - 1.);
        E[2][2] = E[2][2] - 2. * E[2][1];
        E[3][2] = E[3][2] - 2. * E[3][1];
        E[4][2] = E[4][2] - 2. * (E[4][1] - 1.);
        E[5][2] = E[5][2] - 2. * E[5][1];
        E[6][2] = E[6][2] - 2. * (E[6][1] - 1.);

        for (int ii = 1; ii <= 6; ii++)
            for (int kk = 1; kk <= 4; kk++)
                EE[ii][kk] = E[ii][kk];

        DDER[1][1] = EE[1][1] + EE[1][2] + U32;
        DDER[1][2] = EE[2][1] + EE[2][2];
        DDER[1][3] = EE[3][1] + EE[3][2] - U1 * U3;
        DDER[2][2] = EE[4][1] + EE[4][2] + U2;
        DDER[2][3] = EE[5][1] + EE[5][2];
        DDER[3][3] = EE[6][1] + EE[6][2] + U12;
        DDER[2][1] = -DDER[1][2];
        DDER[3][1] = DDER[1][3];
        DDER[3][2] = -DDER[2][3];

        ERG = 0.;
        for (i = 1; i <= 3; i++) {
            for (k = 1; k <= 3; k++) {
                ERG = ERG + conjg(EFL[i]) * DDER[i][k] * EFL[k];
            }
        }

        ENDEN = creal(ERG * conjg(ERG) / (ERG + conjg(ERG))) * 2.0 - ENFIELD;

        DEZ[1][1] = EE[1][3];
        DEZ[1][2] = EE[2][3];
        DEZ[1][3] = EE[3][3];
        DEZ[2][2] = EE[4][3];
        DEZ[2][3] = EE[5][3];
        DEZ[3][3] = EE[6][3];
        DEZ[2][1] = -DEZ[1][2];
        DEZ[3][1] = DEZ[1][3];
        DEZ[3][2] = -DEZ[2][3];

        DEP[1][1] = EE[1][4];
        DEP[1][2] = EE[2][4];
        DEP[1][3] = EE[3][4];
        DEP[2][2] = EE[4][4];
        DEP[2][3] = EE[5][4];
        DEP[3][3] = EE[6][4];
        DEP[2][1] = -DEP[1][2];
        DEP[3][1] = DEP[1][3];
        DEP[3][2] = -DEP[2][3];

        ENFLZ = 0.;
        ENFLP = 0.;
        for (i = 1; i <= 3; i++) {
            for (k = 1; k <= 3; k++) {
                ENFLZ = ENFLZ + conjg(EFL[i]) * DEZ[i][k] * EFL[k];
                ENFLP = ENFLP + conjg(EFL[i]) * DEP[i][k] * EFL[k];
            }
        }

        ENFLUXZ = creal(ENFLZ * conjg(ENFLZ) / (ENFLZ + conjg(ENFLZ))) * 2.0;
        ENFLUXP = creal(ENFLP * conjg(ENFLP) / (ENFLP + conjg(ENFLP))) * 2.0;

        A = F32(2.9979e8) / CV * F32(8.8542e-12) / 4.0;
        ENFLUXZ = -ENFLUXZ * creal(X) / Z * A;
        ENFLUXP = -ENFLUXP * creal(X) / P * A;
        /* Fortran FORMAT: ' enden= ', e9.3, ' enfl_p=', e9.3, ' enfl_z= ', e9.3, ' ', $ */
        printf(" enden= ");
        fE(9, 3, 0, ENDEN);
        printf(" enfl_p=");
        fE(9, 3, 0, ENFLUXP);
        printf(" enfl_z= ");
        fE(9, 3, 0, ENFLUXZ);
        printf(" ");
    }
}
