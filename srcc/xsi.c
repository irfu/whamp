/*
 * xsi.c -- translation of SUBROUTINE CHI (src/xsi.f90)
 *
 *        ARGUMENTS: XSI    CONTAINS THE SUSCEPTIBILITY TENSOR
 *                          ON RETURN.  (6x4, Fortran order, 1-based)
 *                   J      COMPONENT NUMBER.
 *                   IB     INDEX FOR AA (ALPHA)
 *                   KOL    DETERMINES WHETHER DERIVATIVES
 *                          SHOULD BE EVALUATED.
 *                   IERR   ERROR FLAG, IS SET =1 IF DAMPING
 *                          IS TOO STRONG.
 */
#include "whamp.h"
#include "comin.h"

void CHI(cd *XSI, int J, int IB, int KOL, int *IERR)
{
    int L;
    double A, ABZ, AL, ALA, Pv, RX, VZ, Zl, XI, XR;
    cd Xloc, XY, AI, BL, CL, DL, BLY, RC[4];
    cd PS, PSP, PSY, PPY, DPs, DPP, Y, ZY;

/* XSI is a 6x4 column-major array matching Fortran XSI(6,4) */
#define XS(I, K) XSI[(size_t)((I) - 1) + (size_t)((K) - 1) * 6]

    /*          **** RESIDUES FOR PADE APPROXIMANT **** */
    static const cd B[8] = {(-1.734012457471826e-2) + (-4.630639291680322e-2) * I,
                            (-1.734012457471826e-2) + (4.630639291680322e-2) * I,
                            (-7.399169923225014e-1) + (8.395179978099844e-1) * I,
                            (-7.399169923225014e-1) + (-8.395179978099844e-1) * I,
                            (5.840628642184073e0) + (9.536009057643667e-1) * I,
                            (5.840628642184073e0) + (-9.536009057643667e-1) * I,
                            (-5.583371525286853e0) + (-1.120854319126599e1) * I,
                            (-5.583371525286853e0) + (1.120854319126599e1) * I};
    static const cd C[8] = {(2.237687789201900e0) + (-1.625940856173727e0) * I,
                            (-2.237687789201900e0) + (-1.625940856173727e0) * I,
                            (1.465234126106004e0) + (-1.789620129162444e0) * I,
                            (-1.465234126106004e0) + (-1.789620129162444e0) * I,
                            (.8392539817232638e0) + (-1.891995045765206e0) * I,
                            (-.8392539817232638e0) + (-1.891995045765206e0) * I,
                            (.2739362226285564e0) + (-1.941786875844713e0) * I,
                            (-.2739362226285564e0) + (-1.941786875844713e0) * I};

    Xloc = XX[J];
    Zl = ZZ[J];
    Pv = PP[J];
    A = AA[J][IB];
    VZ = VD[J] * Zl;
    AI = I * A;
    if (ASS[J] == 0.)
        AI = -AI;
    AL = .5 * Pv * Pv;
    ALA = A * AL;

    for (int k = 0; k < 24; k++)
        XSI[k] = 0.;

    /*                         TEST FOR STRONG DAMPING */
    /*  Some tests modified to allow Z < 0. (1996-02-20, Kjell R.) */
    XI = cimag(Xloc);
    XR = fabs(creal(Xloc) - Zl * VD[J]);
    ABZ = fabs(Zl);
    for (;;) {
        if (XI >= -ABZ)
            break;
        /*  Allow high frequency, cold plasma waves ...(ADDED 1993-04-02, KJELL R.) */
        if (XR > 1. + 5. * ABZ && XR > 1. + 5. * Pv)
            break;
        RX = XR - fint(XR);
        if (RX > 0.5)
            RX = 1. - RX;
        if (XI >= -(double)(float)0.6 * RX)
            break;
        *IERR = 1;
        return;
    }

    XS(1, 1) = A;
    XS(6, 1) = A + 2 * A * VD[J] * VD[J];
    for (L = 1; L <= 8; L++) {
        BL = B[L - 1];
        /*  Modified to allow Z < 0. (1996-02-20, Kjell R.) */
        /* **** CL=C(L) */
        if (Zl < 0.) {
            CL = -C[L - 1];
        } else {
            CL = C[L - 1];
        }
        DL = CL + VD[J];
        Y = Xloc - DL * Zl;
        BLY = BL / Y;
        /*        ****** EVALUATE THE R-FUNCTION ****** */
        RYLA(Y, ALA, RC);
        XY = 1. + A * Zl * CL / Y;
        PS = XY * RC[0];
        PSP = XY * RC[2];
        /*                     ****** FORM SUCEPTIBILITY TENSOR ****** */
        XS(1, 1) = XS(1, 1) + A * BL * Y * PS;
        XS(2, 1) = XS(2, 1) + AI * BL * PSP;
        XS(3, 1) = XS(3, 1) + A * Pv * BL * DL * PS;
        XS(4, 1) = XS(4, 1) + BLY * PSP;
        /* CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R. */
        /* ****      XSI(5,1)=XSI(5,1)+AI*P*BLY*DL*PSP */
        XS(5, 1) = XS(5, 1) - AI * Pv * BLY * DL * PSP;
        XS(6, 1) = XS(6, 1) + 2. * A * BLY * DL * DL * (Xloc - VZ + AL * PS);

        if (KOL <= 1)
            continue;
        /*                    ****** FORM X-DERIVATIVES OF XSI ****** */
        PSY = XY * RC[1];
        PPY = XY * RC[3];
        XY = Xloc / Y;
        DPs = XY * (PSY - A * Zl * CL / Y * RC[0]);
        DPP = XY * (PPY - A * Zl * CL / Y * RC[2]);
        XS(1, 2) = XS(1, 2) + A * BL * (Y * DPs + Xloc * PS);
        XS(2, 2) = XS(2, 2) + AI * BL * DPP;
        XS(3, 2) = XS(3, 2) + A * Pv * BL * DL * DPs;
        XS(4, 2) = XS(4, 2) + BLY * (DPP - XY * PSP);
        /* CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R. */
        /* ****       XSI(5,2)=XSI(5,2)+AI*P*BLY*DL*(DPP-XY*PSP) */
        XS(5, 2) = XS(5, 2) - AI * Pv * BLY * DL * (DPP - XY * PSP);
        XS(6, 2) = XS(6, 2) + 2. * A * BLY * DL * DL *
                                    (Xloc + AL * DPs - XY * (Xloc - VZ + AL * PS));

        if (KOL <= 2)
            continue;
        /*                    ****** FORM Z-DERIVATIVES OF XSI ****** */
        ZY = Zl / Y;
        DPs = ZY * (A * CL * XY * RC[0] - DL * PSY);
        /* CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R. */
        /* ****       XSI(5,2)=XSI(5,2)+AI*P*BLY*DL*(DPP-XY*PSP) */
        XS(5, 2) = XS(5, 2) - AI * Pv * BLY * DL * (DPP - XY * PSP);

        DPP = ZY * (A * CL * XY * RC[2] - DL * PPY);
        ZY = DL * ZY;
        XS(1, 3) = XS(1, 3) + A * BL * Y * (DPs - ZY * PS);
        XS(2, 3) = XS(2, 3) + AI * BL * DPP;
        XS(3, 3) = XS(3, 3) + A * Pv * BL * DL * DPs;
        XS(4, 3) = XS(4, 3) + BLY * (DPP + ZY * PSP);
        /* CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R. */
        /* ****       XSI(5,3)=XSI(5,3)+AI*P*BLY*DL*(DPP+ZY*PSP) */
        XS(5, 3) = XS(5, 3) - AI * Pv * BLY * DL * (DPP + ZY * PSP);
        XS(6, 3) = XS(6, 3) + 2. * A * BLY * DL * DL *
                                    (AL * DPs - VZ + ZY * (Xloc - VZ + AL * PS));

        if (KOL <= 3)
            continue;
        /*                   ****** FORM P-DERIVATIVES OF XSI ****** */
        RYLA(Y - 1., ALA, RC);
        DPs = 2. * (PSP - PS);
        /*      DPP=2.*AL*((Y/(Y-1.))**2*RC(1,2)-PSP)-Y*DP */
        /* CCCCCCorrection, 1991 - 03 - 13, Kjell R */
        /*      DPP=2.*AL*((Y/(Y-1.))**2*(1.+A*Z*CL/Y)*RC(1,2)-PSP)-Y*DP */
        /*     As pointed out by Scott Boardsen, NASA/MSFC, the AL above should
        *     really be ALA.   Corrected  1991 - 12 - 18, Kjell R */
        {
            cd yratio = Y / (Y - 1.);
            DPP = 2. * ALA * (yratio * yratio * (1. + A * Zl * CL / Y) * RC[2] - PSP) -
                  Y * DPs;
        }
        XS(1, 4) = XS(1, 4) + A * BL * Y * DPs;
        XS(2, 4) = XS(2, 4) + AI * BL * DPP;
        XS(3, 4) = XS(3, 4) + A * Pv * BL * DL * (DPs + PS);
        XS(4, 4) = XS(4, 4) + BLY * (2. * PSP + DPP);
        /* CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R. */
        /* ****       XSI(5,4)=XSI(5,4)+AI*P*BLY*DL*(PSP+DPP) */
        XS(5, 4) = XS(5, 4) - AI * Pv * BLY * DL * (PSP + DPP);
        XS(6, 4) = XS(6, 4) + 4. * A * BLY * DL * DL * AL * PSP;
    }
    /*      **** COMPLETE XSI(4, ) **** */
    XS(4, 1) = XS(1, 1) - 2. * A * A * AL * XS(4, 1);
    XS(4, 2) = XS(1, 2) - 2. * A * A * AL * XS(4, 2);
    XS(4, 3) = XS(1, 3) - 2. * A * A * AL * XS(4, 3);
    XS(4, 4) = XS(1, 4) - 2. * A * A * AL * XS(4, 4);
}
