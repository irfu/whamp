/*
 * rtay.c -- translation of SUBROUTINE RTAY (src/rtay.f90)
 *
 *                  ******** TAYLOR SERIES ********
 */
#include "whamp.h"

void RTAY(cd Y, double AL, cd *RC)
{
    int i;
    double T;
    cd Y2, PN, PYN, COT;

/* RC is a 2x2 column-major array matching Fortran RC(2,2) */
#define RC11 RC[0]
#define RC21 RC[1]
#define RC12 RC[2]
#define RC22 RC[3]

    Y2 = Y * Y;
    PN = Y / (Y2 - 1.);
    PYN = -Y * (Y2 + 1.) / ((Y2 - 1.) * (Y2 - 1.));
    RC11 = PN;
    RC12 = PN;
    RC21 = PYN;
    RC22 = PYN;

    for (i = 2; i <= 100; i++) {
        COT = (double)(2 * i - 1) / (Y2 - (double)(i * i)) * AL;
        PYN = COT * (PYN - 2. * Y2 / (Y2 - (double)(i * i)) * PN);
        PN = COT * PN;
        RC11 = RC11 + PN;
        RC21 = RC21 + PYN;
        RC12 = RC12 + i * PN;
        RC22 = RC22 + i * PYN;
        T = cabs(PN) * 1.E8;
        if (T < cabs(RC11))
            break;
    }
}
