/*
 * ryla.c -- translation of SUBROUTINE RYLA (src/ryla.f90)
 *
 * Select method of evaluation
 */
#include "whamp.h"

void RYLA(cd Y, double AL, cd *RC)
{
    double AY;

    /*           ****  CHOOSE METHOD OF EVALUATION **** */
    if (AL < 4.) {
        RTAY(Y, AL, RC); /* Taylor series */
        return;
    }
    AY = cabs(Y);
    if (AY * AY > 75. * AL) {
        RTAY(Y, AL, RC); /* Taylor series */
        return;
    }
    if (AY > 40. + AL / 3) {
        RTAY(Y, AL, RC); /* Taylor series */
        return;
    }

    if (3. * (AL - 10.) > AY && AY * AY < 15. * AL) {
        RASY(Y, AL, RC); /* asymptotic series */
        return;
    } else {
        RINT(Y, AL, RC); /* numerical integration */
        return;
    }
}
