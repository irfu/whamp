/*
 * rasy.c -- translation of SUBROUTINE RASY (src/rasy.f90)
 *
 * Assymptotic series estimate
 */
#include "whamp.h"

void RASY(cd Y, double AL, cd *RC)
{
    int m, n;
    double A, AY, C;
    /* T is assigned by `T = C` on every iteration with m < AY; since
     * AY = |Y| + 2 >= 2 > 0 = m on the first iteration, T is always
     * written before the `C >= T` test can be reached.  Initialise to
     * silence -Wmaybe-uninitialized (the value is never observed). */
    double T = 0.;
    cd Y2, COT, P, PY, PP, PPY, PN, PYN, QN, QYN;
    static const double PI = 3.14159265358979;

/* RC is a 2x2 column-major array matching Fortran RC(2,2) */
#define RC11 RC[0]
#define RC21 RC[1]
#define RC12 RC[2]
#define RC22 RC[3]

    Y2 = Y * Y;
    COT = cos(PI * Y) / sin(PI * Y);
    /*                  1.E99 IS TOO BIG FOR S/370 HARDWARE. SET TO LARGEST */
    /*                  POSSIBLE FOR IBM MACHINES */
    /*     C=1.E99 */
    C = 7.2e35;
    PN = -Y / AL;
    PYN = PN;
    A = 1. / (AL * sqrt(2. * PI * AL));
    QN = PI * Y2 * COT * A;
    QYN = QN * (2. - Y * PI * COT) - Y * PI * PI * Y2 * A;

    P = PN + QN;
    PY = PYN + QYN;
    PP = -PN - 1.5 * QN;
    PPY = -PYN - 1.5 * QYN;
    AY = cabs(Y) + 2.;

    for (n = 1; n <= 100; n++) {
        m = n - 1;
        PYN = (PYN * ((double)(m * m) - Y2) - 2. * Y2 * PN) / ((double)(2 * m + 1) * AL);
        PN = PN * ((double)(m * m) - Y2) / ((double)(2 * m + 1) * AL);
        QYN = (QYN * ((m + .5) * (m + .5) - Y2) - 2. * Y2 * QN) / (2. * (double)n * AL);
        QN = QN * ((m + .5) * (m + .5) - Y2) / (2. * (double)n * AL);
        if ((double)m < AY) {
            /* fall through and accumulate */
        } else {
            C = (double)n * (cabs(PN) + cabs(QN));
            if (C <= F32(1.E-7) * cabs(PP))
                break;
            if (C >= T)
                break;
        }
        P = P + PN + QN;
        PY = PY + PYN + QYN;
        PP = PP - ((double)n + 1.) * PN - ((double)n + 1.5) * QN;
        PPY = PPY - ((double)n + 1.) * PYN - ((double)n + 1.5) * QYN;
        T = C;
    }

    RC11 = P + PN + QN;
    RC21 = PY + PYN + QYN;
    RC12 = PP + P;
    RC22 = PPY + PY;
}
