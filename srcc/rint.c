/*
 * rint.c -- translation of SUBROUTINE RINT (src/rint.f90)
 *
 *                  ******** NUMERICAL INTEGRATION ********
 */
#include "whamp.h"

void RINT(cd YY, double AL, cd *RC)
{
    int i;
    double B, C, G, SIG, T, UL, Z, Xloc, XO, XY, YR;
    cd YA;
    cd Y, COT, Dv, EXF, F, H, O, Pv, R, RY, RP, RPY, S;
    /* ABSCISSAS FOR GAUSSIAN INTEGRATION */
    static const double A[16] = {-.989400934991649,
                                 -.944575023073232,
                                 -.865631202387831,
                                 -.755404408355003,
                                 -.617876244402643,
                                 -.458016777657227,
                                 -.281603550779258,
                                 -.095012509837637,
                                 .989400934991649,
                                 .944575023073232,
                                 .865631202387831,
                                 .755404408355003,
                                 .617876244402643,
                                 .458016777657227,
                                 .281603550779258,
                                 .095012509837637};
    static const double W[16] = {.027152459411754,
                                 .062253523938647,
                                 .095158511682492,
                                 .124628971255533,
                                 .149595988816576,
                                 .169156519395002,
                                 .182603415044923,
                                 .189450610455068,
                                 .027152459411754,
                                 .062253523938647,
                                 .095158511682492,
                                 .124628971255533,
                                 .149595988816576,
                                 .169156519395002,
                                 .182603415044923,
                                 .189450610455068};
    static const double PI = 3.14159265358979;

/* RC is passed as a 2x2 column-major array matching Fortran RC(2,2) */
#define RC11 RC[0]
#define RC21 RC[1]
#define RC12 RC[2]
#define RC22 RC[3]

    RC11 = 0.;
    RC21 = 0.;
    RC12 = 0.;
    RC22 = 0.;

    if (creal(YY) < 0.) {
        Y = -YY;
        SIG = -1.;
    } else {
        Y = YY;
        SIG = 1.;
    }
    YA = cimag(Y);
    YR = creal(Y);
    /* Fortran: UL = PI - REAL(2.8d0*Y/(36.d0 + Y)) : real part of the complex expr */
    UL = PI - creal(2.8 * Y / (36. + Y));
    COT = cos(PI * Y) / sin(PI * Y);
    Dv = PI * (1. + COT * COT);
    C = YR / AL;
    XO = log(C + sqrt(1. + C * C));

    for (i = 1; i <= 16; i++) {
        Xloc = UL / 2. * (1. + A[i - 1]);
        Z = sin(Xloc);
        C = cos(Xloc);
        G = YR / AL * Xloc / Z;
        T = sqrt(1. + G * G);
        B = log(G + T);
        G = (1. / Xloc - C / Z) * G / T;
        T = AL * (T * C - 1.);
        Z = creal(cexp(Xloc * YA));
        C = .5 * (Z + 1. / Z);
        S = .5 * I * (Z - 1. / Z);
        F = COT + G;
        H = 1. - G * COT;
        EXF = cexp(T - Y * B);
        O = B * C + Xloc * S;
        Pv = Xloc * C - B * S;
        XY = Xloc * YR;
        R = (F * C + H * S) * EXF;
        RY = (F * O - H * Pv + Dv * (C - G * S)) * EXF;
        RP = ((F * T - H * XY) * C + (H * T + F * XY) * S) * EXF;
        RPY = (F * (T * O - XY * Pv) - H * (T * Pv + XY * O) +
               ((T + XY * G) * C - (G * T - XY) * S) * Dv) *
              EXF;

        Xloc = XO / 2. * (1. + A[i - 1]);
        Z = exp(Xloc);
        C = (Z + 1. / Z) / 2. - 1.;
        Pv = exp(AL * C - Y * Xloc);
        RC11 = RC11 + W[i - 1] * (UL * R + XO * Pv);
        RC21 = RC21 - W[i - 1] * (UL * RY + XO * Xloc * Pv);
        RC12 = RC12 + W[i - 1] * (UL * RP + XO * AL * C * Pv);
        RC22 = RC22 - W[i - 1] * (UL * RPY + XO * AL * Xloc * C * Pv);
    }

    O = Y / AL;
    Pv = Y * Y / 2.;
    RC11 = O * (Y * RC11 / 2. - 1.) * SIG;
    RC21 = 2. * RC11 + O * (Pv * RC21 + 1.) * SIG;
    RC12 = Y * O * RC12 / 2. * SIG;
    RC22 = 2. * RC12 + O * Pv * RC22 * SIG;
}
