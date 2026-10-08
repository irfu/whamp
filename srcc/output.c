/*
 * output.c -- translation of SUBROUTINE OUTPT (src/output.f90)
 *
 * The Fortran ENTRY INOUT is translated as a separate function INOUT().
 * All local variables of OUTPT carry a Fortran SAVE statement and are
 * shared with the ENTRY, so they are file-scope statics here.
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"

/* Locals of OUTPT/INOUT with Fortran SAVE semantics */
static double andE, andphi, coef_poynt, dla, dlb, dle, Sx, Sy, Sz;
static cd andbz;
static int IOS, J, KMX, M, N, K;
static char IOU[21];
static char IC;

static const double PI = 3.14159265358979;

/* Helper: emulate Fortran  READ(*,'(A)')  -- pad with blanks, return 0 on
 * success, nonzero on end-of-file.  buf must hold len+1 bytes. */
int fortran_read_line(char *buf, int len)
{
    int c, LoopI = 0;
    while (LoopI < len) {
        c = getchar();
        if (c == EOF) {
            if (LoopI == 0) {
                buf[0] = '\0';
                return 1;
            }
            break;
        }
        if (c == '\n')
            break;
        buf[LoopI++] = (char)c;
    }
    while (LoopI < len)
        buf[LoopI++] = ' ';
    buf[len] = '\0';
    return 0;
}

void OUTPT(void)
{
    K = 0;
    for (;;) {
        K = K + 1;
        if (K > KMX)
            printf("  \n"); /* 6 FORMAT('  ') */
        if (K > KMX)
            return;
        IC = IOU[K - 1];
        if (IC == '/') {
            printf("  \n"); /* 6 FORMAT('  ') */
        } else {
            switch (IC) {
            case 'O':
            case 'o': /* O */
                printf(" elip= "); /* 772 FORMAT(' elip= ', e8.2, ' ', $) */
                fE(8, 2, 0, andE);
                printf(" ");
                /* ** bp/bz amplitude ** */
                break;
            case 'L':
            case 'l': /* L */
                printf(" bp/bz= "); /* 774 FORMAT(' bp/bz= ', e8.2, ' ', $) */
                fE(8, 2, 0, sqrt(cimag(BFL[1]) * cimag(BFL[1]) +
                                 creal(BFL[2]) * creal(BFL[2])) /
                                cabs(BFL[3]));
                printf(" ");
                break;
            case 'M':
            case 'm': /* M */
                /* ** bx/by amplitude ** */
                printf(" bx/by= "); /* 776 format(' bx/by= ', e8.2, ' ', $) */
                fE(8, 2, 0, fabs(cimag(BFL[1]) / creal(BFL[2])));
                printf(" ");
                break;
            case 'N':
            case 'n': /* N */
            {
                /* ** ellipticity min(bx,by)/max(bx.by) ** */
                double a1 = fabs(cimag(BFL[1]));
                double a2 = fabs(creal(BFL[2]));
                double sgn = (creal(BFL[1]) * cimag(BFL[2]) - cimag(BFL[1]) * creal(BFL[2])) /
                             fabs(creal(BFL[1]) * cimag(BFL[2]) -
                                  cimag(BFL[1]) * creal(BFL[2]));
                printf(" e= "); /* 778 FORMAT(' e= ', e8.2, ' ', $) */
                fE(8, 2, 0, (a1 < a2 ? a1 : a2) / (a1 > a2 ? a1 : a2) * sgn);
                printf(" ");
            } break;
            case 'X':
            case 'x': /* X */
            {
                /* ** bz-bx phase angle ** */
                double ar, mx;
                andbz = BFL[3] * conjg(BFL[1]) / cabs(BFL[1]);
                ar = creal(andbz);
                mx = fabs(ar) > cabs(andbz) ? fabs(ar) : cabs(andbz);
                andphi = acos(ar / mx) * 180 / PI;
                andphi = FSIGN(andphi, cimag(andbz));
                /* 792 FORMAT(' phi(bz-bx)= ', f8.2, ' ', $) */
                printf(" phi(bz-bx)= ");
                fF(8, 2, andphi);
                printf(" ");
            } break;
            case 'V':
            case 'v': /* V */
                /*     write Poynting vector uW/m^2 */
                coef_poynt = 10.0 / 4.0 / PI / 2.0;
                Sx = creal(EFL[2] * conjg(BFL[3]) - EFL[3] * conjg(BFL[2])) * coef_poynt;
                Sy = creal(EFL[3] * conjg(BFL[1]) - EFL[1] * conjg(BFL[3])) * coef_poynt;
                Sz = creal(EFL[1] * conjg(BFL[2]) - EFL[2] * conjg(BFL[1])) * coef_poynt;
                /* 794 FORMAT(' Sx= ', e8.2, ' Sy= ', e8.2, ' Sz= ', e9.3, ' ', $) */
                printf(" Sx= ");
                fE(8, 2, 0, Sx);
                printf(" Sy= ");
                fE(8, 2, 0, Sy);
                printf(" Sz= ");
                fE(9, 3, 0, Sz);
                printf(" ");
                break;
            case 'Y':
            case 'y': /* Y */
                /*     * abs(e)/abs(b) * */
                AV();
                break;
            case 'H':
            case 'h': /* H */
                dle = creal(EFL[1] * conjg(EFL[1]) + EFL[2] * conjg(EFL[2]) +
                            EFL[3] * conjg(EFL[3]));
                dlb = creal(BFL[1] * conjg(BFL[1]) + BFL[2] * conjg(BFL[2]) +
                            BFL[3] * conjg(BFL[3]));
                dla = sqrt(dle / dlb);
                printf(" E/B= "); /* 782 FORMAT(' E/B= ', E9.3, ' ', $) */
                fE(9, 3, 0, dla);
                printf(" ");
                printf("  \n"); /* PRINT 6 */
                break;
            case 'F':
            case 'f': /* A,F */
                printf(" "); /* 11 FORMAT(' ', 1pE14.7, 1PE10.2, '  ', $) */
                fE(14, 7, 1, creal(X));
                fE(10, 2, 1, cimag(X));
                printf("  ");
                break;
            case 'P':
            case 'p': /* P */
                printf(" "); /* 13 FORMAT(' ', F12.7, '  ', $) */
                fF(12, 7, P);
                printf("  ");
                break;
            case 'Z':
            case 'z': /* Z */
                printf(" "); /* 15 FORMAT(' ', F12.7, '  ', $) */
                fF(12, 7, Z);
                printf("  ");
                break;
            case 'E':
            case 'e': /* E */
                /* 17 FORMAT(' EX=', F7.4, F8.4, '  EY=', F7.4, F8.4,
                 *          '  EZ=', F7.4, F8.4, ' ', $) */
                printf(" EX=%7.4f%8.4f  EY=%7.4f%8.4f  EZ=%7.4f%8.4f ",
                       creal(EFL[1]), cimag(EFL[1]), creal(EFL[2]), cimag(EFL[2]),
                       creal(EFL[3]), cimag(EFL[3]));
                break;
            case 'B':
            case 'b': /* B */
                /* 19 FORMAT(' BX=', 1PE10.2, 1PE10.2, '  BY=', 1PE10.2, 1PE10.2,
                 *          '  BZ=', 1PE10.2, 1PE10.2, ' ', $) */
                printf(" BX=");
                fE(10, 2, 1, creal(BFL[1]));
                fE(10, 2, 1, cimag(BFL[1]));
                printf("  BY=");
                fE(10, 2, 1, creal(BFL[2]));
                fE(10, 2, 1, cimag(BFL[2]));
                printf("  BZ=");
                fE(10, 2, 1, creal(BFL[3]));
                fE(10, 2, 1, cimag(BFL[3]));
                printf(" ");
                break;
            case 'G':
            case 'g': /* G */
                /* 21 FORMAT(' VGP= ', 1PE9.2, '  VGZ= ', 1PE9.2, '  ', $) */
                printf(" VGP= ");
                fE(9, 2, 1, VG[1]);
                printf("  VGZ= ");
                fE(9, 2, 1, VG[2]);
                printf("  ");
                break;
            case 'S':
            case 's': /* s */
                /* 23 FORMAT(' SGP= ', 1PE9.2, '  SGZ= ', E9.2, '  ', $)
                 * (the 1P scale persists for the following E9.2) */
                printf(" SGP= ");
                fE(9, 2, 1, SG[1]);
                printf("  SGZ= ");
                fE(9, 2, 1, SG[2]);
                printf("  ");
                break;
            case 'D':
            case 'd': /* D */
                /* 25 FORMAT(1P, 'D=', 2E10.2, '  DX=', 2E10.2, '  DZ=', 2E10.2,
                 *          ' DP=', 2E10.2, /, $) */
                printf("D=");
                fE(10, 2, 1, creal(D));
                fE(10, 2, 1, cimag(D));
                printf("  DX=");
                fE(10, 2, 1, creal(DX));
                fE(10, 2, 1, cimag(DX));
                printf("  DZ=");
                fE(10, 2, 1, creal(DZ));
                fE(10, 2, 1, cimag(DZ));
                printf(" DP=");
                fE(10, 2, 1, creal(DP));
                fE(10, 2, 1, cimag(DP));
                printf("\n");
                break;
            case 'R':
            case 'r': /* R */
                /* 27 FORMAT(' RI=', 1P, 2E10.2, $) */
                printf(" RI=");
                fE(10, 2, 1, creal(RI));
                fE(10, 2, 1, cimag(RI));
                printf("  \n"); /* PRINT 6 */
                break;
            case 'U':
            case 'u': /* U */
                /* 276 FORMAT(' ene= ', 1PE10.2, $) */
                printf(" ene= ");
                fE(10, 2, 1, ENE * 2.0);
                /*      PRINT 6 */
                break;
            case 'T':
            case 't': /* T */
                /* 29 FORMAT(' E', 2I1, '=', 2(1PE10.2), '  EX', 2I1, '=', 2(1PE10.2),
                 *     'EZ', 2I1, '=', 2(1PE10.2), '  EP', 2I1, '=', 2(1PE10.2), /) */
                for (J = 1; J <= 6; J++) {
                    N = 1 + J / 4 + J / 6;
                    M = J - J / 4 * 2 - J / 6;
                    printf(" E%1d%1d=", N, M);
                    fE(10, 2, 1, creal(E[J][1]));
                    fE(10, 2, 1, cimag(E[J][1]));
                    printf("  EX%1d%1d=", N, M);
                    fE(10, 2, 1, creal(E[J][2]));
                    fE(10, 2, 1, cimag(E[J][2]));
                    printf("EZ%1d%1d=", N, M);
                    fE(10, 2, 1, creal(E[J][3]));
                    fE(10, 2, 1, cimag(E[J][3]));
                    printf("  EP%1d%1d=", N, M);
                    fE(10, 2, 1, creal(E[J][4]));
                    fE(10, 2, 1, cimag(E[J][4]));
                    /* FORMAT 29 ends with '/': terminate this record and
                     * emit the empty record begun by the trailing slash. */
                    printf("\n\n");
                }
                printf("  \n"); /* PRINT 6 */
                /* 31 FORMAT(1P, ' XX=', 12E12.3/' PP=', 6E12.3/' ZZ=', 6E12.3/)
                 * The item list is XX (complex -> 20 interleaved re/im
                 * values), PP (10) and ZZ (10) = 40 values.  The format
                 * consumes 12/6/6 values per record; after the trailing '/'
                 * the format reverts (emitting an empty record) and restarts
                 * at ' XX=' for the remaining values.  The exact record
                 * layout for the fixed 40-value list is reproduced below. */
                {
                    double items[40];
                    int cnt, k;
                    for (k = 0; k < 10; k++) {
                        items[2 * k] = creal(XX[k + 1]);
                        items[2 * k + 1] = cimag(XX[k + 1]);
                    }
                    for (k = 0; k < 10; k++)
                        items[20 + k] = PP[k + 1];
                    for (k = 0; k < 10; k++)
                        items[30 + k] = ZZ[k + 1];
                    /* record 1: ' XX=' + items[0..11] */
                    printf(" XX=");
                    for (cnt = 0; cnt < 12; cnt++)
                        fE(12, 3, 1, items[cnt]);
                    printf("\n");
                    /* record 2: ' PP=' + items[12..17] */
                    printf(" PP=");
                    for (cnt = 12; cnt < 18; cnt++)
                        fE(12, 3, 1, items[cnt]);
                    printf("\n");
                    /* record 3: ' ZZ=' + items[18..23] */
                    printf(" ZZ=");
                    for (cnt = 18; cnt < 24; cnt++)
                        fE(12, 3, 1, items[cnt]);
                    printf("\n");
                    /* reversion: trailing '/' emits an empty record */
                    printf("\n");
                    /* record 4: ' XX=' + items[24..35] */
                    printf(" XX=");
                    for (cnt = 24; cnt < 36; cnt++)
                        fE(12, 3, 1, items[cnt]);
                    printf("\n");
                    /* record 5: ' PP=' + items[36..39] */
                    printf(" PP=");
                    for (cnt = 36; cnt < 40; cnt++)
                        fE(12, 3, 1, items[cnt]);
                    printf("\n");
                }
                break;
            case ' ': /* ignore spaces */
                break;
            default:
                /* List-directed: 'char:' then char item -> no separator
                 * (character item immediately following a character item). */
                printf(" Unknown output request:%c\n", IC);
                break;
            }
        }
    }
}

void INOUT(void)
{
    int restartFormatLoop;
    for (;;) { /* output_format_loop */
        printf("#OUTPUT: \n");
        printf("\n");
        IOS = fortran_read_line(IOU, 20);
        if (IOS != 0) {
            memset(IOU, ' ', 20);
            IOU[20] = '\0';
        }
        restartFormatLoop = 0;
        for (K = 1; K <= 20; K++) { /* loop through all input characters */
            IC = IOU[K - 1];
            switch (IC) {
            case 'H':
            case 'h':
            case '?': /* help */
                /* 105 FORMAT: print output options, defined in the beginning of file */
                printf(" The output is determined by a string of letters:\n");
                /* 'letters:'// -> the double slash emits an empty record */
                printf("\n");
                printf(" b     wave magnetic field components.\n");
                printf(" d     dispersion function and derivatives.\n");
                printf(" e     wave electric field components.\n");
                printf(" f     frequency <real,imaginery>.\n");
                printf(" g     group velocity components.\n");
                printf(" h     |e|/|b| [mV/nT].\n");
                printf(" l     |bp|/|bz|.\n");
                printf(" m     Im[bx]/Re[by].\n");
                printf(" n     ellipticity bx/by.\n");
                printf(" o     ellipticity general\n");
                printf(" p     perpendicular component of wave vector.\n");
                printf(" r     refractive index.\n");
                printf(" s     spatial growth-rates.\n");
                printf(" t     dielectric tensor and derivatives.\n");
                printf(" z     z-component of wave vector.\n");
                printf(" u     total wave energy / energy in electric field.\n");
                printf(" v     Poynting flux (in uW/m2 for <E^2>=0.5(mV/m)^2).\n");
                printf(" x     phase of bz against bx (/+/ means bz is in front of bx).\n");
                printf(" y     energy density and flux of each plasma component.\n");
                printf(" The results are  normally printed on one line in the order\n");
                printf(" they are specified. A new line is obtained by inserting\n");
                printf(" a \"/\" in the string.\n");
                printf(" Example: output: pzf/e\n");
                printf(" The wave numbers and the frequency are printed on one line,\n");
                printf(" and the electric field components on the next.\n");
                /* '...on the next.'//) -> trailing double slash: one empty
                 * record, and the record itself is also terminated. */
                printf("\n\n");
                /* Fortran: CYCLE output_format_loop -> restart the OUTER loop */
                restartFormatLoop = 1;
                break;
            case ' ':
                break;
            case '+': /* print debug info */
                printDebugInfo = 1;
                IOU[K - 1] = ' ';
                break;
            case '-': /* remove printing debug info */
                printDebugInfo = 0;
                IOU[K - 1] = ' ';
                break;
            default: /* any character */
                KMX = K;
                break;
            }
        }
        if (restartFormatLoop)
            continue; /* CYCLE output_format_loop */
        break; /* EXIT output_format_loop */
    }
}
