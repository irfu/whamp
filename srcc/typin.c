/*
 * typin.c -- translation of SUBROUTINE TYPIN (src/typin.f90)
 *
 *        ARGUMENTS: NPL    NEW PLASMA. WHEN A PLASMA PARAMETER IS
 *                          CHANGED, NPL IS SET TO 1.
 *                   KFS    =1  P SPECIFIED LAST.
 *                          =2  Z SPECIFIED LAST.
 */
#include "whamp.h"
#include "comin.h"

/* Emulate Fortran READ(*,'(A)') (defined in output.c) */
int fortran_read_line(char *buf, int len);

/* Fortran declares IOUT/KP/KZ with initialization, giving them SAVE
 * semantics across calls. */
static int IOUT = 2, KP = 1, KZ = 1;

/* Emulate Fortran INDEX(string, substring) for a single character */
static int findex(const char *str, char c)
{
    const char *p = strchr(str, c);
    return p ? (int)(p - str) + 1 : 0;
}

static void print_help(void)
{
    /* 131 FORMAT */
    printf(" AN INPUT LINE MAY CONSIST OF UP TO 80 CHARACTERS.\n");
    printf(" THE FORMAT IS:\n");
    printf(" NAME1=V11,V12,V13,...NAME2=V21,V22,...NAME\n");
    printf(" THE NAMES ARE CHOSEN FROM THE LIST:\n");
    printf("\n");
    printf(" NAME              PARAMETER\n");
    printf(" A(I)              THE ALPHA1 PARAMETER IN THE DISTRIBUTION.\n");
    printf("                   (I) IS THE COMPONENT NUMBER, I=1 - 6.\n");
    printf(" B(I)              THE ALPHA2 PARAMETER IN THE DISTRIBUTION.\n");
    printf(" C                 THE ELECTRON CYCLOTRON FREQ. IN KHZ.\n");
    printf(" D(I)              THE DELTA PARAMETER IN THE DISTRIBUTION\n");
    printf(" F                 FREQUENCY, START VALUE FOR ITERATION.\n");
    printf(" L            L=1  THE P AND Z PARAMETERS ARE INTERPRETED\n");
    printf("                   AS LOGARITHMS OF THE WAVE NUMBERS. THIS\n");
    printf("                   OPTION ALLOWS FOR LOGARITHMIC STEPS.\n");
    printf("              L=0  DEFAULT VALUE. LINEAR STEPS.\n");
    printf(" M(I)              MASS IN UNITS OF PROTON MASS.\n");
    printf(" N(I)              NUMBER DENSITY IN PART./CUBIC METER\n");
    printf(" P(I)              PERPENDICULAR WAVE VECTOR COMPONENTS.\n");
    printf("                   P(1) IS THE SMALLEST VALUE, P(2) THE\n");
    printf("                   LARGEST VALUE, AND P(3) THE INCREMENT.\n");
    printf(" S                 STOP! TERMINATES THE PROGRAM.\n");
    /* 132 FORMAT */
    printf(" T(I)              TEMPERATURE IN KEV\n");
    printf(" V(I)              DRIFT VELOCITY / THERMAL VELOCITY.\n");
    printf(" Z(I)              Z-COMPONENT OF WAVE VECTOR. I HAS THE\n");
    printf("                   SAME MEANING AS FOR P(I).\n");
    printf(" A NAME WITHOUT INDEX REFERS TO THE FIRST ELEMENT, \"A\" IS \n");
    printf(" THUS EQUIVALENT TO \"A(1)\". THE VALUES V11,V12,.. MAY BE \n");
    printf(" SPECIFIED IN I-, F-, OR E-FORMAT, SEPARATED BY COMMA(,).\n");
    printf(" THE \"=\" IS OPTIONAL, BUT MAKES THE INPUT MORE READABLE.\n");
    printf(" EXAMPLE: INPUT:A1.,2. B(3).5,P=.1,.2,1.E-2\n");
    printf(" THIS SETS A(1)=1., A(2)=2., B(3)=.5, P(1)=.1, P(2)=.2,\n");
    printf(" AND P(3)=.01. IF THE INCREMENT P(3)/Z(3) IS NEGATIVE, P/Z\n");
    printf(" WILL FIRST BE SET TO P(2)/Z(2) AND THEN STEPPED DOWN TO\n");
    printf(" P(1)/Z(1)\n");
    printf(" THE LAST SPECIFIED OF P AND Z WILL VARY FIRST.\n");
    printf(" IF THE LETTER \"O\" (WITHOUT VALUE) IS INCLUDED, YOU WILL\n");
    printf("  BE ASKED TO SPECIFY A NEW OUTPUT FORMAT.\n");
    /* FORMAT ends with '/' after the last string -> one extra empty record */
    printf("\n");
}

void TYPIN(int *NPL, int *KFS)
{
    int IE = 1, IOF = 1, IOS, NC, NV = 0;
    int flagSuccessReading, flagAmbiguousCharacter, flagTooLongNumber, flagNextNumber;
    int helpRequested;
    float DEC, DEK;
    const char UPP_LTRS[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const char DIGITS[] = "0123456789";
    const char LOW_LTRS[] = "abcdefghijklmnopqrstuvwxyz";
    char INP[81], IC, inputVariable;
    float TV[3]; /* Fortran TV(2), 1-based: TV[1], TV[2] */
    float inputNumber;

    *NPL = 0; /* default is that plasma parameters are not changed */
    flagSuccessReading = 0; /* default value that did not succeed to read the line */
    for (;;) {              /* input_loop */
        inputVariable = ' '; /* default empty */
        IOS = 1;

        for (;;) {
            if (IOS == 0)
                break;
            printf("#INPUT: \n");
            printf("  \n");
            IOS = fortran_read_line(INP, 80);
            if (IOS != 0) {
                /* REWIND IFILE (unit 5): no-op for C stdin */
            }
        }
        NC = 0;
        helpRequested = 0;
        for (;;) { /* loop_input_line */
            flagAmbiguousCharacter = 0; /* default assume non-alpha character */
            flagTooLongNumber = 0;      /* default assume input numbers are not too long */
            flagNextNumber = 0; /* default that we enter the first number of variable */
            NC = NC + 1;
            if (NC <= 80) { /* still within input line 80 character limit */
                IC = INP[NC - 1];
                if (IC == ' ')
                    continue;
                if (inputVariable == ' ') {
                    flagNextNumber = 1;
                } else {
                    if (findex(UPP_LTRS, IC) > 0 || findex(LOW_LTRS, IC) > 0) {
                        /* IC is alphabetic character: fall through */
                    } else if (findex(DIGITS, IC) > 0) {
                        /* IC is digit */
                        TV[IE] = TV[IE] * DEK + (float)(findex(DIGITS, IC) - 1) * DEC;
                        DEC = DEC * DEK / 10.;
                        NV = 1;
                        continue;
                    } else if (IC != ',') { /* strange character */
                        if (IC == '-') {
                            if (TV[IE] != 0.) {
                                flagAmbiguousCharacter = 1;
                                break;
                            }
                            DEC = -DEC;
                            continue;
                        }

                        if (IC == '.') {
                            if (fabs(DEC) != 1.) {
                                flagAmbiguousCharacter = 1;
                                break;
                            }
                            DEK = DEK / 10.;
                            DEC = DEC / 10.;
                            continue;
                        }

                        if (IC != ')')
                            continue;
                        if (DEC != 1.) {
                            flagAmbiguousCharacter = 1;
                            break;
                        }
                        if (IE != 1 || TV[1] <= 0.) {
                            flagAmbiguousCharacter = 1;
                            break;
                        }
                        IOF = (int)(TV[1] + .1);
                        TV[1] = 0.;
                        NV = 0;
                        continue;
                    }
                }
            }

            if (flagNextNumber == 0) { /* still working on the same number input */
                if (NV == 1) {
                    if (IOF > 10) { /* more than 10 plasma components */
                        flagTooLongNumber = 1;
                        break;
                    }
                    if (IOF > 3 && (findex("cflpz", inputVariable) > 0)) {
                        /* too many components for variable */
                        flagTooLongNumber = 1;
                        break;
                    }
                    if (IOF > 1 && (findex("cfl", inputVariable) > 0)) {
                        /* scalars have more than 1 component */
                        flagTooLongNumber = 1;
                        break;
                    }
                    if (IC == 'E' || IC == 'e') {
                        IE = 2;
                        DEK = 10.;
                        DEC = 1.;
                        continue;
                    }
                    inputNumber = TV[1] * (float)pow(10., (double)TV[2]);

                    switch (inputVariable) {
                    case 'A':
                    case 'a':
                        AA[IOF][1] = inputNumber;
                        break;
                    case 'B':
                    case 'b':
                        AA[IOF][2] = inputNumber;
                        break;
                    case 'C':
                    case 'c':
                        XC = inputNumber;
                        break;
                    case 'D':
                    case 'd':
                        DD[IOF] = inputNumber;
                        break;
                    case 'F':
                    case 'f':
                        XOI = inputNumber;
                        break;
                    case 'L':
                    case 'l':
                        PZL = inputNumber;
                        break;
                    case 'M':
                    case 'm':
                        ASS[IOF] = inputNumber;
                        break;
                    case 'N':
                    case 'n':
                        DN[IOF] = inputNumber;
                        break;
                    case 'P':
                    case 'p':
                        PM[IOF] = inputNumber;
                        break;
                    case 'T':
                    case 't':
                        TA[IOF] = inputNumber;
                        break;
                    case 'V':
                    case 'v':
                        VD[IOF] = inputNumber;
                        break;
                    case 'Z':
                    case 'z':
                        ZM[IOF] = inputNumber;
                        break;
                    }
                    IOF = IOF + 1;
                    if (inputVariable == 'p')
                        KP = IOF;
                    if (inputVariable == 'z')
                        KZ = IOF;
                }
            }

            if (NC > 80) {
                flagSuccessReading = 1;
                break;
            }
            DEK = 10.;
            DEC = 1.;
            TV[1] = 0.;
            TV[2] = 0.;
            NV = 0;
            IE = 1;
            if (IC == ',')
                continue;
            inputVariable = ' ';
            IOF = 1;
            if (IC == 'A' || IC == 'a')
                inputVariable = 'a';
            if (IC == 'B' || IC == 'b')
                inputVariable = 'b';
            if (IC == 'C' || IC == 'c')
                inputVariable = 'c';
            if (IC == 'D' || IC == 'd')
                inputVariable = 'd';
            if (IC == 'F' || IC == 'f')
                inputVariable = 'f';
            if (IC == 'H' || IC == 'h') {
                print_help();
                helpRequested = 1; /* cycle input_loop */
                break;
            }
            if (IC == 'L' || IC == 'l')
                inputVariable = 'l';
            if (IC == 'M' || IC == 'm')
                inputVariable = 'm';
            if (IC == 'N' || IC == 'n')
                inputVariable = 'n';
            if (IC == 'O' || IC == 'o') {
                IOUT = 0;
                continue;
            }
            if (IC == 'P' || IC == 'p')
                inputVariable = 'p';
            if (IC == 'S' || IC == 's')
                exit(0); /* STOP */
            if (IC == 'T' || IC == 't')
                inputVariable = 't';
            if (IC == 'V' || IC == 'v')
                inputVariable = 'v';
            if (IC == 'Z' || IC == 'z')
                inputVariable = 'z';

            if (inputVariable == ' ')
                break;
            if (findex("abcdmntv", inputVariable) > 0)
                *NPL = 1; /* plasma parameters changed */
            if (inputVariable == 'p')
                KP = 1;
            if (inputVariable == 'p')
                *KFS = 1;
            if (inputVariable == 'z')
                KZ = 1;
            if (inputVariable == 'z')
                *KFS = 2;
        }
        if (helpRequested)
            continue; /* cycle input_loop */
        if (flagAmbiguousCharacter == 1) {
            printf(" AMBIGUITY CAUSED BY THE CHARACTER \"%c\"\n", IC);
            printf(" THE REST OF THE LINE IS IGNORED. PLEASE TRY AGAIN!\n");
            continue;
        }

        if (flagTooLongNumber == 1) {
            TV[1] = TV[1] * (float)pow(10., (double)TV[2]);
            printf(" THE VALUE");
            fE(11, 3, 0, (double)TV[1]);
            printf(" WILL NOT FIT IN THE VARIABLE FIELD\n");
            printf(" THE REST OF THE LINE IS IGNORED. PLEASE TRY AGAIN!\n");
            continue;
        }
        if (flagSuccessReading == 0) {
            IOS = 1;
            for (;;) {
                if (IOS == 0)
                    break;
                printf("HELP, YES OR NO?");
                {
                    char tmp[81];
                    IOS = fortran_read_line(tmp, 80);
                    IC = tmp[0];
                    if (IOS != 0) {
                        /* REWIND IFILE: no-op for C stdin */
                    }
                }
            }
        }
        if (IC == 'N' || IC == 'n')
            continue;
        if (flagSuccessReading == 0) {
            print_help();
            continue;
        }
        if (fabs(XOI) <= 0.) {
            printf(" START FREQUENCY");
            scanf("%lf", &XOI);
        }
        if (KP == 1) {
            printf(" PERP. WAVE VECTOR UNDEFINED!");
            continue;
        }

        if (KP == 2) {
            PM[2] = PM[1];
            KP = 3;
        }
        if (KP == 3) {
            PM[3] = PM[2] - PM[1];
            KP = 4;
        }
        if (KP == 4) {
            if (PM[3] == 0)
                PM[3] = 10;
        }

        if (KZ == 1) {
            printf(" PARALLEL WAVE VECTOR UNDEFINED!");
            continue;
        }

        if (KZ == 2) {
            ZM[2] = ZM[1];
            KZ = 3;
        }
        if (KZ == 3) {
            ZM[3] = ZM[2] - ZM[1];
            KZ = 4;
        }
        if (KZ == 4) {
            if (ZM[3] == 0)
                ZM[3] = 10;
        }
        break;
    }

    if (IOUT != 1)
        INOUT();
    IOUT = 1;
    return;
}
