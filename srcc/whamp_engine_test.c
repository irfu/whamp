/*
 * whamp_engine_test.c -- translation of PROGRAM WHAMP (src/whamp_engine_test.f90)
 *
 * Same host program as whamp.c but driving whamp_engine() instead of the
 * inline solver, and dumping the engine output arrays after each solve.
 *
 * The Fortran host association of the contained subroutines
 * (print_plasma_parameters, species_symbol) is translated by making them
 * file-scope statics.
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"
#include "comoutput.h"

static int J, KFS;
static int isChangedPlasmaModel;

static void print_plasma_parameters(void);
static void species_symbol(double mass, char *symbol);

/* Fortran: pure function species_symbol(mass) result(symbol)
 * character(5) :: symbol */
static void species_symbol(double mass, char *symbol)
{
    if (mass == 0)
        strcpy(symbol, "e-");
    else if (mass == 1)
        strcpy(symbol, "H+");
    else if (mass == 2)
        strcpy(symbol, "He++");
    else if (mass == 4)
        strcpy(symbol, "He+");
    else if (mass == 16)
        strcpy(symbol, "O+");
    else
        sprintf(symbol, "m=%3d", (int)fint(mass)); /* write(symbol,'(a,I3)') 'm=',mass */
}

int main(int argc, char **argv)
{
    int narg, iarg;
    char inputParameter[21], modelFilename[21];

    /* Default plasma model */
    DN[1] = 1.0e6;
    DN[2] = 1.0e6;
    DN[3] = DN[4] = DN[5] = DN[6] = DN[7] = DN[8] = DN[9] = DN[10] = 0.0;
    TA[1] = F32(0.01);
    TA[2] = F32(0.001);
    TA[3] = TA[4] = TA[5] = TA[6] = TA[7] = TA[8] = TA[9] = TA[10] = 0.0;
    DD[1] = 1.0;
    DD[2] = 1.0;
    DD[3] = DD[4] = DD[5] = DD[6] = DD[7] = DD[8] = DD[9] = DD[10] = 0.0;
    AA[1][1] = 5.0;
    AA[2][1] = 1.0;
    for (J = 3; J <= 10; J++)
        AA[J][1] = 0.0;
    AA[1][2] = F32(0.1);
    AA[2][2] = F32(0.1);
    for (J = 3; J <= 10; J++)
        AA[J][2] = 0.0;
    ASS[1] = 16.0;
    for (J = 2; J <= 10; J++)
        ASS[J] = 0.0;
    VD[1] = 1.0;
    for (J = 2; J <= 10; J++)
        VD[J] = 0.0;
    XC = F32(2.79928);
    PZL = 0.0;
    cycleZFirst = 1;
    PM[1] = 0.0;
    PM[2] = 0.0;
    PM[3] = 10.0;
    ZM[1] = 0.0;
    ZM[2] = 0.0;
    ZM[3] = 10.0;
    XOI = F32(.1);

    /* Check command line input parameters */
    narg = argc - 1;
    iarg = 1;
    for (;;) {
        if (iarg > narg)
            break;
        /* get_command_argument(iArg, inputParameter): pad/truncate to 20 */
        {
            const char *a = argv[iarg];
            size_t n = strlen(a);
            if (n > 20)
                n = 20;
            memcpy(inputParameter, a, n);
            memset(inputParameter + n, ' ', 20 - n);
            inputParameter[20] = '\0';
        }
        if (printDebugInfo) {
            printf(" Input parameter:");
            {
                int k; /* write(*,*) pads character to length 20 */
                for (k = 0; k < 20; k++)
                    putchar(inputParameter[k]);
            }
            putchar('\n');
        }
        /* select case (adjustl(inputParameter)) -- Fortran string compare
         * ignores trailing blanks */
        {
            char adj[21];
            int k = 0, s = 0, n;
            while (inputParameter[s] == ' ' && s < 20)
                s++;
            while (s + k < 20 && k < 20) {
                adj[k] = inputParameter[s + k];
                k++;
            }
            n = k;
            while (n > 0 && adj[n - 1] == ' ')
                n--;
            adj[n] = '\0';

            if (strcmp(adj, "-help") == 0 || strcmp(adj, "-h") == 0 ||
                strcmp(adj, "--help") == 0) {
                printf(" usage: whamp [-help] [-debug] [-maxiterations <number>] "
                       "[-file <modelFilename>] \n");
                return 0;
            } else if (strcmp(adj, "-debug") == 0) {
                printDebugInfo = 1;
                if (printDebugInfo)
                    printf(" Enable debugging\n");
            } else if (strcmp(adj, "-file") == 0) {
                if (iarg == narg) {
                    printf(" ERROR: File name not given\n");
                    return 0;
                }
                iarg = iarg + 1;
                {
                    const char *a = argv[iarg];
                    size_t len = strlen(a);
                    if (len > 20)
                        len = 20;
                    memcpy(modelFilename, a, len);
                    memset(modelFilename + len, ' ', 20 - len);
                    modelFilename[20] = '\0';
                }
                if (printDebugInfo) {
                    printf(" Reading file: ");
                    {
                        int k2;
                        for (k2 = 0; k2 < 20; k2++)
                            putchar(modelFilename[k2]);
                    }
                    putchar('\n');
                }
                read_input_file(modelFilename);
            } else if (strcmp(adj, "-maxiterations") == 0) {
                if (iarg == narg) {
                    printf(" ERROR: maxiterations is not specified\n");
                    return 0;
                }
                iarg = iarg + 1;
                {
                    const char *a = argv[iarg];
                    size_t len = strlen(a);
                    if (len > 20)
                        len = 20;
                    memcpy(inputParameter, a, len);
                    memset(inputParameter + len, ' ', 20 - len);
                    inputParameter[20] = '\0';
                }
                sscanf(inputParameter, "%d", &maxIterations);
                if (printDebugInfo) {
                    printf(" Max iterations: ");
                    ld_int(maxIterations);
                    putchar('\n');
                }
            } else {
                if (printDebugInfo) {
                    printf(" Option '");
                    printf("%s", adj);
                    printf("' is unknown\n");
                }
            }
        }
        iarg = iarg + 1;
    }

    whamp_engine(); /* call to get plasma parameters only */

    for (;;) { /* plasma_update_loop */
        print_plasma_parameters();
        isChangedPlasmaModel = 0;
        /*                  ****  ASK FOR INPUT!  **** */
        for (;;) { /* typin_loop */
            /* for new plasma skip calling typin until convergence checked */
            if (!isChangedPlasmaModel)
                TYPIN(&isChangedPlasmaModel, &KFS);
            if (KFS == 1)
                cycleZFirst = 0;
            if (KFS == 2)
                cycleZFirst = 1;
            if (isChangedPlasmaModel)
                break; /* cycle plasma_update_loop */
            isChangedPlasmaModel = 0;
            whamp_engine();
            /* write (*, *) 'PM=', PM -- ld_real already emits the 5 trailing
             * pad spaces that act as the list-directed separator. */
            printf(" PM=");
            ld_real(PM[1]);
            ld_real(PM[2]);
            ld_real(PM[3]);
            printf("\n");
            /* write (*, *) 'cycleZfirst=', cycleZfirst */
            printf(" cycleZfirst=");
            ld_int(cycleZFirst);
            printf("\n");
            /* write (*, *) 'fOUT=', fOUT -- list-directed over the whole
             * allocated array in Fortran (column-major) order.  The C buffer
             * is already stored column-major, so walk it linearly with the
             * 0-based F2D offset.  ld_complex's internal padding supplies the
             * list-directed separator. */
            printf(" fOUT=");
            {
                int n = kperpSize * kparSize;
                int i;
                /* Each complex item is preceded by a single list-directed
                 * separator space (a complex item starts with '(' which is
                 * non-blank, unlike a right-justified real field). */
                for (i = 0; i < n; i++) {
                    putchar(' ');
                    ld_complex(fOUT[i]);
                }
            }
            printf("\n");
        }
    }
    return 0;
}

static void print_plasma_parameters(void)
{
    char symbol[8];
    /* 101 FORMAT('# PLASMA FREQ.:',F11.4,'KHZ GYRO FREQ.:',F10.4,'KHZ   ',
       'ELECTRON DENSITY:',1PE11.5,'M-3') */
    printf("# PLASMA FREQ.:");
    fF(11, 4, PX);
    printf("KHZ GYRO FREQ.:");
    fF(10, 4, XC);
    printf("KHZ   ELECTRON DENSITY:");
    fE(11, 5, 1, DEN);
    printf("M-3\n");
    for (J = 1; J <= JMA; J++) {
        /* 102 FORMAT('# ',A3,'  DN=',1PE12.5,'  T=',0PF9.5,'  D=',F4.2,
           '  A=',F4.2,'  B=',F4.2,' VD=',F5.2) */
        species_symbol(ASS[J], symbol);
        /* A3: left-justified, blank padded to width 3 */
        printf("# %-3.3s  DN=", symbol);
        fE(12, 5, 1, DN[J]);
        printf("  T=");
        fF(9, 5, TA[J]);
        printf("  D=");
        fF(4, 2, DD[J]);
        printf("  A=");
        fF(4, 2, AA[J][1]);
        printf("  B=");
        fF(4, 2, AA[J][2]);
        printf(" VD=");
        fF(5, 2, VD[J]);
        printf("\n");
    }
}
