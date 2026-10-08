/*
 * whamp.c -- translation of PROGRAM WHAMP (src/whamp.f90)
 *
 * The Fortran host association of the contained subroutines
 * (print_plasma_parameters, root_finding) is translated by making the
 * program variables file-scope statics.
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"

static int LoopI, IERR, IRK, J, KFS;
static int rootFindingConverged;
static int solutionIsTooHeavilyDamped;
static int isChangedPlasmaModel;
static double REN[11]; /* particle mass expressed in masses of first particles */
static double RN;      /* mass of first particle in electron masses */
static double ADIR;    /* abs(D) */
static double DEK, DKP, DKZ, KV, PFQ, PLG, PLO, PO, PVO;
static double RED, TR, XA, XI, ZLG, ZLO, ZO, ZVO;
static cd XO, XVO, ddDX, OME, FPX, DOX, DOZ, DOP;
static double T[11], ST[11];

static void print_plasma_parameters(void);
static void root_finding(void);

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
    int argi;

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
                    char trimmed[21];
                    int n2 = (int)strlen(adj);
                    (void)n2;
                    strcpy(trimmed, adj);
                    printf(" Option '");
                    printf("%s", trimmed);
                    printf("' is unknown\n");
                }
            }
        }
        iarg = iarg + 1;
        (void)argi;
    }

    IERR = 0;

loop_plasma_update:
    for (;;) {
        isChangedPlasmaModel = 0; /* changed to .true. in code when new plasma
                                     parameters are entered */
        DEN = 0.;
        RED = 0.;
        for (J = 1; J <= 10; J++) { /* loop_species */
            REN[J] = F32(1836.1) * ASS[J];
            if (REN[J] == 0.)
                REN[J] = 1.;
            T[J] = TA[J] / TA[1];
            if (DN[J] == 0.)
                continue;
            JMA = J;
            RED = RED + DN[J] / REN[J];
            if (ASS[J] == 0.)
                DEN = DEN + DN[J];
        }

        RN = REN[1];
        /*                  ****  NORMALIZED TEMPERATURES AND VELOCITIES.  **** */
        for (J = 1; J <= JMA; J++) {
            REN[J] = REN[J] / RN;
            T[J] = T[J] * REN[J];
            ST[J] = sqrt(T[J]);
        }

        DEK = 12405.;
        PFQ = RED / DEK;
        PX = sqrt(PFQ);
        XA = XC / RN;
        TR = TA[1] / RN;
        CV = TR * (1022. + TR) / ((511. + TR) * (511. + TR));
        CV = 1. / sqrt(CV);
        DEK = DEK * RN;

        print_plasma_parameters();
        /*                  ****  ASK FOR INPUT!  **** */
    loop_typin:
        for (;;) {
            /* for new plasma skip calling typin until convergence checked */
            if (!isChangedPlasmaModel)
                TYPIN(&isChangedPlasmaModel, &KFS);
            if (isChangedPlasmaModel)
                goto cycle_loop_plasma_update;
            isChangedPlasmaModel = 0;
            KV = 1;
            PLG = PM[1];
            if (PM[3] < 0.)
                PLG = PM[2];
            ZLG = ZM[1];
            if (ZM[3] < 0.)
                ZLG = ZM[2];
            if (PZL == 1.) {
                P = pow(10., PLG);
                Z = pow(10., ZLG);
            } else {
                P = PLG;
                Z = ZLG;
            }
            X = XOI;
            for (;;) { /* loop_z_p */
                OME = (X * XA) * (X * XA);
                FPX = PFQ / OME;
                for (J = 1; J <= JMA; J++) {
                    XX[J] = X * REN[J];
                    PP[J] = P * ST[J];
                    ZZ[J] = Z * ST[J];
                    XP[J] = DN[J] / DEK / REN[J] / OME;
                }

                solutionIsTooHeavilyDamped = 0; /* default not heaviily damped */
                rootFindingConverged = 0;       /* default no convergence */
                root_finding();

                if (!rootFindingConverged) {
                    /* 125 FORMAT(2X,'NO CONVERGENCE!'/'  KP=',F6.3,'  KZ=',F6.4,
                       '  X=',E12.2,E12.2/'  I=',I3,'  IRK=',I3/) */
                    printf("  NO CONVERGENCE!\n");
                    printf("  KP=");
                    fF(6, 3, P);
                    printf("  KZ=");
                    fF(6, 4, Z);
                    printf("  X=");
                    fE(12, 2, 0, creal(X));
                    fE(12, 2, 0, cimag(X));
                    printf("\n");
                    printf("  I=%3d  IRK=%3d\n", LoopI, IRK);
                    printf("\n");
                    if (KFS == 1)
                        PLG = 1e99; /* end cycling in P */
                    if (KFS == 2)
                        ZLG = 1e99; /* end cycling in Z */
                }
                if (rootFindingConverged && !solutionIsTooHeavilyDamped) {
                    /*                  ****  CONVERGENCE!  **** */
                    DIFU(4, JMA, &IERR);

                    XI = cimag(X);
                    VG[1] = -creal(DP / DX);
                    VG[2] = -creal(DZ / DX);
                    RI = sqrt(P * P + Z * Z) * CV / X; /* refractive index */
                    if (VG[1] != 0.)
                        SG[1] = XI / VG[1];
                    if (VG[2] != 0.)
                        SG[2] = XI / VG[2];
                    /*          ****  PRINT THE RESULTS.  **** */
                    OUTPT();
                    PO = P;
                    ZO = Z;
                    XO = X;
                    if (KV != 0) {
                        XVO = X;
                        ZVO = Z;
                        ZLO = ZLG;
                        PVO = P;
                        PLO = PLG;
                        DOX = DX;
                        DOZ = DZ;
                        DOP = DP;
                        KV = 0;
                    }
                }
                if (solutionIsTooHeavilyDamped) {
                    printf("  TOO HEAVILY DAMPED!\n");
                    printf("     \n");
                    IERR = 0;
                    OUTPT();
                    if (KFS == 1)
                        PLG = 1e99;
                    if (KFS == 2)
                        ZLG = 1e99;
                }

                if (KFS == 1) { /* cycle first P */
                    PLG = PLG + PM[3];
                    /*                   **** UPDATE P AND Z.  **** */
                    if (PLG >= PM[1] && PLG <= PM[2]) {
                        P = PLG + PZL * (pow(10., PLG) - PLG);
                    } else {
                        ZLG = ZLG + ZM[3];
                        printf("\n"); /* print * (no item list) */
                        if (ZLG < ZM[1] || ZLG > ZM[2])
                            goto cycle_loop_typin;
                        KV = 1;
                        PLG = PLO;
                        P = PVO;
                        Z = ZLG + PZL * (pow(10., ZLG) - ZLG);
                    }
                } else if (KFS == 2) { /* cycle first Z */
                    ZLG = ZLG + ZM[3];
                    if (ZLG >= ZM[1] && ZLG <= ZM[2]) {
                        Z = ZLG + PZL * (pow(10., ZLG) - ZLG);
                    } else {
                        PLG = PLG + PM[3];
                        printf("\n"); /* print * (no item list) */
                        if (PLG < PM[1] || PLG > PM[2])
                            goto cycle_loop_typin;
                        KV = 1;
                        ZLG = ZLO;
                        Z = ZVO;
                        P = PLG + PZL * (pow(10., PLG) - PLG);
                    }
                }
                /*                    ****  NEW START FREQUENCY.  **** */
                if (KV != 0) {
                    DKP = P - PVO;
                    DKZ = Z - ZVO;
                    ddDX = (DKP * DOP + DKZ * DOZ) / DOX;
                    X = XVO - ddDX;
                } else {
                    DKP = P - PO;
                    DKZ = Z - ZO;
                    ddDX = (DKP * DP + DKZ * DZ) / DX;
                    X = XO - ddDX;
                }
                continue; /* loop_z_p */
            }
        cycle_loop_typin:;
        }
    cycle_loop_plasma_update:;
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

static void root_finding(void)
{
    /* Newton's iteration method with small adjustments:
     * 1) convergence criteria on relative and absolute size of CX
     * 2) convergence criteria on relative change in D */
    cd CX; /* correction in Newton's iteration method */

    DIFU(2, JMA, &IERR);
    if (IERR != 0)
        solutionIsTooHeavilyDamped = 1;
    /*                  ****  START OF ITERATION.  **** */
    if (printDebugInfo) {
        printf(" START:. X= ");
        ld_complex(X);
        printf(" D= ");
        ld_complex(D);
        printf(" DX= ");
        ld_complex(DX);
        printf("\n");
    }
    for (LoopI = 1; LoopI <= maxIterations; LoopI++) { /* loop_iteration */
        ADIR = cabs(D);
        IRK = 1;
        CX = D / DX;
        /*CX=CX*X/(2*CX+X) ! finding zero of (w^2 D), faster convergence */
        for (;;) { /* irk_loop */
            X = X - CX;
            OME = (X * XA) * (X * XA);
            FPX = PFQ / OME;
            for (J = 1; J <= JMA; J++) {
                XP[J] = DN[J] / DEK / REN[J] / OME;
                XX[J] = X * REN[J];
            }
            DIFU(2, JMA, &IERR);
            if (printDebugInfo) {
                /* write(*,'(I2,A ,I2 ,A,2E16.8,A,2E16.8 ,A,2E16.8,A,2E16.8)')
                 * I2=' ',A='.',I2='.',A=' X=' then 2E16.8 (P=0) per complex */
                printf("%2d.%2d. X=", LoopI, IRK);
                fE(16, 8, 0, creal(X));
                fE(16, 8, 0, cimag(X));
                printf(" CX=");
                fE(16, 8, 0, creal(CX));
                fE(16, 8, 0, cimag(CX));
                printf(" D=");
                fE(16, 8, 0, creal(D));
                fE(16, 8, 0, cimag(D));
                printf(" DX=");
                fE(16, 8, 0, creal(DX));
                fE(16, 8, 0, cimag(DX));
                printf("\n");
            }
            if (IERR != 0) {
                solutionIsTooHeavilyDamped = 1;
                goto end_loop_iteration; /* exit loop_iteration */
            }
            if (cabs(D) < ADIR) {
                if ((cabs(CX) <= F32(1.E-6) * cabs(X)) /* relative frequency precision */
                    || (cabs(CX) < F32(1e-6))) {  /* absolute precision */
                    rootFindingConverged = 1;
                    if (LoopI >= 2) { /* at least 2 steps have been made */
                        goto end_loop_iteration; /* exit loop_iteration */
                    } else {
                        goto cycle_loop_iteration; /* cycle loop_iteration */
                    }
                } else {
                    goto cycle_loop_iteration; /* cycle loop_iteration */
                }
            } else {
                X = X + CX;
                CX = CX / 2.;
                if (IRK > 3) { /* check for sitting at local minima */
                    if ((cabs(D) - ADIR) / ADIR < F32(1e-5)) {
                        printf("   Local minima!\n"); /* print *,' Local minima!' */
                        rootFindingConverged = 0;
                        goto end_loop_iteration; /* exit loop_iteration */
                    }
                }
                IRK = IRK + 1;
                if (IRK > maxIterations)
                    goto end_loop_iteration; /* exit loop_iteration */
            }
        }
    cycle_loop_iteration:;
    }
end_loop_iteration:;
}
