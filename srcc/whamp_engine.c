/*
 * whamp_engine.c -- translation of SUBROUTINE WHAMP_ENGINE (src/whamp_engine.f90)
 *
 * The Fortran contained subroutines (plasma_setup, input_setup,
 * root_finding, allocate_output_matrices, save_output) access the host
 * variables through host association; they are translated as static
 * file-scope functions with the host variables as file-scope statics.
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"
#include "comoutput.h"

/* Host variables of WHAMP_ENGINE (Fortran locals, re-created on each
 * call; zero-initialised here which is deterministic and matches the
 * Fortran behaviour for all defined uses). */
static int LoopI, IERR, IRK, J;
static int rootFindingConverged;
static int solutionIsTooHeavilyDamped;
static double REN[11]; /* particle mass expressed in masses of first particles */
static double RN;      /* mass of first particle in electron masses */
static double ADIR;    /* abs(D) */
static double DEK, DKP, DKZ, KV, PFQ, PLG, PLO, PO, PVO;
static double RED, TR, XA, XI, ZLG, ZLO, ZO, ZVO;
static cd XO, XVO, ddDX, OME, FPX, DOX, DOZ, DOP;
static double T[11], ST[11];
static double coef_poynt; /* coefficient used to simplify Poynting flux estimate */
static const double PI = 3.14159265358979;

static void plasma_setup(void)
{
    DEN = 0.;
    RED = 0.;
    for (J = 1; J <= 10; J++) {
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
}

static void input_setup(void)
{
    /* empty in the Fortran source */
}

static void root_finding(void)
{
    cd CX; /* correction in Newton's iteration method */

    DIFU(2, JMA, &IERR);
    if (IERR != 0)
        solutionIsTooHeavilyDamped = 1;
    /*                  ****  START OF ITERATION.  **** */
    for (LoopI = 1; LoopI <= maxIterations; LoopI++) {
        ADIR = cabs(D);
        IRK = 1;
        CX = D / DX;
        /*CX=CX*X/(2*CX+X) ! finding zero of (w^2 D), faster convergence */
        for (;;) {
            X = X - CX;
            OME = (X * XA) * (X * XA);
            FPX = PFQ / OME;
            for (J = 1; J <= JMA; J++) {
                XP[J] = DN[J] / DEK / REN[J] / OME;
                XX[J] = X * REN[J];
            }
            DIFU(2, JMA, &IERR);
            if (IERR != 0) {
                solutionIsTooHeavilyDamped = 1;
                goto end_loop_iteration;
            }
            if (cabs(D) < ADIR) {
                if ((cabs(CX) <= F32(1.E-6) * cabs(X)) /* relative frequency precision */
                    || (cabs(CX) < F32(1e-6))) {       /* absolute precision */
                    rootFindingConverged = 1;
                    if (LoopI >= 2) { /* at least 2 steps have been made */
                        goto end_loop_iteration;
                    } else {
                        goto cycle_loop_iteration;
                    }
                } else {
                    goto cycle_loop_iteration;
                }
            } else {
                X = X + CX;
                CX = CX / 2.;
                if (IRK > 3) { /* check for sitting at local minima */
                    if ((cabs(D) - ADIR) / ADIR < F32(1e-5)) {
                        /*PRINT*,' Local minima!' */
                        rootFindingConverged = 0;
                        goto end_loop_iteration;
                    }
                }
                IRK = IRK + 1;
                if (IRK > maxIterations)
                    goto end_loop_iteration;
            }
        }
    cycle_loop_iteration:;
    }
end_loop_iteration:;
}

static void allocate_output_matrices(void)
{
    /* estimate the size of matrices */
    int LoopI;
    if (kperpOUT)
        free(kperpOUT);
    if (kparOUT)
        free(kparOUT);
    if (fOUT)
        free(fOUT);
    if (ExOUT)
        free(ExOUT);
    if (EyOUT)
        free(EyOUT);
    if (EzOUT)
        free(EzOUT);
    if (BxOUT)
        free(BxOUT);
    if (ByOUT)
        free(ByOUT);
    if (BzOUT)
        free(BzOUT);
    if (SxOUT)
        free(SxOUT);
    if (SyOUT)
        free(SyOUT);
    if (SzOUT)
        free(SzOUT);
    if (EBOUT)
        free(EBOUT);
    if (VGPOUT)
        free(VGPOUT);
    if (VGZOUT)
        free(VGZOUT);
    if (SGPOUT)
        free(SGPOUT);
    if (SGZOUT)
        free(SGZOUT);
    if (uOUT)
        free(uOUT);
    if (flagSolutionFoundOUT)
        free(flagSolutionFoundOUT);
    if (flagTooHeavilyDampedOUT)
        free(flagTooHeavilyDampedOUT);
    if (flagNoConvergenceOUT)
        free(flagNoConvergenceOUT);

    if (PM[1] == PM[2]) { /* kperp is one value */
        kperpSize = 1;
        kperpOUT = (double *)calloc((size_t)kperpSize + 1, sizeof(double));
        for (LoopI = 1; LoopI <= kperpSize; LoopI++)
            kperpOUT[LoopI] = PM[1]; /* Fortran array assignment kperpOUT = PM(1) */
    } else {                     /* kperp is vector */
        kperpSize = 1 + (int)floor(((PM[1] > PM[2] ? PM[1] : PM[2]) -
                                    (PM[1] < PM[2] ? PM[1] : PM[2])) *
                                   copysign(1.0, PM[3]) / PM[3]);
        kperpOUT = (double *)calloc((size_t)kperpSize + 1, sizeof(double));
        for (LoopI = 1; LoopI <= kperpSize; LoopI++) {
            kperpOUT[LoopI] = PM[1] + (LoopI - 1.0) * PM[3];
        }
    }

    if (ZM[1] == ZM[2]) { /* scalar */
        kparSize = 1;
        kparOUT = (double *)calloc((size_t)kparSize + 1, sizeof(double));
        for (LoopI = 1; LoopI <= kparSize; LoopI++)
            kparOUT[LoopI] = ZM[1]; /* Fortran array assignment kparOUT = ZM(1) */
    } else { /* vector */
        kparSize = 1 + (int)floor(((ZM[1] > ZM[2] ? ZM[1] : ZM[2]) -
                                   (ZM[1] < ZM[2] ? ZM[1] : ZM[2])) *
                                  copysign(1.0, ZM[3]) / ZM[3]);
        kparOUT = (double *)calloc((size_t)kparSize + 1, sizeof(double));
        for (LoopI = 1; LoopI <= kparSize; LoopI++) {
            kparOUT[LoopI] = ZM[1] + (LoopI - 1.0) * ZM[3];
        }
    }

    fOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    ExOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    EyOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    EzOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    BxOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    ByOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    BzOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    SxOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    SyOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    SzOUT = (cd *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(cd));
    EBOUT = (double *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(double));
    VGPOUT = (double *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(double));
    VGZOUT = (double *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(double));
    SGPOUT = (double *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(double));
    SGZOUT = (double *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(double));
    uOUT = (double *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(double));
    flagSolutionFoundOUT =
        (int *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(int));
    flagTooHeavilyDampedOUT =
        (int *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(int));
    flagNoConvergenceOUT =
        (int *)calloc((size_t)kperpSize * (size_t)kparSize, sizeof(int));
    /* Fortran: flagSolutionFoundOUT = 0 etc. (calloc already zeroes) */
}

static void save_output(void)
{
    int indexKperp, indexKpar;
    if (PLG == PM[1]) {
        indexKperp = 1;
    } else {
        indexKperp = 1 + fnint((PLG - PM[1]) / PM[3]);
    }
    if (ZLG == ZM[1]) {
        indexKpar = 1;
    } else {
        indexKpar = 1 + fnint((ZLG - ZM[1]) / ZM[3]);
    }

    if (rootFindingConverged) {
        F2D(flagSolutionFoundOUT, indexKperp, indexKpar) = 1;
        F2D(fOUT, indexKperp, indexKpar) = X;
        F2D(ExOUT, indexKperp, indexKpar) = EFL[1];
        F2D(EyOUT, indexKperp, indexKpar) = EFL[2];
        F2D(EzOUT, indexKperp, indexKpar) = EFL[3];
        F2D(BxOUT, indexKperp, indexKpar) = BFL[1];
        F2D(ByOUT, indexKperp, indexKpar) = BFL[2];
        F2D(BzOUT, indexKperp, indexKpar) = BFL[3];
        /*     write Poynting vector uW/m^2 */
        coef_poynt = 10.0 / 4.0 / PI / 2.0;
        F2D(SxOUT, indexKperp, indexKpar) =
            creal(EFL[2] * conjg(BFL[3]) - EFL[3] * conjg(BFL[2])) * coef_poynt;
        F2D(SyOUT, indexKperp, indexKpar) =
            creal(EFL[3] * conjg(BFL[1]) - EFL[1] * conjg(BFL[3])) * coef_poynt;
        F2D(SzOUT, indexKperp, indexKpar) =
            creal(EFL[1] * conjg(BFL[2]) - EFL[2] * conjg(BFL[1])) * coef_poynt;
        F2D(EBOUT, indexKperp, indexKpar) = sqrt(
            (creal(EFL[1] * conjg(EFL[1]) + EFL[2] * conjg(EFL[2]) +
                   EFL[3] * conjg(EFL[3]))) /
            (creal(BFL[1] * conjg(BFL[1]) + BFL[2] * conjg(BFL[2]) +
                   BFL[3] * conjg(BFL[3]))));
        F2D(VGPOUT, indexKperp, indexKpar) = VG[1];
        F2D(VGZOUT, indexKperp, indexKpar) = VG[2];
        F2D(SGPOUT, indexKperp, indexKpar) = SG[1];
        F2D(SGZOUT, indexKperp, indexKpar) = SG[2];
        F2D(uOUT, indexKperp, indexKpar) = ENE;
    } else if (solutionIsTooHeavilyDamped) {
        F2D(flagTooHeavilyDampedOUT, indexKperp, indexKpar) = 1;
    } else {
        F2D(flagNoConvergenceOUT, indexKperp, indexKpar) = 1;
    }
}

void whamp_engine(void)
{
    IERR = 0;
    allocate_output_matrices();
    plasma_setup();
    if ((PM[1] == 0.0) && (ZM[1] == 0.0))
        return;
    input_setup();

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

        solutionIsTooHeavilyDamped = 0; /* default not heavily damped */
        rootFindingConverged = 0;       /* default no convergence */
        root_finding();

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
            /*if (printDebugInfo)              CALL OUTPT */
            save_output();
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
        } else {
            save_output();
        }
        if (!rootFindingConverged) {
            /*if (printDebugInfo) PRINT 125,P,Z,X,LoopI,IRK */
            if (cycleZFirst == 1)
                PLG = 1e99; /* end cycling in P */
            if (cycleZFirst == 2)
                ZLG = 1e99; /* end cycling in Z */
        }
        if (solutionIsTooHeavilyDamped) {
            /*if (printDebugInfo) PRINT*,' TOO HEAVILY DAMPED!' */
            IERR = 0;
            /*if (printDebugInfo) CALL OUTPT */
            if (cycleZFirst == 1)
                PLG = 1e99;
            if (cycleZFirst == 2)
                ZLG = 1e99;
        }
        if (cycleZFirst == 0) { /* cycle first P */
            PLG = PLG + PM[3];
            /*                   ****  UPDATE P AND Z.  **** */
            if (PLG >= PM[1] && PLG <= PM[2]) {
                P = PLG + PZL * (pow(10., PLG) - PLG);
            } else {
                ZLG = ZLG + ZM[3];
                if (ZLG < ZM[1] || ZLG > ZM[2]) {
                    break; /* exit loop_z_p */
                }
                KV = 1;
                PLG = PLO;
                P = PVO;
                Z = ZLG + PZL * (pow(10., ZLG) - ZLG);
            }
        } else if (cycleZFirst == 1) { /* cycle first Z */
            ZLG = ZLG + ZM[3];
            if (ZLG >= ZM[1] && ZLG <= ZM[2]) {
                Z = ZLG + PZL * (pow(10., ZLG) - ZLG);
            } else {
                PLG = PLG + PM[3];
                if (PLG < PM[1] || PLG > PM[2])
                    break; /* exit loop_z_p */
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
    }
    return;
}
