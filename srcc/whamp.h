/*
 * whamp.h -- shared definitions for the C translation of the WHAMP
 * Fortran sources in ../src.
 *
 * Conventions used throughout the translation:
 *   - Fortran `real(kind(1.0d0))` (d2p) -> C `double`
 *   - Fortran `complex(kind(1.0d0))`    -> C `double complex` (typedef cd)
 *   - Fortran `real` (default)          -> C `float`
 *   - Fortran arrays keep 1-based indexing: an array declared
 *     A(10) in Fortran becomes A[11] in C and is used with indices 1..10.
 *   - Fortran module variables become globals defined in comin.c,
 *     comcout.c and comoutput.c.
 */
#ifndef WHAMP_H
#define WHAMP_H

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double complex cd;

/* Fortran unsuffixed real literals (e.g. 1836.1, 0.01, 1e-6) are SINGLE
 * precision.  When such a literal participates in a double-precision
 * expression it is first rounded to float and then widened.  Wrap every
 * unsuffixed literal that is not exactly representable in float with F32()
 * so the C matches gfortran bit-for-bit. */
#define F32(x) ((double)(float)(x))

/* Fortran intrinsic helpers ------------------------------------------- */

/* Fortran INT(x): truncation toward zero */
static inline double fint(double x) { return trunc(x); }
/* Fortran NINT(x): nearest integer, ties away from zero */
static inline int fnint(double x)
{
    return (int)(x >= 0.0 ? floor(x + 0.5) : ceil(x - 0.5));
}
/* Fortran SIGN(A,B): magnitude of A, sign of B */
#define FSIGN(A, B) copysign((A), (B))

/* Fortran CONJG() */
#define conjg(z) conj(z)

/* Emit an edited value right-justified in a field of width w, applying
 * gfortran's fallbacks: when the value does not fit, the leading zero of
 * "0."/"-0." is dropped; if it still does not fit the field is filled
 * with '*'. */
static inline void femit(int w, const char *s)
{
    int len = (int)strlen(s);
    if (len <= w) {
        printf("%*s", w, s);
        return;
    }
    {
        const char *p = (*s == '-') ? s + 1 : s;
        if (p[0] == '0' && p[1] == '.') {
            char alt[512];
            int n = 0;
            if (*s == '-')
                alt[n++] = '-';
            strcpy(alt + n, p + 1);
            if ((int)strlen(alt) <= w) {
                printf("%*s", w, alt);
                return;
            }
        }
    }
    {
        int i;
        for (i = 0; i < w; i++)
            putchar('*');
    }
}

/* Emulate Fortran Ew.d editing with P scaling of 0 or 1.
 * p=0 prints 0.d1d2...d E±ee (Fortran default), p=1 prints d.d1...d E±ee.
 * gfortran leaves out the 'E' when the exponent needs three digits.
 * Rounding is done by the C library on the exact binary value, which is
 * what gfortran does as well. */
static inline void fE(int w, int d, int p, double v)
{
    char tmp[512], out[512];
    int expo, neg = 0;
    char *ep;

    /* gfortran edits IEEE special values as fixed strings right-justified
     * in the field: NaN, Infinity, -Infinity. */
    if (isnan(v)) {
        femit(w, "NaN");
        return;
    }
    if (isinf(v)) {
        femit(w, v < 0.0 ? "-Infinity" : "Infinity");
        return;
    }
    if (v < 0.0 || (v == 0.0 && signbit(v))) {
        neg = 1;
        v = -v;
    }
    if (v == 0.0) {
        int i, o = 0;
        if (neg)
            out[o++] = '-';
        out[o++] = '0';
        out[o++] = '.';
        for (i = 0; i < d; i++)
            out[o++] = '0';
        o += snprintf(out + o, sizeof out - (size_t)o, "E+00");
        (void)o;
        femit(w, out);
        return;
    }
    if (p == 1) {
        int ex, o = 0;
        snprintf(tmp, sizeof tmp, "%.*e", d, v);
        ep = strchr(tmp, 'e');
        if (!ep)
            ep = strchr(tmp, 'E');
        ex = atoi(ep + 1);
        if (neg)
            out[o++] = '-';
        memcpy(out + o, tmp, (size_t)(ep - tmp));
        o += (int)(ep - tmp);
        o += snprintf(out + o, sizeof out - (size_t)o, "%s%c%02d",
                      (ex > 99 || ex < -99) ? "" : "E", ex < 0 ? '-' : '+',
                      ex < 0 ? -ex : ex);
        out[o] = '\0';
        femit(w, out);
        return;
    }
    /* p == 0: exactly d significant digits, printed as 0.d1..d E±ee */
    snprintf(tmp, sizeof tmp, "%.*e", d > 0 ? d - 1 : 0, v);
    ep = strchr(tmp, 'e');
    if (!ep)
        ep = strchr(tmp, 'E');
    expo = atoi(ep + 1) + 1;
    {
        char digits[512];
        int nd = 0, i, o = 0;
        const char *q = tmp;
        for (; q != ep; q++)
            if (*q >= '0' && *q <= '9')
                digits[nd++] = *q;
        if (neg)
            out[o++] = '-';
        out[o++] = '0';
        out[o++] = '.';
        for (i = 0; i < d; i++)
            out[o++] = digits[i < nd ? i : nd - 1];
        snprintf(out + o, sizeof out - (size_t)o, "%s%c%02d",
                 (expo > 99 || expo < -99) ? "" : "E", expo < 0 ? '-' : '+',
                 expo < 0 ? -expo : expo);
        femit(w, out);
    }
}

/* Emulate Fortran Fw.d editing (round-half-even on the binary value). */
static inline void fF(int w, int d, double v)
{
    char out[512];
    /* gfortran edits IEEE special values as fixed strings right-justified
     * in the field: NaN, Infinity, -Infinity. */
    if (isnan(v)) {
        femit(w, "NaN");
        return;
    }
    if (isinf(v)) {
        femit(w, v < 0.0 ? "-Infinity" : "Infinity");
        return;
    }
    snprintf(out, sizeof out, "%.*f", d, v);
    femit(w, out);
}

/* Emulate Fortran list-directed editing of a real(8) item.
 * gfortran writes a 26 character field:
 *   - fixed notation when the adjusted exponent E
 *     (value = 0.d1d2... x 10^E) satisfies 0 <= E <= 17, with 17
 *     significant digits, right justified in 21 characters and the field
 *     padded to 26 with trailing blanks;
 *   - otherwise exponential P=1 notation d.dddddE+ee (3 digit exponent),
 *     right justified in the full 26 characters. */
static inline void ld_real_w(int w, double v)
{
    char out[64];
    int expo;
    if (v == 0.0 || !isfinite(v)) {
        snprintf(out, sizeof out, "%.*f", 16, v);
        printf("%*s%5s", w - 5, out, "");
        return;
    }
    expo = (int)floor(log10(fabs(v))) + 1;
    if (expo >= 0 && expo <= 17) {
        char fixed[64];
        int decimals = 17 - expo;
        if (decimals > 0)
            snprintf(fixed, sizeof fixed, "%.*f", decimals, v);
        else
            snprintf(fixed, sizeof fixed, "%.0f.", v); /* Fortran keeps the point */
        printf("%*s%5s", w - 5, fixed, "");
    } else {
        char tmp[64], mant[64];
        int ex;
        snprintf(tmp, sizeof tmp, "%.16e", v);
        {
            char *ep = strchr(tmp, 'e');
            if (!ep)
                ep = strchr(tmp, 'E');
            ex = atoi(ep + 1);
            *ep = '\0';
            snprintf(mant, sizeof mant, "%s", tmp);
        }
        snprintf(out, sizeof out, "%sE%c%03d", mant, ex < 0 ? '-' : '+', ex < 0 ? -ex : ex);
        printf("%*s", w, out);
    }
}
/* Fortran list-directed write of one real(8) item (field width 26). */
static inline void ld_real(double v) { ld_real_w(26, v); }
/* Fortran list-directed write of one integer item (field width 12). */
static inline void ld_int(int v) { printf("%12d", v); }
/* Fortran list-directed write of one logical item. */
static inline void ld_logical(int v) { printf("%12s", v ? "T" : "F"); }
/* Fortran list-directed write of one complex(8) item. */
static inline void ld_complex(cd z)
{
    printf("(");
    ld_real_w(25, creal(z));
    printf(",");
    ld_real_w(25, cimag(z));
    printf(")");
}

/* Prototypes ------------------------------------------------------------ */

/* energy.c */
void ENERGY(cd U1, cd U3, cd U2, cd U12, cd U32);
/* rint.c */
void RINT(cd YY, double AL, cd *RC);
/* rtay.c */
void RTAY(cd Y, double AL, cd *RC);
/* rasy.c */
void RASY(cd Y, double AL, cd *RC);
/* ryla.c */
void RYLA(cd Y, double AL, cd *RC);
/* xsi.c */
void CHI(cd *XSI, int J, int IB, int KOL, int *IERR);
/* difu.c */
void DIFU(int KOL, int JMAX, int *IERR);
/* av.c */
void AV(void);
/* ctof.c */
void CTOF_CONVERT_STRING(const char *CSTRING, char *FSTRING);
/* input.c */
void read_input_file(const char *FILENAME);
/* output.f90 */
void OUTPT(void);
void INOUT(void);
/* typin.c */
void TYPIN(int *NPL, int *KFS);
/* whamp_engine.c */
void whamp_engine(void);

#endif /* WHAMP_H */
