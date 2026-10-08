/*
 * ctof.c -- translation of SUBROUTINE CTOF_CONVERT_STRING (src/ctof.f90)
 */
#include "whamp.h"

void CTOF_CONVERT_STRING(const char *CSTRING, char *FSTRING)
{
    enum
    {
        MAX_LENGTH = 80
    };
    int i;
    char C = ' ';

    /*
     *       CALL CEOS(C)
     */
    for (i = 1; i <= MAX_LENGTH; i++) {
        if (CSTRING[i - 1] == C) {
            /* Fortran: FSTRING = CSTRING(1:I-1) */
            memcpy(FSTRING, CSTRING, (size_t)(i - 1));
            FSTRING[i - 1] = '\0';
            return;
        }
    }
}
