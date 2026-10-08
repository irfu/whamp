/*
 * input.c -- translation of SUBROUTINE read_input_file (src/input.f90)
 */
#include "whamp.h"
#include "comin.h"
#include "comcout.h"

void read_input_file(const char *FILENAME)
{
    int i, ioerr;
    FILE *LU = NULL;

    if (FILENAME[0] == '\0') {
        printf(" ERROR model filename not given!\n");
        exit(0);
    } else {
        if (printDebugInfo)
            printf(" # read_input_file: file =  %s\n", FILENAME);
        LU = fopen(FILENAME, "r");
        ioerr = (LU == NULL) ? 1 : 0;
        if (ioerr != 0) { /* error opening file */
            printf(" READ_INPUT_FILE: ERROR IN OPEN, FILE =  %s\n", FILENAME);
            return;
        }

        /* Fortran list-directed reads of 10 values per line */
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &DN[i]) != 1)
                goto read_error;
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &TA[i]) != 1)
                goto read_error;
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &DD[i]) != 1)
                goto read_error;
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &AA[i][1]) != 1)
                goto read_error;
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &AA[i][2]) != 1)
                goto read_error;
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &ASS[i]) != 1)
                goto read_error;
        for (i = 1; i <= 10; i++)
            if (fscanf(LU, "%lf", &VD[i]) != 1)
                goto read_error;
        if (fscanf(LU, "%lf", &XC) != 1)
            goto read_error;
        if (fscanf(LU, "%lf", &PZL) != 1)
            goto read_error;

        fclose(LU);
        return;

    read_error:
        printf(" READ_INPUT_FILE: ERROR IN READ, FILE =  %s\n", FILENAME);
        fclose(LU);
        return;
    }
}
