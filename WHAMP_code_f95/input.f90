      SUBROUTINE READ_INPUT_FILE(FILENAME)
      use comin
      use comcout
      implicit none
        integer :: I,LU=16
      CHARACTER*(*) FILENAME
      CHARACTER FFILENAME*80

      
      write(*,*) '# input file name: '
      read(*,101)FILENAME
  101 format(A)
!      CALL CTOF_CONVERT_STRING(FILENAME,FFILENAME) 
      if (FILENAME .EQ.'1') return
      write(*,*)'# read_input_file: file = ', filename
!      OPEN(UNIT=LU,FILE=FFILENAME,ERR=500)      
      OPEN(UNIT=LU,FILE=FILENAME,ERR=500)
      GOTO 501
!  500 WRITE(*,*)'READ_INPUT_FILE: ERROR IN OPEN, FILE = ', FFILENAME
  500 WRITE(*,*)'READ_INPUT_FILE: ERROR IN OPEN, FILE = ', FILENAME
      RETURN
  501 CONTINUE

      READ(LU,*,ERR=999)(DN(I), I=1,10)
      READ(LU,*,ERR=999)(TA(I), I=1,10)
      READ(LU,*,ERR=999)(DD(I), I=1,10)
      READ(LU,*,ERR=999)(AA(I,1), I=1,10)
      READ(LU,*,ERR=999)(AA(I,2), I=1,10)
      READ(LU,*,ERR=999)(ASS(I), I=1,10)
      READ(LU,*,ERR=999)(VD(I), I=1,10)
      READ(LU,*)XC
      READ(LU,*)PZL

      CLOSE(LU)
      RETURN

  999 WRITE(*,*)'READ_INPUT_FILE: ERROR IN READ, FILE = ',FFILENAME
      CLOSE(LU)
      RETURN
      END
