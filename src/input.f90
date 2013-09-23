      SUBROUTINE READ_INPUT_FILE(FILENAME)
      use comin
      use comcout
      implicit none
      integer :: I,ioerr,LU=16
      CHARACTER*(*) FILENAME
      CHARACTER FFILENAME*80
      
      write(*,*) '# input file name: '
      read(*,101)FILENAME
  101 format(A)
      if (FILENAME .EQ.'1') return
      write(*,*)'# read_input_file: file = ', filename
      OPEN(UNIT=LU,FILE=FILENAME,iostat=ioerr)
      if (ioerr /= 0) then ! error opening file
          WRITE(*,*)'READ_INPUT_FILE: ERROR IN OPEN, FILE = ', FILENAME
          RETURN
      end if

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
