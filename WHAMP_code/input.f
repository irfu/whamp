      SUBROUTINE READ_INPUT_FILE(FILENAME)

      IMPLICIT REAL*8 ( A - H, O - Z ) 
      CHARACTER*(*) FILENAME
      COMMON /XPZ/ XX,PP,ZZ,A,B,D,ASS,VD,
     C   DN,TA,XP,CV,PM,ZM,XOI,XC,PZL
      COMPLEX*16 XX( 1 : 10 ), XP( 1 : 10 )
      REAL*8 PP( 1 : 10 ), ZZ( 1 : 10 ), A( 1 : 10 ), B( 1 : 10 )
      REAL*8 D( 1 : 10 )
      REAL*8 ASS( 1 : 10 ), VD( 1 : 10 ), DN( 1 : 10 ), TA( 1 : 10 )
      REAL*8 PM( 1 : 3 ), ZM( 1 : 3 )
      REAL*8 CV, XOI, XC, PZL
      CHARACTER FFILENAME*80,C*1

      DATA LU/16/
      
      write(*,*) '# input file name: '
      read(*,101)FILENAME
  101 format(A)
c      CALL CTOF_CONVERT_STRING(FILENAME,FFILENAME) 
      if (FILENAME .EQ.'1') return
      write(*,*)'# read_input_file: file = ', filename
c      OPEN(UNIT=LU,FILE=FFILENAME,ERR=500)      
      OPEN(UNIT=LU,FILE=FILENAME,ERR=500)
      GOTO 501
c  500 WRITE(*,*)'READ_INPUT_FILE: ERROR IN OPEN, FILE = ', FFILENAME
  500 WRITE(*,*)'READ_INPUT_FILE: ERROR IN OPEN, FILE = ', FILENAME
      RETURN
  501 CONTINUE

      READ(LU,*,ERR=999)(DN(I), I=1,10)
      READ(LU,*,ERR=999)(TA(I), I=1,10)
      READ(LU,*,ERR=999)(D(I), I=1,10)
      READ(LU,*,ERR=999)(A(I), I=1,10)
      READ(LU,*,ERR=999)(B(I), I=1,10)
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
