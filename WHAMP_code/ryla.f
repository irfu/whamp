***************************************************************
*
* File ryla.f 

      SUBROUTINE RYLA(Y,AL,RC)
      IMPLICIT REAL*8 ( A - H, O - Z )
      COMPLEX*16 Y,RC(2,2)
C
C           ****  CHOOSE HETHOD OF EVALUATION ****
      IF(AL.LT.4.) GOTO 1
      AY=ABS(Y)
      IF(AY**2.GT.75.*AL) GOTO 1
      IF(AY.GT.40.+AL/3) GOTO 1
      IF(3.*(AL-10.).GT.AY.AND.AY**2.LT.15.*AL) GOTO 3
C          ******** NUMERICAL INTEGRATION ********
      CALL RINT(Y,AL,RC)
      RETURN
C        ******** TAYLOR SERIES ********
    1 CALL RTAY(Y,AL,RC)
      RETURN
C         ******** ASYMPTOTIC SERIES ********
    3 CALL RASY(Y,AL,RC)
      RETURN
      END 
*
