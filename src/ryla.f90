!**************************************************************
!
! File ryla.f 

      SUBROUTINE RYLA(Y,AL,RC)
      implicit none
! find the kind of a high precision variable, by finding 
! the kind of 1.0d0
integer, parameter :: d2p=kind(1.0d0)
        real(kind=d2p) :: AL,AY
      COMPLEX(kind=d2p) :: Y,RC(2,2)
!
!           ****  CHOOSE HETHOD OF EVALUATION ****
      IF(AL.LT.4.) GOTO 1
      AY=ABS(Y)
      IF(AY**2.GT.75.*AL) GOTO 1
      IF(AY.GT.40.+AL/3) GOTO 1
      IF(3.*(AL-10.).GT.AY.AND.AY**2.LT.15.*AL) GOTO 3
!          ******** NUMERICAL INTEGRATION ********
      CALL RINT(Y,AL,RC)
      RETURN
!        ******** TAYLOR SERIES ********
    1 CALL RTAY(Y,AL,RC)
      RETURN
!         ******** ASYMPTOTIC SERIES ********
    3 CALL RASY(Y,AL,RC)
      RETURN
      END 
!
