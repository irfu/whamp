! Select method of evaluation
SUBROUTINE RYLA(Y, AL, RC)
   implicit none
   integer, parameter :: d2p = kind(1.0d0)
   real(kind=d2p) :: AL, AY
   COMPLEX(kind=d2p) :: Y, RC(2, 2)
   !
   !           ****  CHOOSE HETHOD OF EVALUATION ****
   IF (AL .LT. 4.) then
      CALL RTAY(Y, AL, RC) ! Taylor series
      RETURN
   end if
   AY = ABS(Y)
   IF (AY**2 .GT. 75.*AL) then
      CALL RTAY(Y, AL, RC) ! Taylor series
      RETURN
   end if
   IF (AY .GT. 40.+AL/3) then
      CALL RTAY(Y, AL, RC) ! Taylor series
      RETURN
   end if

   IF (3.*(AL - 10.) .GT. AY .AND. AY**2 .LT. 15.*AL) then
      CALL RASY(Y, AL, RC) ! asymptotic series
      RETURN
   else
      CALL RINT(Y, AL, RC) ! numerical integration
      RETURN
   end if
END SUBROUTINE
