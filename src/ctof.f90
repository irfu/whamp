!
!***********************************************************************
! 
SUBROUTINE CTOF_CONVERT_STRING(CSTRING,FSTRING)
  ! 
  implicit none
  integer,parameter :: MAX_LENGTH = 80
  integer :: I
  CHARACTER C*1 
  CHARACTER*(*) CSTRING,FSTRING
  !
  !       CALL CEOS(C)
  C=' '
  DO I=1,MAX_LENGTH
     IF(CSTRING(I:I).EQ.C)THEN
        FSTRING=CSTRING(1:I-1)
        RETURN
     ENDIF

  END DO
end SUBROUTINE CTOF_CONVERT_STRING
