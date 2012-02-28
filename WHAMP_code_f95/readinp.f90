subroutine  READ_INPUT_FILE(FILENAME)
  use comin
  implicit none
  integer :: i
  CHARACTER FILENAME*(80)

  OPEN (3, FILE = FILENAME)  
  DO i = 1,10
     READ (3,*) DN(i)
  end do
  DO i = 1,10
     READ (3,*) TA(i)
  end do
END subroutine READ_INPUT_FILE
