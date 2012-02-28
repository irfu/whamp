!**************************************************************
!
! File comcout.h 

module comcout
implicit none
! find the kind of a high precision variable, by finding 
! the kind of 1.0d0
integer, parameter,private :: d2p=kind(1.0d0)

complex(kind=d2p) :: X, EFL(3), BFL(3), D ,DX, DZ, DP, E(6,4)
      
real(kind=d2p) :: P,Z, VG(2),SG(2),RI(2),ENE

end module comcout
