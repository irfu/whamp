!**************************************************************
!
! File rtay.f

SUBROUTINE RTAY(Y, AL, RC)
   implicit none
   ! find the kind of a high precision variable, by finding the kind of 1.0d0
   integer, parameter :: d2p = kind(1.0d0)
   integer :: I
   real(kind=d2p) :: AL, T
   COMPLEX(kind=d2p) :: Y, Y2, RC(2, 2), PN, PYN, COT
   !                  ******** TAYLOR SERIES ********
   Y2 = Y*Y
   PN = Y/(Y2 - 1.)
   PYN = -Y*(Y2 + 1.)/(Y2 - 1.)**2
   RC(1, 1:2) = PN
   RC(2, 1:2) = PYN
   !
   DO I = 2, 100
      COT = (2*I - 1)/(Y2 - I**2)*AL
      PYN = COT*(PYN - 2.*Y2/(Y2 - I**2)*PN)
      PN = COT*PN
      RC(1, 1) = RC(1, 1) + PN
      RC(2, 1) = RC(2, 1) + PYN
      RC(1, 2) = RC(1, 2) + I*PN
      RC(2, 2) = RC(2, 2) + I*PYN
      T = ABS(PN)*1.E8
      IF (T .LT. ABS(RC(1, 1))) exit
   end DO
end SUBROUTINE RTAY
