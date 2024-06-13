!**************************************************************
!
! File difu.f

SUBROUTINE DIFU(KOL, JMAX, IERR)
   !        ARGUMENTS: KOL     =2 D AND ITS X-DERIVATIVE COMPUTED.
   !                           =4 ALL DERIVATIVES AND WAVE FIELDS.
   !                   JMAX    =1 - 6, NUMBER OF COMPONENTS.
   !                   IERR    ERROR FLAG
   use comin
   use comcout
   implicit none
   integer, parameter :: d2p = kind(1.0d0)

   integer :: IB, IERR, J, JMAX, KOL, KOLLA
   real(kind=d2p) :: Q
   COMPLEX(kind=d2p) :: XSI(6, 4), DF, U1, U2, U3, U12, U13, U32,&
        & A, B, C, DA, DB, DC

   !              *********** FORM DIELECTRIC TENSOR ************

   E(1:6, 1:4) = (0.d0, 0.d0)
   E(1, 1) = 1.
   E(4, 1) = 1.
   E(6, 1) = 1.
   !
   JMA = JMAX
   KOLLA = MIN(KOL, 4)
   species_loop: DO J = 1, JMA
      E(1, 1) = E(1, 1) - XP(J)
      E(4, 1) = E(4, 1) - XP(J)
      E(6, 1) = E(6, 1) - XP(J)
      IF (AA(J, 1) .EQ. AA(J, 2)) then
         AA(J, 2) = 0.
         DD(J) = 1.
      end IF
      IB = 1
      DF = XP(J)/(AA(J, 1)*(AA(J, 1) - AA(J, 2)))
      Q = AA(J, 1) - DD(J)*AA(J, 2)
      chi_loop: do
         CALL CHI(XSI, J, IB, KOL, IERR)
         IF (IERR .NE. 0) RETURN
         E(1:6, 1:KOLLA) = E(1:6, 1:KOLLA) + DF*Q*XSI(1:6, 1:KOLLA)
         !
         IF (IB .EQ. 2) cycle species_loop
         Q = (DD(J) - 1.)*AA(J, 1)
         IF (Q .EQ. 0.) cycle species_loop
         IB = 2
      end do chi_loop
   end DO species_loop
   !                       *** DIELECTRIC TENSOR COMPUTED ***
   !
   !       ******* FORM REFRACTIVE INDEX, CV=SPEED OF LIGHT/THERM. SPEED. *
   U1 = PP(1)*CV/XX(1)
   U3 = ZZ(1)*CV/XX(1)
   U12 = U1*U1
   U32 = U3*U3
   U2 = U12 + U32
   U13 = 2.*U1*U3
   !
   !                      ******** FORM DISPERSION FUNCTION ********
   A = U12*E(1, 1) + U13*E(3, 1) + U32*E(6, 1)
   !*************  Sign error corrected in Feb. 1989.   Kjell R.  ******
   !***  B=U2*(E(1,1)*E(6,1)-E(3,1)**2)+(U3*E(5,1)+U1*E(2,1))**2
   B = U2*(E(1, 1)*E(6, 1) - E(3, 1)**2) + (U3*E(5, 1) - U1*E(2, 1))**2
   C = (E(1, 1)*E(6, 1) - E(3, 1)**2)*E(4, 1) + E(6, 1)*E(2, 1)**2
   C = C + (E(1, 1)*E(5, 1) + 2.*E(2, 1)*E(3, 1))*E(5, 1)
   !
   D = (U2 - E(4, 1))*A - B + C
   IF (KOL .LE. 1) RETURN
   !        ****** COMPLETE X-DERIVATIVE OF DIELECTRIC TENSOR ******
   E(1, 2) = E(1, 2) - 2.*(E(1, 1) - 1.)
   E(2, 2) = E(2, 2) - 2.*E(2, 1)
   E(3, 2) = E(3, 2) - 2.*E(3, 1)
   E(4, 2) = E(4, 2) - 2.*(E(4, 1) - 1.)
   E(5, 2) = E(5, 2) - 2.*E(5, 1)
   E(6, 2) = E(6, 2) - 2.*(E(6, 1) - 1.)
   !          ****** X-DERIVATIVE OF DISPERSION FUNCTION *******
   DA = (E(1, 2) - 2.*E(1, 1))*U12 + U13*(E(3, 2) - 2.*E(3, 1)) +&
        & (E(6, 2) - 2.*E(6, 1))*U32
   DB = 2.*(U3*E(5, 1) - U1*E(2, 1))*(U3*(E(5, 2) - E(5, 1)) - U1*(E(2, 2) - E(2, 1)&
        &)) + U2*((E(1, 2) - E(1, 1))*E(6, 1) - 2.*(E(3, 2) - E(3, 1))*E(3, 1) +&
        &(E(6, 2) - E(6, 1))*E(1, 1))
   DC = (E(1, 2)*E(6, 1) + E(1, 1)*E(6, 2) - 2.*E(3, 1)*E(3, 2))*E(4, 1)
   DC = DC + (E(1, 1)*E(6, 1) - E(3, 1)**2)*E(4, 2)
   DC = DC + E(5, 1)*(E(1, 2)*E(5, 1) + 2.*E(1, 1)*E(5, 2) +&
        &E(2, 1)*E(3, 2) + 2.*E(2, 2)*E(3, 1))
   DC = DC + E(2, 1)*(E(6, 2)*E(2, 1) + 2.*E(6, 1)*E(2, 2) +&
        &E(5, 1)*E(3, 2) + 2.*E(5, 2)*E(3, 1))
   !
   DX = ((U2 - E(4, 1))*DA - (2.*U2 + E(4, 2))*A - DB + DC)/XX(1)
   IF (KOL .LE. 2) RETURN
   DZ = (0., 0.)
   IF (ZZ(1) .NE. 0.) then
      !         ****** Z-DERIVATIVE OF DISPERSION FUNCTION ******
      DA = U12*E(1, 3) + U13*(E(3, 3) + E(3, 1)) + U32*(E(6, 3) + 2.*E(6, 1))
      !      DB=2.*(U3*E(5,1)-U1*E(3,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,1))+
      !CCCCCCorrection 1991- 03 - 11, Kjell R.
      !      DB=2.*(U3*E(5,1)-U1*E(2,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,1))+
      !CCCCCCorrection of error found by Scott Boardsen, NASA/MSFC. 1992-03-12, Kjell R.
      DB = 2.*(U3*E(5, 1) - U1*E(2, 1))*(U3*(E(5, 3) + E(5, 1)) - U1*E(2, 3)) + &
           &2.*U32*(E(1, 1)*E(6, 1) - E(3, 1)**2) +&
           &U2*(E(1, 3)*E(6, 1) + E(1, 1)*E(6, 3) - 2.*E(3, 1)*E(3, 3))
      DC = (E(1, 3)*E(6, 1) + E(1, 1)*E(6, 3) - 2.*E(3, 1)*E(3, 3))*E(4, 1)
      DC = DC + (E(1, 1)*E(6, 1) - E(3, 1)**2)*E(4, 3)
      DC = DC + E(5, 1)*(E(1, 3)*E(5, 1) + 2.*E(1, 1)*E(5, 3) +&
           &E(2, 1)*E(3, 3) + 2.*E(2, 3)*E(3, 1))
      DC = DC + E(2, 1)*(E(6, 3)*E(2, 1) + 2.*E(6, 1)*E(2, 3) +&
           &E(5, 1)*E(3, 3) + 2.*E(5, 3)*E(3, 1))
      !
      DZ = ((U2 - E(4, 1))*DA + (2.*U32 - E(4, 3))*A - DB + DC)/ZZ(1)
   end IF
   IF (KOL .LE. 3) RETURN
   !        ****** P-DERIVATIVE OF DISPERSION FUNCTION ******
   DP = (0., 0.)
   IF (PP(1) .NE. 0.) then
      !
      DA = U12*(E(1, 4) + 2.*E(1, 1)) + U13*(E(3, 4) + E(3, 1)) + U32*E(6, 4)
      DB = 2.*(U3*E(5, 1) - U1*E(2, 1))*(U3*E(5, 4) - U1*(E(2, 4) + E(2, 1))) +&
           &2.*U12*(E(1, 1)*E(6, 1) - E(3, 1)**2) +&
           &U2*(E(1, 4)*E(6, 1) + E(1, 1)*E(6, 4) - 2.*E(3, 1)*E(3, 4))
      DC = (E(1, 4)*E(6, 1) + E(1, 1)*E(6, 4) - 2.*E(3, 1)*E(3, 4))*E(4, 1)
      DC = DC + (E(1, 1)*E(6, 1) - E(3, 1)**2)*E(4, 4)
      DC = DC + E(5, 1)*(E(1, 4)*E(5, 1) + 2.*E(1, 1)*E(5, 4) +&
           &E(2, 1)*E(3, 4) + 2.*E(2, 4)*E(3, 1))
      DC = DC + E(2, 1)*(E(6, 4)*E(2, 1) + 2.*E(6, 1)*E(2, 4) +&
           &E(5, 1)*E(3, 4) + 2.*E(5, 4)*E(3, 1))
      !
      DP = ((U2 - E(4, 1))*DA + (2.*U12 - E(4, 4))*A - DB + DC)/PP(1)
   end IF
   !                     ******** COMPUTE ELECTRIC FIELD ********
   CALL ENERGY(U1, U3, U2, U12, U32)
   !       THE ELECTRIC FIELD IS 1 MV/M.
   !       THE MAGNETIC FIELD WILL BE IN GAMMA
END SUBROUTINE DIFU

