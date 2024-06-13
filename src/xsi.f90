!**************************************************************
!
! File xsi.f

SUBROUTINE CHI(XSI, J, IB, KOL, IERR)
   !        ARGUMENTS: XSI    CONTAINS THE SUSCEPTIBILITY TENSOR
   !                          ON RETURN.
   !                   J      COMPONENT NUMBER.
   !                   IB     INDEX FOR AA (ALPHA)
   !                   KOL    DETERMINES WHETHER DERIVATIVES
   !                          SHOULD BE EVALUATED.
   !                   IERR   ERROR FLAG, IS SET =1 IF DAMPING
   !                          IS TOO STRONG.
   use comin
   implicit none
   ! find the kind of a high precision variable, by finding
   ! the kind of 1.0d0
   integer, parameter :: d2p = kind(1.0d0)
   integer :: IB, J, KOL, IERR, L
   real(kind=d2p) :: A, ABZ, AL, ALA, P, RX, VZ, Z, XI, XR
   COMPLEX(kind=d2p) X, XY, AI, BL, CL, DL, BLY, RC(2, 2), XSI(6, 4),&
        & PS, PSP, PSY, PPY, DP, DPP, Y, ZY

   !          **** RESIDUES FOR PADE APPROXIMANT ****
   complex(kind=d2p), dimension(8) :: &
        & B = [(-1.734012457471826d-2, -4.630639291680322d-2),&
        &      (-1.734012457471826d-2,  4.630639291680322d-2),&
        &      (-7.399169923225014d-1,  8.395179978099844d-1),&
        &      (-7.399169923225014d-1, -8.395179978099844d-1),&
        &      ( 5.840628642184073d0,   9.536009057643667d-1),&
        &      ( 5.840628642184073d0,  -9.536009057643667d-1),&
        &      (-5.583371525286853d0,  -1.120854319126599d1),&
        &      (-5.583371525286853d0,   1.120854319126599d1)],&
        & C = [( 2.237687789201900d0,  -1.625940856173727d0),&
        &      (-2.237687789201900d0,  -1.625940856173727d0),&
        &      ( 1.465234126106004d0,  -1.789620129162444d0),&
        &      (-1.465234126106004d0,  -1.789620129162444d0),&
        &      (  .8392539817232638d0, -1.891995045765206d0),&
        &      ( -.8392539817232638d0, -1.891995045765206d0),&
        &      (  .2739362226285564d0, -1.941786875844713d0),&
        &      ( -.2739362226285564d0, -1.941786875844713d0)]
   !
   X = XX(J)
   Z = ZZ(J)
   P = PP(J)
   A = AA(J, IB)
   VZ = VD(J)*Z
   AI = (0.d0, 1.d0)*A
   IF (ASS(J) .EQ. 0.) AI = -AI
   AL = .5*P*P
   ALA = A*AL
   !
   XSI(1:6, 1:4) = (0.d0, 0.d0)

   !                         TEST FOR STRONG DAMPING
   !  Some tests modified to allow Z < 0. (1996-02-20, Kjell R.)
   XI = DIMAG(X)
   XR = ABS(REAL(X) - Z*VD(J))
   ABZ = ABS(Z)
   do
      IF (XI .GE. -ABZ) exit
      !  Allow high frequency, cold plasma waves ...(ADDED 1993-04-02, KJELL R.)
      IF (XR .GT. 1.+5.*ABZ .AND. XR .GT. 1.+5.*P) exit
      RX = XR - INT(XR)
      IF (RX .GT. 0.5) RX = 1.-RX
      IF (XI .GE. -0.6*RX) exit
      IERR = 1
      RETURN
   end do

   XSI(1, 1) = A
   XSI(6, 1) = A + 2*A*VD(J)**2
   loop: DO L = 1, 8
      BL = B(L)
      !  Modified to allow Z < 0. (1996-02-20, Kjell R.)
      !**** CL=C(L)
      if (z .lt. 0.) then
         CL = -C(L)
      else
         CL = C(L)
      end if
      DL = CL + VD(J)
      Y = X - DL*Z
      BLY = BL/Y
      !        ****** EVALUATE THE R-FUNCTION ******
      CALL RYLA(Y, ALA, RC)
      XY = 1.+A*Z*CL/Y
      PS = XY*RC(1, 1)
      PSP = XY*RC(1, 2)
      !                     ****** FORM SUCEPTIBILITY TENSOR ******
      XSI(1, 1) = XSI(1, 1) + A*BL*Y*PS
      XSI(2, 1) = XSI(2, 1) + AI*BL*PSP
      XSI(3, 1) = XSI(3, 1) + A*P*BL*DL*PS
      XSI(4, 1) = XSI(4, 1) + BLY*PSP
      !CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
      !****      XSI(5,1)=XSI(5,1)+AI*P*BLY*DL*PSP
      XSI(5, 1) = XSI(5, 1) - AI*P*BLY*DL*PSP
      XSI(6, 1) = XSI(6, 1) + 2.*A*BLY*DL*DL*(X - VZ + AL*PS)
      !
      IF (KOL .LE. 1) cycle loop
      !                    ****** FORM X-DERIVATIVES OF XSI ******
      PSY = XY*RC(2, 1)
      PPY = XY*RC(2, 2)
      XY = X/Y
      DP = XY*(PSY - A*Z*CL/Y*RC(1, 1))
      DPP = XY*(PPY - A*Z*CL/Y*RC(1, 2))
      XSI(1, 2) = XSI(1, 2) + A*BL*(Y*DP + X*PS)
      XSI(2, 2) = XSI(2, 2) + AI*BL*DPP
      XSI(3, 2) = XSI(3, 2) + A*P*BL*DL*DP
      XSI(4, 2) = XSI(4, 2) + BLY*(DPP - XY*PSP)
      !CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
      !****       XSI(5,2)=XSI(5,2)+AI*P*BLY*DL*(DPP-XY*PSP)
      XSI(5, 2) = XSI(5, 2) - AI*P*BLY*DL*(DPP - XY*PSP)
      XSI(6, 2) = XSI(6, 2) + 2.*A*BLY*DL*DL*(X + AL*DP - XY*(X - VZ + AL*PS))
      !
      IF (KOL .LE. 2) cycle loop
      !                    ****** FORM Z-DERIVATIVES OF XSI ******
      ZY = Z/Y
      DP = ZY*(A*CL*XY*RC(1, 1) - DL*PSY)
      !CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
      !****       XSI(5,2)=XSI(5,2)+AI*P*BLY*DL*(DPP-XY*PSP)
      XSI(5, 2) = XSI(5, 2) - AI*P*BLY*DL*(DPP - XY*PSP)

      DPP = ZY*(A*CL*XY*RC(1, 2) - DL*PPY)
      ZY = DL*ZY
      XSI(1, 3) = XSI(1, 3) + A*BL*Y*(DP - ZY*PS)
      XSI(2, 3) = XSI(2, 3) + AI*BL*DPP
      XSI(3, 3) = XSI(3, 3) + A*P*BL*DL*DP
      XSI(4, 3) = XSI(4, 3) + BLY*(DPP + ZY*PSP)
      !CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
      !****       XSI(5,3)=XSI(5,3)+AI*P*BLY*DL*(DPP+ZY*PSP)
      XSI(5, 3) = XSI(5, 3) - AI*P*BLY*DL*(DPP + ZY*PSP)
      XSI(6, 3) = XSI(6, 3) + 2.*A*BLY*DL*DL*(AL*DP - VZ + ZY*(X - VZ + AL*PS))
      !
      IF (KOL .LE. 3) cycle loop
      !                   ****** FORM P-DERIVATIVES OF XSI ******
      CALL RYLA(Y - 1., ALA, RC)
      DP = 2.*(PSP - PS)
      !      DPP=2.*AL*((Y/(Y-1.))**2*RC(1,2)-PSP)-Y*DP
      !CCCCCCorrection, 1991 - 03 - 13, Kjell R
      !      DPP=2.*AL*((Y/(Y-1.))**2*(1.+A*Z*CL/Y)*RC(1,2)-PSP)-Y*DP
      !     As pointed out by Scott Boardsen, NASA/MSFC, the AL above should
      !     really be ALA.   Corrected  1991 - 12 - 18, Kjell R
      DPP = 2.*ALA*((Y/(Y - 1.))**2*(1.+A*Z*CL/Y)*RC(1, 2) - PSP) - Y*DP
      XSI(1, 4) = XSI(1, 4) + A*BL*Y*DP
      XSI(2, 4) = XSI(2, 4) + AI*BL*DPP
      XSI(3, 4) = XSI(3, 4) + A*P*BL*DL*(DP + PS)
      XSI(4, 4) = XSI(4, 4) + BLY*(2.*PSP + DPP)
      !CCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
      !****       XSI(5,4)=XSI(5,4)+AI*P*BLY*DL*(PSP+DPP)
      XSI(5, 4) = XSI(5, 4) - AI*P*BLY*DL*(PSP + DPP)
      XSI(6, 4) = XSI(6, 4) + 4.*A*BLY*DL*DL*AL*PSP
   end DO loop
   !      **** COMPLETE XSI(4, ) ****
   XSI(4, 1) = XSI(1, 1) - 2.*A**2*AL*XSI(4, 1)
   XSI(4, 2) = XSI(1, 2) - 2.*A**2*AL*XSI(4, 2)
   XSI(4, 3) = XSI(1, 3) - 2.*A**2*AL*XSI(4, 3)
   XSI(4, 4) = XSI(1, 4) - 2.*A**2*AL*XSI(4, 4)
end SUBROUTINE CHI
