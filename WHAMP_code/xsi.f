***************************************************************
*
* File xsi.f 

      SUBROUTINE CHI(XSI,J,IB,KOL,IERR)
C        ARGUMENTS: XSI    CONTAINS THE SUSCEPTIBILITY TENSOR
C                          ON RETURN.
C                   J      COMPONENT NUMBER.
C                   IB     INDEX FOR AA (ALPHA)
C                   KOL    DETERMINES WHETHER DERIVATIVES
C                          SHOULD BE EVALUATED.
C                   IERR   ERROR FLAG, IS SET =1 IF DAMPING
C                          IS TOO STRONG.

      include 'comin.h'
      COMPLEX*16 X,XY,AI,BL,CL,DL,BLY,RC(2,2),XSI(6,4),
     2  B(8),C(8),PS,PSP,PSY,PPY,DP,DPP,Y,ZY

C          **** RESIDUES FOR PADE APPROXIMANT ****
      DATA B/(-1.734012457471826E-2,-4.630639291680322E-2),
     B       (-1.734012457471826E-2, 4.630639291680322E-2),
     B       (-7.399169923225014E-1, 8.395179978099844E-1),
     B       (-7.399169923225014E-1,-8.395179978099844E-1),
     B       (5.840628642184073    , 9.536009057643667E-1),
     B       (5.840628642184073    ,-9.536009057643667E-1),
     B       (-5.583371525286853   ,-1.120854319126599 E1),
     B       (-5.583371525286853   , 1.120854319126599 E1)/,
     C    C/ ( 2.237687789201900,  -1.625940856173727),
     C       (-2.237687789201900,  -1.625940856173727),
     C       ( 1.465234126106004,  -1.789620129162444),
     C       (-1.465234126106004,  -1.789620129162444),
     C       ( .8392539817232638,  -1.891995045765206),
     C       (-.8392539817232638,  -1.891995045765206),
     C       ( .2739362226285564,  -1.941786875844713),
     C       (-.2739362226285564,  -1.941786875844713)/
C
      X=XX(J)
      Z=ZZ(J)
      P=PP(J)
      A=AA(J,IB)
      VZ=VD(J)*Z
      AI=(0.d0,1.d0)*A
      IF(ASS(J).EQ.0.) AI=-AI
      AL=.5*P*P
      ALA=A*AL
C
      CALL ZEROC2( XSI, 1, 6, 1, 4 )

C                         TEST FOR STRONG DAMPING
c  Some tests modified to allow Z < 0. (1996-02-20, Kjell R.)
      XI=DIMAG(X)
      XR = ABS(REAL(X)-Z*VD(J))
      ABZ = ABS(Z)
      IF(XI .GE. -ABZ) GOTO 3
C  Allow high frequency, cold plasma waves ...(ADDED 1993-04-02, KJELL R.)
      IF(XR .GT. 1.+5.*ABZ .AND. XR .GT. 1.+5.*P) GOTO 3
      RX = XR - INT(XR)
      IF(RX .GT. 0.5) RX = 1.-RX
      IF(XI .GE. -0.6*RX) GOTO 3
      IERR=1
      RETURN

    3 XSI(1,1)=A
      XSI(6,1)=A+2*A*VD(J)**2
      DO 4 L=1,8
      BL=B(L)
c  Modified to allow Z < 0. (1996-02-20, Kjell R.)
***** CL=C(L)
      if (z .lt. 0.) then
        CL=-C(L)
      else
        CL= C(L)
      endif
      DL=CL+VD(J)
      Y=X-DL*Z
      BLY=BL/Y
C        ****** EVALUATE THE R-FUNCTION ******
      CALL RYLA(Y,ALA,RC)
      XY=1.+A*Z*CL/Y
      PS=XY*RC(1,1)
      PSP=XY*RC(1,2)
C                     ****** FORM SUCEPTIBILITY TENSOR ******
      XSI(1,1)=XSI(1,1)+A*BL*Y*PS
      XSI(2,1)=XSI(2,1)+AI*BL*PSP
      XSI(3,1)=XSI(3,1)+A*P*BL*DL*PS
      XSI(4,1)=XSI(4,1)+BLY*PSP
CCCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
*****      XSI(5,1)=XSI(5,1)+AI*P*BLY*DL*PSP
      XSI(5,1)=XSI(5,1)-AI*P*BLY*DL*PSP
      XSI(6,1)=XSI(6,1)+2.*A*BLY*DL*DL*(X-VZ+AL*PS)
C
      IF(KOL.LE.1) GOTO 4
C                    ****** FORM X-DERIVATIVES OF XSI ******
      PSY=XY*RC(2,1)
      PPY=XY*RC(2,2)
      XY=X/Y
      DP =XY*(PSY-A*Z*CL/Y*RC(1,1))
      DPP=XY*(PPY-A*Z*CL/Y*RC(1,2))
      XSI(1,2)=XSI(1,2)+A*BL*(Y*DP+X*PS)
      XSI(2,2)=XSI(2,2)+AI*BL*DPP
      XSI(3,2)=XSI(3,2)+A*P*BL*DL*DP
      XSI(4,2)=XSI(4,2)+BLY*(DPP-XY*PSP)
CCCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
*****       XSI(5,2)=XSI(5,2)+AI*P*BLY*DL*(DPP-XY*PSP)
      XSI(5,2)=XSI(5,2)-AI*P*BLY*DL*(DPP-XY*PSP)
      XSI(6,2)=XSI(6,2)+2.*A*BLY*DL*DL*(X+AL*DP-XY*(X-VZ+AL*PS))
C
      IF(KOL.LE.2) GOTO 4
C                    ****** FORM Z-DERIVATIVES OF XSI ******
      ZY=Z/Y
      DP =ZY*(A*CL*XY*RC(1,1)-DL*PSY)
CCCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
*****       XSI(5,2)=XSI(5,2)+AI*P*BLY*DL*(DPP-XY*PSP)
      XSI(5,2)=XSI(5,2)-AI*P*BLY*DL*(DPP-XY*PSP)

      DPP=ZY*(A*CL*XY*RC(1,2)-DL*PPY)
      ZY=DL*ZY
      XSI(1,3)=XSI(1,3)+A*BL*Y*(DP-ZY*PS)
      XSI(2,3)=XSI(2,3)+AI*BL*DPP
      XSI(3,3)=XSI(3,3)+A*P*BL*DL*DP
      XSI(4,3)=XSI(4,3)+BLY*(DPP+ZY*PSP)
CCCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
*****       XSI(5,3)=XSI(5,3)+AI*P*BLY*DL*(DPP+ZY*PSP)
      XSI(5,3)=XSI(5,3)-AI*P*BLY*DL*(DPP+ZY*PSP)
      XSI(6,3)=XSI(6,3)+2.*A*BLY*DL*DL*(AL*DP-VZ+ZY*(X-VZ+AL*PS))
C
      IF(KOL.LE.3) GOTO 4
C                   ****** FORM P-DERIVATIVES OF XSI ******
      CALL RYLA(Y-1.,ALA,RC)
      DP=2.*(PSP-PS)
c      DPP=2.*AL*((Y/(Y-1.))**2*RC(1,2)-PSP)-Y*DP
CCCCCCCorrection, 1991 - 03 - 13, Kjell R 
c      DPP=2.*AL*((Y/(Y-1.))**2*(1.+A*Z*CL/Y)*RC(1,2)-PSP)-Y*DP
C     As pointed out by Scott Boardsen, NASA/MSFC, the AL above should
C     really be ALA.   Corrected  1991 - 12 - 18, Kjell R 
      DPP=2.*ALA*((Y/(Y-1.))**2*(1.+A*Z*CL/Y)*RC(1,2)-PSP)-Y*DP
      XSI(1,4)=XSI(1,4)+A*BL*Y*DP 
      XSI(2,4)=XSI(2,4)+AI*BL*DPP
      XSI(3,4)=XSI(3,4)+A*P*BL*DL*(DP+PS)
      XSI(4,4)=XSI(4,4)+BLY*(2.*PSP+DPP)
CCCCCCCCCCCCCCorrected sign error 1991 - 03 - 28, Kjell R.
*****       XSI(5,4)=XSI(5,4)+AI*P*BLY*DL*(PSP+DPP)
      XSI(5,4)=XSI(5,4)-AI*P*BLY*DL*(PSP+DPP)
      XSI(6,4)=XSI(6,4)+4.*A*BLY*DL*DL*AL*PSP
    4 CONTINUE
C      **** COMPLETE XSI(4, ) ****
      XSI(4,1)=XSI(1,1)-2.*A**2*AL*XSI(4,1)
      XSI(4,2)=XSI(1,2)-2.*A**2*AL*XSI(4,2)
      XSI(4,3)=XSI(1,3)-2.*A**2*AL*XSI(4,3)
      XSI(4,4)=XSI(1,4)-2.*A**2*AL*XSI(4,4)
      END
      