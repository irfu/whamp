***************************************************************
*
* File difu.f 

       SUBROUTINE DIFU(KOL,JMAX,IERR)
C        ARGUMENTS: KOL     =2 D AND ITS X-DERIVATIVE COMPUTED.
C                           =4 ALL DERIVATIVES AND WAVE FIELDS.
C                   JMAX    =1 - 6, NUMBER OF COMPONENTS.
C                   IERR    ERROR FLAG
      include 'comin.h'
C
      include 'comcout.h'
      COMPLEX*16 XSI(6,4), DF, U1, U2, U3, U12, U13, U32,
     #   A, B, C, DA, DB, DC, EE

C              *********** FORM DIELECTRIC TENSOR ************
      DO 1 K=1,4
      DO 1 I=1,6
    1 E(I,K)=(0.,0.)
      E(1,1)=1.
      E(4,1)=1.
      E(6,1)=1.
C
      JMA = JMAX 
      KOLLA =MIN(KOL, 4)      
      DO 5 J=1,JMA
      E(1,1)=E(1,1)-XP(J)
      E(4,1)=E(4,1)-XP(J)
      E(6,1)=E(6,1)-XP(J)
      IF(AA(J,1).NE.AA(J,2)) GOTO 2
      AA(J,2)=0.
      DD(J)=1.
    2 IB=1
      DF=XP(J)/(AA(J,1)*(AA(J,1)-AA(J,2)))
      Q=AA(J,1)-DD(J)*AA(J,2)
    3 CALL CHI(XSI,J,IB,KOL,IERR)
      IF(IERR.NE.0) RETURN 
      if (kol .eq. 5) call epsgrad(xsi, df, q, j, ib)
      DO 4 K=1,KOLLA
      DO 4 I=1,6
      E(I,K)=E(I,K)+DF*Q*XSI(I,K)
    4 CONTINUE
C
      IF(IB.EQ.2) GOTO 5
      Q=(DD(J)-1.)*AA(J,1)
      IF(Q.EQ.0.) GOTO 5
      IB=2
      GOTO 3
    5 CONTINUE
C                       *** DIELECTRIC TENSOR COMPUTED ***
C
C       ******* FORM REFRACTIVE INDEX, CV=SPEED OF LIGHT/THERM. SPEED. *
      U1=PP(1)*CV/XX(1)
      U3=ZZ(1)*CV/XX(1)
      U12=U1*U1
      U32=U3*U3
      U2=U12+U32
      U13=2.*U1*U3
C
C                      ******** FORM DISPERSION FUNCTION ********
      A=U12*E(1,1)+U13*E(3,1)+U32*E(6,1)
**************  Sign error corrected in Feb. 1989.   Kjell R.  ******
****  B=U2*(E(1,1)*E(6,1)-E(3,1)**2)+(U3*E(5,1)+U1*E(2,1))**2
      B=U2*(E(1,1)*E(6,1)-E(3,1)**2)+(U3*E(5,1)-U1*E(2,1))**2
      C=(E(1,1)*E(6,1)-E(3,1)**2)*E(4,1)+E(6,1)*E(2,1)**2
      C=C+(E(1,1)*E(5,1)+2.*E(2,1)*E(3,1))*E(5,1)
C
      D=(U2-E(4,1))*A-B+C
      IF(KOL.LE.1) RETURN
C        ****** COMPLETE X-DERIVATIVE OF DIELECTRIC TENSOR ******
      E(1,2)=E(1,2)-2.*(E(1,1)-1.)
      E(2,2)=E(2,2)-2.* E(2,1)
      E(3,2)=E(3,2)-2.* E(3,1)
      E(4,2)=E(4,2)-2.*(E(4,1)-1.)
      E(5,2)=E(5,2)-2.* E(5,1)
      E(6,2)=E(6,2)-2.*(E(6,1)-1.)
C          ****** X-DERIVATIVE OF DISPERSION FUNCTION *******
      DA=(E(1,2)-2.*E(1,1))*U12+U13*(E(3,2)-2.*E(3,1))+
     +   (E(6,2)-2.*E(6,1))*U32
      DB=2.*(U3*E(5,1)-U1*E(2,1))*(U3*(E(5,2)-E(5,1))-U1*(E(2,2)-E(2,1)
     +))+U2*((E(1,2)-E(1,1))*E(6,1)-2.*(E(3,2)-E(3,1))*E(3,1)+
     +(E(6,2)-E(6,1))*E(1,1))
      DC=(E(1,2)*E(6,1)+E(1,1)*E(6,2)-2.*E(3,1)*E(3,2))*E(4,1)
      DC=DC+(E(1,1)*E(6,1)-E(3,1)**2)*E(4,2)
      DC=DC+E(5,1)*(E(1,2)*E(5,1)+2.*E(1,1)*E(5,2)+
     +E(2,1)*E(3,2)+2.*E(2,2)*E(3,1))
      DC=DC+E(2,1)*(E(6,2)*E(2,1)+2.*E(6,1)*E(2,2)+
     +E(5,1)*E(3,2)+2.*E(5,2)*E(3,1))
C
      DX=((U2-E(4,1))*DA-(2.*U2+E(4,2))*A-DB+DC)/XX(1)
      IF(KOL.LE.2) RETURN
      DZ=(0.,0.)
      IF(ZZ(1).EQ.0.) GOTO 6
C         ****** Z-DERIVATIVE OF DISPERSION FUNCTION ******
      DA=U12*E(1,3)+U13*(E(3,3)+E(3,1))+U32*(E(6,3)+2.*E(6,1))
C      DB=2.*(U3*E(5,1)-U1*E(3,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,1))+
CCCCCCCorrection 1991- 03 - 11, Kjell R.
c      DB=2.*(U3*E(5,1)-U1*E(2,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,1))+ 
CCCCCCCorrection of error found by Scott Boardsen, NASA/MSFC. 1992-03-12, Kjell R.
      DB=2.*(U3*E(5,1)-U1*E(2,1))*(U3*(E(5,3)+E(5,1))-U1*E(2,3))+ 
     +2.*U32*(E(1,1)*E(6,1)-E(3,1)**2)+
     +U2*(E(1,3)*E(6,1)+E(1,1)*E(6,3)-2.*E(3,1)*E(3,3))
      DC=(E(1,3)*E(6,1)+E(1,1)*E(6,3)-2.*E(3,1)*E(3,3))*E(4,1)
      DC=DC+(E(1,1)*E(6,1)-E(3,1)**2)*E(4,3)
      DC=DC+E(5,1)*(E(1,3)*E(5,1)+2.*E(1,1)*E(5,3)+
     +E(2,1)*E(3,3)+2.*E(2,3)*E(3,1))
      DC=DC+E(2,1)*(E(6,3)*E(2,1)+2.*E(6,1)*E(2,3)+
     +E(5,1)*E(3,3)+2.*E(5,3)*E(3,1))
C
      DZ=((U2-E(4,1))*DA+(2.*U32-E(4,3))*A-DB+DC)/ZZ(1)
    6 IF(KOL.LE.3) RETURN
C        ****** P-DERIVATIVE OF DISPERSION FUNCTION ******
      DP=(0.,0.)
      IF(PP(1).EQ.0.) GOTO 7
C
      DA=U12*(E(1,4)+2.*E(1,1))+U13*(E(3,4)+E(3,1))+U32*E(6,4)
      DB=2.*(U3*E(5,1)-U1*E(2,1))*(U3*E(5,4)-U1*(E(2,4)+E(2,1)))+
     +2.*U12*(E(1,1)*E(6,1)-E(3,1)**2)+
     +U2*(E(1,4)*E(6,1)+E(1,1)*E(6,4)-2.*E(3,1)*E(3,4))
      DC=(E(1,4)*E(6,1)+E(1,1)*E(6,4)-2.*E(3,1)*E(3,4))*E(4,1)
      DC=DC+(E(1,1)*E(6,1)-E(3,1)**2)*E(4,4)
      DC=DC+E(5,1)*(E(1,4)*E(5,1)+2.*E(1,1)*E(5,4)+
     +E(2,1)*E(3,4)+2.*E(2,4)*E(3,1))
      DC=DC+E(2,1)*(E(6,4)*E(2,1)+2.*E(6,1)*E(2,4)+
     +E(5,1)*E(3,4)+2.*E(5,4)*E(3,1))
C
      DP=((U2-E(4,1))*DA+(2.*U12-E(4,4))*A-DB+DC)/PP(1)
C                     ******** COMPUTE ELECTRIC FIELD ********
    7 CALL ENERGY(U1,U3,U2,U12,U32)
C       THE ELECTRIC FIELD IS 1 MV/M.
C       THE MAGNETIC FIELD WILL BE IN GAMMA
      END












