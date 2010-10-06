***************************************************************
*
* File whamp.f 
c       SUBROUTINE WHAMP(IREAD_FILE, FILENAME)
      PROGRAM WHAMP 
      IMPLICIT REAL*8 ( A - H, O - Z )
      CHARACTER FILENAME*(80)
      COMPLEX*16 X,XO,XVO,XX(10),XP(10),DX,OME,FPX,DIR,DIX,DIZ,DIP,
     * EPS(6,4),DOX,DOZ,DOP,CX,EFL(3),BFL(3),RI
      DIMENSION REN(10),T(10),ST(10),ISP(10),ITID(7)
      CHARACTER SPE(5)*3
      COMMON /XPZ/ XX,PP(10),ZZ(10),A(10),B(10),D(10),ASS(10),VD(10),
     C   DN(10),TA(10),XP,CV,PM(3),ZM(3),XOI,XC,PZL
      COMMON /COUT/ X,P,Z,EFL,BFL,DIR,DIX,DIZ,DIP,EPS,VG(2),SG(2),RI,ENE
      DATA SPE/'E- ','H+ ','HE+','O+ ','   '/
c     REAL CHARGE(6)
c      DATA CHARGE/1.,1.,1.,1.,1.,1./
C 
c      IF(IREAD_FILE.EQ.1)CALL READ_INPUT_FILE(FILENAME)   
      CALL READ_INPUT_FILE(FILENAME)
C
      NPL=0 

      NPL=0
    1 DEN=0.
      RED=0.
       DO 2 J=1,10
      REN(J)=1836.1*ASS(J)
      IF(REN(J).EQ.0.) REN(J)=1.
      T(J)=TA(J)/TA(1)
      ISP(J)=SQRT(ASS(J))
      IF(ISP(J).LT.4) ISP(J)=ISP(J)+1
      IF(DN(J).EQ.0.) GOTO 2
      JMA=J
      RED=RED+DN(J)/REN(J)
      IF( ISP(J).EQ.1.) DEN=DEN+DN(J)
    2 CONTINUE
C
      RN=REN(1)
C                  ****  NORMALIZED TEMPERATURES AND VELOCITIES.  ****
      DO 3 J=1,JMA
      REN(J)=REN(J)/RN
      T(J)=T(J)*REN(J)
    3 ST(J)=SQRT(T(J))
C
      DEK=12405.
      PFQ=RED/DEK
      PX=SQRT(PFQ)
      XA=XC/RN
      TR=TA(1)/RN
      CV=TR*(1022.+TR)/(511.+TR)**2
      CV=1./SQRT(CV)
      DEK=DEK*RN
C                  ****  PRINT PLASMA PARAMETERS.  ****
C      CALL CLOCK(ITID)
C      PRINT 100,(ITID(I), I=7,2,-1)
  100 FORMAT( 1X,  'DATE ',  I4,  2 ( '-',  I2.2 ),
     A        2X,  'TIME:',  I2.2,  2 ( '.',  I2.2 )  /  )
      PRINT 101,PX,XC,DEN
  101 FORMAT('# PLASMA FREQ.:',  F11.4,
     #       'KHZ GYRO FREQ.:',  F10.4,  'KHZ   ',
     #       'ELECTRON DENSITY:',1PE11.5,  'M-3'   )
      DO 4 J=1,JMA
  102 FORMAT('# ',  A3,  '  DN=',1PE12.5,  '  T=',0PF9.5,  '  D=',  F4.2,
     #'  A=',  F4.2,  '  B=',  F4.2,  ' VD=',  F5.2)
    4 PRINT 102,SPE(ISP(J)),DN(J),TA(J),D(J),A(J),B(J),VD(J)
C
      IF(NPL.EQ.1) GOTO 6
C                  ****  ASK FOR INPUT!  ****
    5 CALL TYPIN(NPL,KFS)
      IF(NPL.EQ.1) GOTO 1
    6 NPL=0
      KV=1
      PLG=PM(1)
      IF(PM(3).LT.0.) PLG=PM(2)
      P=PLG
      ZLG=ZM(1)
      IF(ZM(3).LT.0.) ZLG=ZM(2)
      Z=ZLG
      IF(PZL.NE.1.) GOTO 7
      P=10.**PLG
      Z=10.**ZLG
    7 X=XOI
   10 OME=(X*XA)**2
      FPX=PFQ/OME
      DO 11 J=1,JMA
      XX(J)=X*REN(J)
      PP(J)=P*ST(J)
      ZZ(J)=Z*ST(J)
   11 XP(J)=DN(J)/DEK/REN(J)/OME
C
      CALL DIFU(2,JMA,IERR)
      IF(IERR.NE.0) GOTO 50
C                  ****  START OF ITERATION.  ****
      DO 20 I=1,20
      ADIR=ABS(DIR)
      IRK=0
      CX=DIR/DIX
   15 X=X-CX
      OME=(X*XA)**2
      FPX=PFQ/OME
      DO 16 J=1,JMA
      XP(J)=DN(J)/DEK/REN(J)/OME
   16 XX(J)=X*REN(J)
      IF(ABS(CX).LE.1.E-6*ABS(X)) GOTO 30
      CALL DIFU(2,JMA,IERR)
      IF(IERR.NE.0) GOTO 50
      IF(ABS(DIR).LT.ADIR) GOTO 20
      X=X+CX
      CX=CX/2.
      IRK=IRK+1
      IF(IRK.GT.20) GOTO 25
      GOTO 15
   20 CONTINUE
C
   25 PRINT 125,P,Z,X,I,IRK
  125 FORMAT(2X,'NO CONVERGENCE!'/'  KP=',F6.3,'  KZ=',
     +  F6.4,'  X=',E12.2,E12.2/'  I=',I3,'  IRK=',I3/)
      GOTO 55 
C                  ****  CONVERGENCE!  ****
   30 CALL DIFU(4,JMA,IERR)
      IF(IERR.NE.0) GOTO 50
      X=X-DIR/DIX
C
      XI=DIMAG(X)
      VG(1)=-DIP/DIX
      VG(2)=-DIZ/DIX
      RI=SQRT(P**2+Z**2)*CV/X
      IF(VG(1).NE.0.) SG(1)=XI/VG(1)
      IF(VG(2).NE.0.) SG(2)=XI/VG(2)
C          ****  PRINT THE RESULTS.  ****
   34 CALL OUTPT
      PO=P
      ZO=Z
      XO=X
      IF(KV.EQ.0) GOTO 35
      XVO=X
      ZVO=Z
      ZLO=ZLG
      PVO=P
      PLO=PLG
      DOX=DIX
      DOZ=DIZ
      DOP=DIP
      KV=0
   35 GOTO(36,38) KFS
   36 PLG=PLG+PM(3)
C                   **** UPDATE P AND Z.  ****
      IF(PLG.GE.PM(1).AND.PLG.LE.PM(2)) GOTO 39
      ZLG=ZLG+ZM(3)
      print*
      IF(ZLG.LT.ZM(1).OR.ZLG.GT.ZM(2)) GOTO 5
      KV=1
      PLG=PLO
      P=PVO
   37 Z=ZLG+PZL*(10.**ZLG-ZLG)
      GOTO 40
C
   38 ZLG=ZLG+ZM(3)
      IF(ZLG.GE.ZM(1).AND.ZLG.LE.ZM(2)) GOTO 37
      PLG=PLG+PM(3)
      print*
      IF(PLG.LT.PM(1).OR.PLG.GT.PM(2)) GOTO 5
      KV=1
      ZLG=ZLO
      Z=ZVO
   39 P=PLG+PZL*(10.**PLG-PLG)
C                    ****  NEW START FREQUENCY.  ****
   40 IF(KV.EQ.0) GOTO 41
      DKP=P-PVO
      DKZ=Z-ZVO
      DX=(DKP*DOP+DKZ*DOZ)/DOX
      X=XVO-DX
      GOTO 10
   41 DKP=P-PO
      DKZ=Z-ZO
      DX=(DKP*DIP+DKZ*DIZ)/DIX
      X=XO-DX
      GOTO 10
   50 PRINT*,' TOO HEAVILY DAMPED!'
      PRINT*,'   '
      IERR=0
      CALL OUTPT
   55 IF(KFS .EQ. 1) PLG = 1.D99
      IF(KFS .EQ. 2) ZLG = 1.D99
      GOTO 35
   60 CONTINUE
      END 



      SUBROUTINE EPSGRAD(XSI, DF, Q, J, IB) 
      include 'comin.h' 

******  This dummy routine, which is called from DIFU, replaces
******  a subroutine needed in the ray-tracing version of the code.
******  
C 
       RETURN
       END




