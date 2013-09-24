!**************************************************************
!
! File whamp.f 
!       SUBROUTINE WHAMP(IREAD_FILE, FILENAME)
PROGRAM WHAMP 
  use comin
  use comcout
  implicit none
  integer, parameter :: d2p=kind(1.0d0)

  integer :: I,IERR,IRK,ISP,J,KFS
  integer :: flag_convergence,flag_too_heavily_damped,flag_new_plasma
  real(kind=d2p) :: ADIR,DEK,DEN,DKP,DKZ,KV,PFQ,PLG,PLO,PO,PVO,PX
  real(kind=d2p) :: RED,REN, RN, ST, T, TR, XA, XI, ZLG, ZLO, ZO,ZVO
  CHARACTER FILENAME*(80)
  COMPLEX(kind=d2p) :: XO,XVO,ddDX,OME,FPX, DOX,DOZ,DOP,CX
  DIMENSION REN(10),T(10),ST(10),ISP(10)
  CHARACTER SPE(5)*3

  DATA SPE/'E- ','H+ ','HE+','O+ ','   '/
  !     REAL CHARGE(6)
  !      DATA CHARGE/1.,1.,1.,1.,1.,1./
  ! 
  !      IF(IREAD_FILE.EQ.1)CALL READ_INPUT_FILE(FILENAME)   
  CALL READ_INPUT_FILE(FILENAME)
  !
  IERR=0

  plasma_update_loop: do
     flag_new_plasma=0 ! changed to 1 in code when new plasma parameters are entered
     DEN=0.d+0
     RED=0.d+0
     species_loop: DO J=1,10
        REN(J)=1836.1*ASS(J)
        IF(REN(J).EQ.0.) REN(J)=1.
        T(J)=TA(J)/TA(1)
        ISP(J)=SQRT(ASS(J))
        IF(ISP(J).LT.4) ISP(J)=ISP(J)+1
        IF(DN(J).EQ.0.) cycle species_loop
        JMA=J
        RED=RED+DN(J)/REN(J)
        IF( ISP(J).EQ.1.) DEN=DEN+DN(J)
     end do species_loop
     !
     RN=REN(1)
     !                  ****  NORMALIZED TEMPERATURES AND VELOCITIES.  ****
     DO J=1,JMA
        REN(J)=REN(J)/RN
        T(J)=T(J)*REN(J)
        ST(J)=SQRT(T(J))
     end do
     !
     DEK=12405._d2p
     PFQ=RED/DEK
     PX=SQRT(PFQ)
     XA=XC/RN
     TR=TA(1)/RN
     CV=TR*(1022._d2p+TR)/(511._d2p+TR)**2
     CV=1./SQRT(CV)
     DEK=DEK*RN
     !                  ****  PRINT PLASMA PARAMETERS.  ****
     PRINT 101,PX,XC,DEN
101  FORMAT('# PLASMA FREQ.:',  F11.4,&
          &       'KHZ GYRO FREQ.:',  F10.4,  'KHZ   ',&
          &       'ELECTRON DENSITY:',1PE11.5,  'M-3'   )
     DO  J=1,JMA
102     FORMAT('# ',  A3,  '  DN=',1PE12.5,  '  T=',0PF9.5,  '  D=',  F4.2,&
             &'  A=',  F4.2,  '  B=',  F4.2,  ' VD=',  F5.2)
        PRINT 102,SPE(ISP(J)),DN(J),TA(J),DD(J),AA(J,1),AA(J,2),VD(J)
     end do
     !
     !                  ****  ASK FOR INPUT!  **** 
     typin_loop: do 
        !for new plasma skip calling typin until convergence checked
        if(flag_new_plasma==0) CALL TYPIN(flag_new_plasma,KFS) 
        if(flag_new_plasma==1) cycle plasma_update_loop
        flag_new_plasma=0
        KV=1
        PLG=PM(1)
        IF(PM(3).LT.0.) PLG=PM(2)
        P=PLG
        ZLG=ZM(1)
        IF(ZM(3).LT.0.) ZLG=ZM(2)
        Z=ZLG
        IF(PZL==1.) then
           P=10.**PLG
           Z=10.**ZLG
        end if
        X=XOI
        z_p_loop: do
           OME=(X*XA)**2
           FPX=PFQ/OME
           DO  J=1,JMA
              XX(J)=X*REN(J)
              PP(J)=P*ST(J)
              ZZ(J)=Z*ST(J)
              XP(J)=DN(J)/DEK/REN(J)/OME
           end do
           !
           flag_too_heavily_damped=0;
           CALL DIFU(2,JMA,IERR)
           IF(IERR.NE.0) flag_too_heavily_damped=1
           !                  ****  START OF ITERATION.  ****
           flag_convergence=0 ! default assume no convergence       
           iteration_loop: DO I=1,20
              ADIR=ABS(D)
              IRK=0
              CX=D/DX
              irk_loop: do
                 X=X-CX
                 OME=(X*XA)**2
                 FPX=PFQ/OME
                 DO  J=1,JMA
                    XP(J)=DN(J)/DEK/REN(J)/OME
                    XX(J)=X*REN(J)
                 end do
                 IF(ABS(CX).LE.1.E-6*ABS(X)) then
                    flag_convergence=1
                    exit iteration_loop
                 else
                    CALL DIFU(2,JMA,IERR)
                    IF(IERR.NE.0) then
                       flag_too_heavily_damped=1
                       exit iteration_loop
                    end if
                    IF(ABS(D).LT.ADIR) cycle iteration_loop
                    X=X+CX
                    CX=CX/2.
                    IRK=IRK+1
                    IF(IRK.GT.20) exit iteration_loop
                 end if
              end do irk_loop
           end do iteration_loop
           !
           if (flag_convergence==0) then
              PRINT 125,P,Z,X,I,IRK
125           FORMAT(2X,'NO CONVERGENCE!'/'  KP=',F6.3,'  KZ=',&
                   &  F6.4,'  X=',E12.2,E12.2/'  I=',I3,'  IRK=',I3/)
              IF(KFS .EQ. 1) PLG = 1.D99 ! end cycling in P
              IF(KFS .EQ. 2) ZLG = 1.D99 ! end cycling in Z
           end if
           if ((flag_convergence==1) .AND. (flag_too_heavily_damped==0)) then
              !                  ****  CONVERGENCE!  ****
              CALL DIFU(4,JMA,IERR)
           end if
           if  ((flag_convergence==1) .AND. (flag_too_heavily_damped==0)) then
              X=X-D/DX
              !
              XI=DIMAG(X)
              VG(1)=-DP/DX
              VG(2)=-DZ/DX
              RI=SQRT(P**2+Z**2)*CV/X
              IF(VG(1).NE.0.) SG(1)=XI/VG(1)
              IF(VG(2).NE.0.) SG(2)=XI/VG(2)
              !          ****  PRINT THE RESULTS.  ****
              CALL OUTPT
              PO=P
              ZO=Z
              XO=X
              IF(KV /= 0) then
                 XVO=X
                 ZVO=Z
                 ZLO=ZLG
                 PVO=P
                 PLO=PLG
                 DOX=DX
                 DOZ=DZ
                 DOP=DP
                 KV=0
              end if
           end if
           if (flag_too_heavily_damped==1) then
              PRINT*,' TOO HEAVILY DAMPED!'
              PRINT*,'   '
              IERR=0
              CALL OUTPT
              IF(KFS .EQ. 1) PLG = 1.D99
              IF(KFS .EQ. 2) ZLG = 1.D99
           end if
           
           if (KFS == 1) then        ! cycle first P
              PLG=PLG+PM(3)
              !                   **** UPDATE P AND Z.  ****
              if(PLG.GE.PM(1).AND.PLG.LE.PM(2)) then
                 P=PLG+PZL*(10.**PLG-PLG)
              else 
                 ZLG=ZLG+ZM(3)
                 print*
                 if(ZLG.LT.ZM(1).OR.ZLG.GT.ZM(2)) then 
                    cycle typin_loop 
                 end if
                 KV=1
                 PLG=PLO
                 P=PVO
                 Z=ZLG+PZL*(10.**ZLG-ZLG)
              end if
           else if (KFS == 2) then ! cycle first Z
              ZLG=ZLG+ZM(3)
              IF(ZLG.GE.ZM(1).AND.ZLG.LE.ZM(2)) then
                 Z=ZLG+PZL*(10.**ZLG-ZLG)
              else 
                 PLG=PLG+PM(3)
                 print*
                 IF(PLG.LT.PM(1).OR.PLG.GT.PM(2)) cycle typin_loop
                 KV=1
                 ZLG=ZLO
                 Z=ZVO
                 P=PLG+PZL*(10.**PLG-PLG)
              end if
           end if
           !                    ****  NEW START FREQUENCY.  ****
           IF(KV /= 0) then 
              DKP=P-PVO
              DKZ=Z-ZVO
              ddDX=(DKP*DOP+DKZ*DOZ)/DOX
              X=XVO-ddDX
           else
              DKP=P-PO
              DKZ=Z-ZO
              ddDX=(DKP*DP+DKZ*DZ)/DX
              X=XO-ddDX
           end if
        end do z_p_loop
     end do typin_loop
  end do plasma_update_loop
end program WHAMP

SUBROUTINE EPSGRAD(XSI, DF, Q, J, IB) 
  use comin 

  !*****  This dummy routine, which is called from DIFU, replaces
  !*****  a subroutine needed in the ray-tracing version of the code.
  !*****  
  ! 
  RETURN
END SUBROUTINE EPSGRAD




