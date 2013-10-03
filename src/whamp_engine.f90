SUBROUTINE WHAMP_ENGINE
  use comin
  use comcout
  use comoutput
  implicit none
  !integer, parameter :: d2p=kind(1.0d0)
  integer, parameter :: d2p=8

  integer :: I,IERR,IRK,J,KFS
  logical :: rootFindingConverged 
  logical :: solutionIsTooHeavilyDamped
  real(kind=d2p) :: DEN        ! total electron density
  real(kind=d2p) :: REN(10)    ! particle mass expressed in masses of first particles
  real(kind=d2p) :: RN         ! mass of first particle in electron masses
  real(kind=d2p) :: ADIR       ! abs(D)
  real(kind=d2p) :: DEK,DKP,DKZ,KV,PFQ,PLG,PLO,PO,PVO,PX
  real(kind=d2p) :: RED, ST, T, TR, XA, XI, ZLG, ZLO, ZO,ZVO
  COMPLEX(kind=d2p) :: XO,XVO,ddDX,OME,FPX, DOX,DOZ,DOP
  DIMENSION T(10),ST(10)

  !CALL READ_INPUT_FILE(FILENAME)
  !
  IERR=0

  call allocate_output_matrices
  call plasma_setup
  call input_setup
     typin_loop: do 
        KV=1
        PLG=PM(1)
        IF(PM(3).LT.0.) PLG=PM(2)
        ZLG=ZM(1)
        IF(ZM(3).LT.0.) ZLG=ZM(2)
        IF(PZL==1.) then
           P=10.**PLG
           Z=10.**ZLG
        else
           P=PLG
           Z=ZLG
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
           solutionIsTooHeavilyDamped = .false. ! default not heaviily damped
           rootFindingConverged       = .false. ! default no convergence       
           call root_finding
           !
           if (.not. rootFindingConverged ) then
              PRINT 125,P,Z,X,I,IRK
125           FORMAT(2X,'NO CONVERGENCE!'/'  KP=',F6.3,'  KZ=',&
                   &  F6.4,'  X=',E12.2,E12.2/'  I=',I3,'  IRK=',I3/)
              IF(KFS .EQ. 1) PLG = 1.D99 ! end cycling in P
              IF(KFS .EQ. 2) ZLG = 1.D99 ! end cycling in Z
           end if
           if ((rootFindingConverged) .AND. (.not.solutionIsTooHeavilyDamped)) then
              !                  ****  CONVERGENCE!  ****
              CALL DIFU(4,JMA,IERR)
              !
              XI=DIMAG(X)
              VG(1)=-real(DP/DX)
              VG(2)=-real(DZ/DX)
              RI=SQRT(P**2+Z**2)*CV/real(X)  ! refractive index
              IF(VG(1).NE.0.) SG(1)=XI/VG(1)
              IF(VG(2).NE.0.) SG(2)=XI/VG(2)
              !          ****  PRINT THE RESULTS.  ****
              !              CALL OUTPT
              call save_output
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
           if (solutionIsTooHeavilyDamped) then
              PRINT*,' TOO HEAVILY DAMPED!'
              PRINT*,'   '
              IERR=0
              !CALL OUTPT
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
  contains
  subroutine print_plasma_parameters
 !                  ****  PRINT PLASMA PARAMETERS.  ****
     PRINT 101,PX,XC,DEN
101  FORMAT('# PLASMA FREQ.:',  F11.4,&
          &       'KHZ GYRO FREQ.:',  F10.4,  'KHZ   ',&
          &       'ELECTRON DENSITY:',1PE11.5,  'M-3'   )
     DO  J=1,JMA
102     FORMAT('# ',  A3,  '  DN=',1PE12.5,  '  T=',0PF9.5,  '  D=',  F4.2,&
             &'  A=',  F4.2,  '  B=',  F4.2,  ' VD=',  F5.2)
        PRINT 102,species_symbol(ASS(J)),DN(J),TA(J),DD(J),AA(J,1),AA(J,2),VD(J)
     end do
     !
  end subroutine
  subroutine plasma_setup
          DEN=0.d+0
          RED=0.d+0
          species_loop: DO J=1,10
              REN(J)=1836.1*ASS(J)
              IF(REN(J).EQ.0.) REN(J)=1.
              T(J)=TA(J)/TA(1)
              IF(DN(J).EQ.0.) cycle species_loop
              JMA=J
              RED=RED+DN(J)/REN(J)
              IF( ASS(J).EQ.0.) DEN=DEN+DN(J) 
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
          !
          call print_plasma_parameters()
          !
  end subroutine
  subroutine input_setup
          KFS = zfirst

  end subroutine
  subroutine root_finding
      ! Newton's iteration method with small adjustments:
      ! 1) First 2 steps are taken only along real direction to avoid
      ! convergence to some large imaginary frequency solutions
      ! 2) find roots of (w^2 D) instead of D, thus making slightly faster
      ! convergence
      ! 3) convergence criteria both on relative and absolute size of CX and on
      ! relative change in D
      integer(kind=4),parameter :: maxIterations=50
      complex(kind=d2p) :: CX ! correction in Newton's iteration method
      
      CALL DIFU(2,JMA,IERR)
      IF(IERR.NE.0) solutionIsTooHeavilyDamped = .true.
      !                  ****  START OF ITERATION.  ****
      if (printDebugInfo) write(*,*) 'START:','. X=',X,'D=',D,'DX=',DX ! DEBUG
      iteration_loop: DO I=1,maxIterations
          ADIR=ABS(D)
          IRK=1
          if (I < 2) then ! first step make only in real direction
              CX=REALPART(D)/REALPART(DX)
          else
              CX=D/DX
          end if
          CX=CX*X/(2*CX+X) ! finding zero of (w^2 D), faster convergence
          irk_loop: do
              X=X-CX
              OME=(X*XA)**2
              FPX=PFQ/OME
              DO  J=1,JMA
                  XP(J)=DN(J)/DEK/REN(J)/OME
                  XX(J)=X*REN(J)
              end do
              CALL DIFU(2,JMA,IERR)
              if (printDebugInfo) &
                  & write(*,'(I2,A ,I2 ,A,2E16.8,A,2E16.8 ,A,2E16.8,A,2E16.8)')&
                            & I,'.',IRK,'. X=',X,' CX=',CX,' D=',D ,' DX=',DX ! DEBUG
              IF(IERR.NE.0) then
                  solutionIsTooHeavilyDamped = .true.
                  exit iteration_loop
              end if
              IF(ABS(D).LT.ADIR) then
                  if( (ABS(CX).LE.1.E-6*ABS(X)) &    ! relative frequency precision
                      & .or. (ABS(CX) < 1e-6) ) then ! absolute precision
                      rootFindingConverged = .true.
                      if (I >= 2) then ! at least 2 steps have been made
                          exit iteration_loop
                      else
                          cycle iteration_loop
                      end if
                  else
                      cycle iteration_loop
                  end if 
              else
                  X=X+CX
                  CX=CX/2.
                  if (IRK > 3) then ! check for sitting at local minima
                      if ((ABS(D)-ADIR)/ADIR<1e-5) then
                          rootFindingConverged = .false.
                          exit iteration_loop
                      end if
                  end if
                  IRK=IRK+1
                  IF(IRK.GT.maxIterations) exit iteration_loop
              end if
          end do irk_loop
      end do iteration_loop
  end subroutine
  subroutine allocate_output_matrices
          ! estimate the size of matrices
          integer :: perpSize
          integer :: parSize
          integer :: i
          if (size(PM) == 1) then ! scalar
              perpSize = 1
              allocate (kperpOUT(1))
              kperpOUT = PM
          elseif (size(PM) == 3) then ! vector
              perpSize = 1 + floor((max(PM(1),PM(2))-min(PM(1),PM(2)))*sign(PM(3),1.0d0)/PM(3))
              allocate (kperpOUT(perpSize))
              do i=1,perpSize
                  kperpOUT(i) = PM(1)+(i-1.0)*PM(3)
              end do
          end if
          if (size(ZM) == 1) then ! scalar
              parSize = 1
              kparOUT = ZM
          elseif (size(ZM) == 3) then ! vector
              parSize = 1 + floor((max(ZM(1),ZM(2))-min(ZM(1),ZM(2)))*sign(ZM(3),1.0d0)/ZM(3))
              allocate (kparOUT(parSize))
              do i=1,parSize
                  kparOUT(i) = ZM(1)+(i-1.0)*ZM(3)
              end do
          end if
          allocate (fOUT(perpSize,parSize))
          allocate (ExOUT(perpSize,parSize))
          allocate (EyOUT(perpSize,parSize))
          allocate (EzOUT(perpSize,parSize))
          allocate (BxOUT(perpSize,parSize))
          allocate (ByOUT(perpSize,parSize))
          allocate (BzOUT(perpSize,parSize))
          allocate (SxOUT(perpSize,parSize))
          allocate (SyOUT(perpSize,parSize))
          allocate (SzOUT(perpSize,parSize))
          allocate (EBOUT(perpSize,parSize))
          allocate (VGPOUT(perpSize,parSize))
          allocate (VGZOUT(perpSize,parSize))
          allocate (SGPOUT(perpSize,parSize))
          allocate (SGZOUT(perpSize,parSize))
          allocate (uOUT(perpSize,parSize))
          allocate (flagSolutionFoundOUT(perpSize,parSize))
          allocate (flagTooHeavilyDampedOUT(perpSize,parSize))
          allocate (flagNoConvergenceOUT(perpSize,parSize))
          flagSolutionFoundOUT = 0
  end subroutine
  subroutine save_output
          integer :: indexKperp, indexKpar
          if (PLG == PM(1)) then
              indexKperp = 1
          else
              indexKperp = 1 + nint((PLG-PM(1))/PM(3))
          endif
          if (ZLG == ZM(1)) then
              indexKpar = 1
          else
              indexKpar = 1 + nint((ZLG-ZM(1))/ZM(3))
          endif

          if (rootFindingConverged) then 
              flagSolutionFoundOUT(indexKperp,indexKpar) = 1
              fOUT(indexKperp,indexKpar) = X
              ExOUT(indexKperp,indexKpar) = EFL(1) 
              EyOUT(indexKperp,indexKpar) = EFL(2) 
              EyOUT(indexKperp,indexKpar) = EFL(3) 
              BxOUT(indexKperp,indexKpar) = BFL(1) 
              ByOUT(indexKperp,indexKpar) = BFL(2) 
              ByOUT(indexKperp,indexKpar) = BFL(3) 
              VGPOUT(indexKperp,indexKpar) = VG(1) 
              VGZOUT(indexKperp,indexKpar) = VG(2) 
              SGPOUT(indexKperp,indexKpar) = SG(1) 
              SGZOUT(indexKperp,indexKpar) = SG(2) 
              uOUT(indexKperp,indexKpar) = ENE
          elseif (solutionIsTooHeavilyDamped) then
              flagTooHeavilyDampedOUT(indexKperp,indexKpar) = 1
          else
              flagNoConvergenceOUT(indexKperp,indexKpar) = 1
          endif
  end subroutine
  pure function species_symbol(mass) result(symbol)
          real(kind=d2p),intent(in) :: mass
          character(5)   :: symbol
          if (mass==0) then
              symbol = 'e-'
          else if (mass==1) then
              symbol = 'H+'
          else if (mass==2) then
              symbol = 'He++'
          else if (mass==4) then
              symbol = 'He+'
          else if (mass==16) then
              symbol = 'O+'
          else
             write(symbol,'(a,I3)') 'm=',mass 
          endif
  end function
end subroutine WHAMP_ENGINE
