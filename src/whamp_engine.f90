SUBROUTINE WHAMP_ENGINE
   use comin
   use comcout
   use comoutput
   implicit none
   !integer, parameter :: d2p=kind(1.0d0)
   integer, parameter :: d2p = 8

   integer :: I, IERR, IRK, J
   logical :: rootFindingConverged
   logical :: solutionIsTooHeavilyDamped
!  real(kind=d2p) :: DEN        ! total electron density
   real(kind=d2p) :: REN(10)    ! particle mass expressed in masses of first particles
   real(kind=d2p) :: RN         ! mass of first particle in electron masses
   real(kind=d2p) :: ADIR       ! abs(D)
   real(kind=d2p) :: DEK, DKP, DKZ, KV, PFQ, PLG, PLO, PO, PVO
   !real(kind=d2p) :: PX
   real(kind=d2p) :: RED, ST, T, TR, XA, XI, ZLG, ZLO, ZO, ZVO
   COMPLEX(kind=d2p) :: XO, XVO, ddDX, OME, FPX, DOX, DOZ, DOP
   DIMENSION T(10), ST(10)
   real(kind=d2p) :: coef_poynt ! coefficient used to simplify Poynting flux estimate
   real(kind=d2p), parameter :: PI = 3.14159265358979_d2p

   !CALL READ_INPUT_FILE(FILENAME)
   !
   IERR = 0
   call allocate_output_matrices
   call plasma_setup
   if ((PM(1) == 0.0) .and. (ZM(1) == 0.0)) return
   call input_setup

   KV = 1
   PLG = PM(1)
   IF (PM(3) .LT. 0.) PLG = PM(2)
   ZLG = ZM(1)
   IF (ZM(3) .LT. 0.) ZLG = ZM(2)
   IF (PZL == 1.) then
      P = 10.**PLG
      Z = 10.**ZLG
   else
      P = PLG
      Z = ZLG
   end if
   X = XOI
   loop_z_p: do
      OME = (X*XA)**2
      FPX = PFQ/OME
      DO J = 1, JMA
         XX(J) = X*REN(J)
         PP(J) = P*ST(J)
         ZZ(J) = Z*ST(J)
         XP(J) = DN(J)/DEK/REN(J)/OME
      end do
      !
      solutionIsTooHeavilyDamped = .false. ! default not heavily damped
      rootFindingConverged = .false. ! default no convergence
      call root_finding
      !
      if ((rootFindingConverged) .AND. (.not. solutionIsTooHeavilyDamped)) then
         !                  ****  CONVERGENCE!  ****
         CALL DIFU(4, JMA, IERR)
         !
         XI = DIMAG(X)
         VG(1) = -real(DP/DX)
         VG(2) = -real(DZ/DX)
         RI = SQRT(P**2 + Z**2)*CV/X  ! refractive index
         IF (VG(1) .NE. 0.) SG(1) = XI/VG(1)
         IF (VG(2) .NE. 0.) SG(2) = XI/VG(2)
         !          ****  PRINT THE RESULTS.  ****
         !if (printDebugInfo)              CALL OUTPT
         call save_output
         PO = P
         ZO = Z
         XO = X
         IF (KV /= 0) then
            XVO = X
            ZVO = Z
            ZLO = ZLG
            PVO = P
            PLO = PLG
            DOX = DX
            DOZ = DZ
            DOP = DP
            KV = 0
         end if
      else
         call save_output
      end if
      if (.not. rootFindingConverged) then
         !if (printDebugInfo) PRINT 125,P,Z,X,I,IRK
         !125           FORMAT(2X,'NO CONVERGENCE!'/'  KP=',F6.3,'  KZ=',&
         !                     & F6.4,'  X=',E12.2,E12.2/'  I=',I3,'  IRK=',I3/)
         IF (cycleZFirst == 1) PLG = 1.D99 ! end cycling in P
         IF (cycleZFirst == 2) ZLG = 1.D99 ! end cycling in Z
      end if
      if (solutionIsTooHeavilyDamped) then
         !if (printDebugInfo) PRINT*,' TOO HEAVILY DAMPED!'
         IERR = 0
         !if (printDebugInfo) CALL OUTPT
         IF (cycleZFirst .EQ. 1) PLG = 1.D99
         IF (cycleZFirst .EQ. 2) ZLG = 1.D99
      end if
      if (cycleZFirst == 0) then        ! cycle first P
         PLG = PLG + PM(3)
         !                   **** UPDATE P AND Z.  ****
         if (PLG .GE. PM(1) .AND. PLG .LE. PM(2)) then
            P = PLG + PZL*(10.**PLG - PLG)
         else
            ZLG = ZLG + ZM(3)
            !print*
            if (ZLG .LT. ZM(1) .OR. ZLG .GT. ZM(2)) then
               exit loop_z_p
            end if
            KV = 1
            PLG = PLO
            P = PVO
            Z = ZLG + PZL*(10.**ZLG - ZLG)
         end if
      else if (cycleZFirst == 1) then ! cycle first Z
         ZLG = ZLG + ZM(3)
         IF (ZLG .GE. ZM(1) .AND. ZLG .LE. ZM(2)) then
            Z = ZLG + PZL*(10.**ZLG - ZLG)
         else
            PLG = PLG + PM(3)
            !print*
            IF (PLG .LT. PM(1) .OR. PLG .GT. PM(2)) exit loop_z_p
            KV = 1
            ZLG = ZLO
            Z = ZVO
            P = PLG + PZL*(10.**PLG - PLG)
         end if
      end if
      !                    ****  NEW START FREQUENCY.  ****
      IF (KV /= 0) then
         DKP = P - PVO
         DKZ = Z - ZVO
         ddDX = (DKP*DOP + DKZ*DOZ)/DOX
         X = XVO - ddDX
      else
         DKP = P - PO
         DKZ = Z - ZO
         ddDX = (DKP*DP + DKZ*DZ)/DX
         X = XO - ddDX
      end if
   end do loop_z_p
   return
contains
   subroutine plasma_setup
      DEN = 0.d+0
      RED = 0.d+0
      loop_species: DO J = 1, 10
         REN(J) = mi_o_me*ASS(J)
         IF (REN(J) .EQ. 0.) REN(J) = 1.
         T(J) = TA(J)/TA(1)
         IF (DN(J) .EQ. 0.) cycle loop_species
         JMA = J
         RED = RED + DN(J)/REN(J)
         IF (ASS(J) .EQ. 0.) DEN = DEN + DN(J)
      end do loop_species
      !
      RN = REN(1)
      !                  ****  NORMALIZED TEMPERATURES AND VELOCITIES.  ****
      DO J = 1, JMA
         REN(J) = REN(J)/RN
         T(J) = T(J)*REN(J)
         ST(J) = SQRT(T(J))
      end do
      !
      DEK = 12405._d2p
      PFQ = RED/DEK
      PX = SQRT(PFQ)
      XA = XC/RN
      TR = TA(1)/RN
      CV = TR*(1022._d2p + TR)/(511._d2p + TR)**2
      CV = 1./SQRT(CV)
      DEK = DEK*RN
      !
   end subroutine
   subroutine input_setup

   end subroutine
   subroutine root_finding
      ! Newton's iteration method with small adjustments:
      ! 1) convergence criteria on relative and absolute size of CX
      ! 2) convergence criteria on relative change in D

      complex(kind=d2p)         :: CX                ! correction in Newton's iteration method

      CALL DIFU(2, JMA, IERR)
      IF (IERR .NE. 0) solutionIsTooHeavilyDamped = .true.
      !                  ****  START OF ITERATION.  ****
      !if (printDebugInfo) write(*,*) 'START:','. X=',X,'D=',D,'DX=',DX ! DEBUG
      loop_iteration: DO I = 1, maxIterations
         ADIR = ABS(D)
         IRK = 1
         CX = D/DX
         !CX=CX*X/(2*CX+X) ! finding zero of (w^2 D), faster convergence
         irk_loop: do
            X = X - CX
            OME = (X*XA)**2
            FPX = PFQ/OME
            DO J = 1, JMA
               XP(J) = DN(J)/DEK/REN(J)/OME
               XX(J) = X*REN(J)
            end do
            CALL DIFU(2, JMA, IERR)
            !if (printDebugInfo) &
            !    & write(*,'(I2,A ,I2 ,A,2E16.8,A,2E16.8 ,A,2E16.8,A,2E16.8)')&
            !            & I,'.',IRK,'. X=',X,' CX=',CX,' D=',D ,' DX=',DX
            IF (IERR .NE. 0) then
               solutionIsTooHeavilyDamped = .true.
               exit loop_iteration
            end if
            IF (ABS(D) .LT. ADIR) then
               if ((ABS(CX) .LE. 1.E-6*ABS(X)) & ! relative frequency precision
                   & .or. (ABS(CX) < 1e-6)) then ! absolute precision
                  rootFindingConverged = .true.
                  if (I >= 2) then ! at least 2 steps have been made
                     exit loop_iteration
                  else
                     cycle loop_iteration
                  end if
               else
                  cycle loop_iteration
               end if
            else
               X = X + CX
               CX = CX/2.
               if (IRK > 3) then ! check for sitting at local minima
                  if ((ABS(D) - ADIR)/ADIR < 1e-5) then
                     !PRINT*,' Local minima!'
                     rootFindingConverged = .false.
                     exit loop_iteration
                  end if
               end if
               IRK = IRK + 1
               IF (IRK .GT. maxIterations) exit loop_iteration
            end if
         end do irk_loop
      end do loop_iteration
   end subroutine
   subroutine allocate_output_matrices
      ! estimate the size of matrices
      integer :: i
      if (allocated(kperpOUT)) deallocate (kperpOUT)
      if (allocated(kparOUT)) deallocate (kparOUT)
      if (allocated(fOUT)) deallocate (fOUT)
      if (allocated(ExOUT)) deallocate (ExOUT)
      if (allocated(EyOUT)) deallocate (EyOUT)
      if (allocated(EzOUT)) deallocate (EzOUT)
      if (allocated(BxOUT)) deallocate (BxOUT)
      if (allocated(ByOUT)) deallocate (ByOUT)
      if (allocated(BzOUT)) deallocate (BzOUT)
      if (allocated(SxOUT)) deallocate (SxOUT)
      if (allocated(SyOUT)) deallocate (SyOUT)
      if (allocated(SzOUT)) deallocate (SzOUT)
      if (allocated(EBOUT)) deallocate (EBOUT)
      if (allocated(VGPOUT)) deallocate (VGPOUT)
      if (allocated(VGZOUT)) deallocate (VGZOUT)
      if (allocated(SGPOUT)) deallocate (SGPOUT)
      if (allocated(SGZOUT)) deallocate (SGZOUT)
      if (allocated(uOUT)) deallocate (uOUT)
      if (allocated(flagSolutionFoundOUT)) deallocate (flagSolutionFoundOUT)
      if (allocated(flagTooHeavilyDampedOUT)) deallocate (flagTooHeavilyDampedOUT)
      if (allocated(flagNoConvergenceOUT)) deallocate (flagNoConvergenceOUT)

      if (PM(1) == PM(2)) then ! kperp is one value
         kperpSize = 1
         allocate (kperpOUT(1))
         kperpOUT = PM(1)
      else                     ! kperp is vector
         kperpSize = 1 + floor((max(PM(1), PM(2)) - min(PM(1), PM(2)))*sign(1.0d0, PM(3))/PM(3))
         allocate (kperpOUT(kperpSize))
         do i = 1, kperpSize
            kperpOUT(i) = PM(1) + (i - 1.0)*PM(3)
         end do
      end if

      if (ZM(1) == ZM(2)) then ! scalar
         kparSize = 1
         allocate (kparOUT(1))
         kparOUT = ZM(1)
      else ! vector
         kparSize = 1 + floor((max(ZM(1), ZM(2)) - min(ZM(1), ZM(2)))*sign(1.0d0, ZM(3))/ZM(3))
         allocate (kparOUT(kparSize))
         do i = 1, kparSize
            kparOUT(i) = ZM(1) + (i - 1.0)*ZM(3)
         end do
      end if

      allocate (fOUT(kperpSize, kparSize))
      allocate (ExOUT(kperpSize, kparSize))
      allocate (EyOUT(kperpSize, kparSize))
      allocate (EzOUT(kperpSize, kparSize))
      allocate (BxOUT(kperpSize, kparSize))
      allocate (ByOUT(kperpSize, kparSize))
      allocate (BzOUT(kperpSize, kparSize))
      allocate (SxOUT(kperpSize, kparSize))
      allocate (SyOUT(kperpSize, kparSize))
      allocate (SzOUT(kperpSize, kparSize))
      allocate (EBOUT(kperpSize, kparSize))
      allocate (VGPOUT(kperpSize, kparSize))
      allocate (VGZOUT(kperpSize, kparSize))
      allocate (SGPOUT(kperpSize, kparSize))
      allocate (SGZOUT(kperpSize, kparSize))
      allocate (uOUT(kperpSize, kparSize))
      allocate (flagSolutionFoundOUT(kperpSize, kparSize))
      allocate (flagTooHeavilyDampedOUT(kperpSize, kparSize))
      allocate (flagNoConvergenceOUT(kperpSize, kparSize))
      flagSolutionFoundOUT = 0
      flagTooHeavilyDampedOUT = 0
      flagNoConvergenceOUT = 0
   end subroutine
   subroutine save_output
      integer :: indexKperp, indexKpar
      if (PLG == PM(1)) then
         indexKperp = 1
      else
         indexKperp = 1 + nint((PLG - PM(1))/PM(3))
      end if
      if (ZLG == ZM(1)) then
         indexKpar = 1
      else
         indexKpar = 1 + nint((ZLG - ZM(1))/ZM(3))
      end if

      if (rootFindingConverged) then
         flagSolutionFoundOUT(indexKperp, indexKpar) = 1
         fOUT(indexKperp, indexKpar) = X
         ExOUT(indexKperp, indexKpar) = EFL(1)
         EyOUT(indexKperp, indexKpar) = EFL(2)
         EzOUT(indexKperp, indexKpar) = EFL(3)
         BxOUT(indexKperp, indexKpar) = BFL(1)
         ByOUT(indexKperp, indexKpar) = BFL(2)
         BzOUT(indexKperp, indexKpar) = BFL(3)
         !     write Poynting vector uW/m^2
         coef_poynt = 10.0/4.0/PI/2.0
         SxOUT(indexKperp, indexKpar) = real(efl(2)*conjg(bfl(3)) - efl(3)*conjg(bfl(2)))*coef_poynt
         SyOUT(indexKperp, indexKpar) = real(efl(3)*conjg(bfl(1)) - efl(1)*conjg(bfl(3)))*coef_poynt
         SzOUT(indexKperp, indexKpar) = real(efl(1)*conjg(bfl(2)) - efl(2)*conjg(bfl(1)))*coef_poynt
         EBOUT(indexKperp, indexKpar) = sqrt( &
             & (real(EFL(1)*CONJG(EFL(1)) + EFL(2)*CONJG(EFL(2)) + EFL(3)*CONJG(EFL(3)))) &
             & /(real(BFL(1)*CONJG(BFL(1)) + BFL(2)*CONJG(BFL(2)) + BFL(3)*CONJG(BFL(3)))) &
             & )
         VGPOUT(indexKperp, indexKpar) = VG(1)
         VGZOUT(indexKperp, indexKpar) = VG(2)
         SGPOUT(indexKperp, indexKpar) = SG(1)
         SGZOUT(indexKperp, indexKpar) = SG(2)
         uOUT(indexKperp, indexKpar) = ENE
      elseif (solutionIsTooHeavilyDamped) then
         flagTooHeavilyDampedOUT(indexKperp, indexKpar) = 1
      else
         flagNoConvergenceOUT(indexKperp, indexKpar) = 1
      end if
   end subroutine
end subroutine WHAMP_ENGINE
