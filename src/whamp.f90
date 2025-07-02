PROGRAM WHAMP
   use comin
   use comcout
   implicit none
   integer, parameter :: d2p = 8   !d2p=kind(1.0d0)

   integer           :: I, IERR, IRK, J, KFS
   logical           :: rootFindingConverged
   logical           :: solutionIsTooHeavilyDamped
   logical           :: isChangedPlasmaModel
   real(kind=d2p)    :: REN(10)    ! particle mass expressed in masses of first particles
   real(kind=d2p)    :: RN         ! mass of first particle in electron masses
   real(kind=d2p)    :: PXN        ! plasma frequency of species 1
   real(kind=d2p)    :: VTH        ! thermal velocity of the first species over speed of light
   real(kind=d2p)    :: REDN       ! density of first particle times m_e/m_s
   real(kind=d2p)    :: ADIR       ! abs(D)
   real(kind=d2p)    :: DEK, DET, DKP, DKZ, KV, PFQ, PLG, PLO, PO, PVO
   real(kind=d2p)    :: RED, ST, T, TR, XA, XI, ZLG, ZLO, ZO, ZVO
   COMPLEX(kind=d2p) :: XO, XVO, ddDX, OME, FPX, DOX, DOZ, DOP
   DIMENSION T(10), ST(10)
   integer           :: narg, iarg
   character(len=20) :: inputParameter, modelFilename

! Default plasma model
   DN = [1.0e6, 1.0e6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] !  density in m-3
   TA = [0.01, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ! temperature in eV
   DD = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ! loss cone parameter, default 1.0 (no loss cone)
   AA(:, 1) = [5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ! t_perp/t_par ratio, default 1.0
   AA(:, 2) = [0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ! default 0, i.e. no loss cone
   ASS = [16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] !/ 0-electrons, 1-protons, 16-oxygen
   VD = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ! v_drift/v_term
   XC = 2.79928 ! electron gyrofrequency in kHz 
   PZL = 0.0 ! 1 - log scale, 0 - linear scale
   mi_o_me = 1836.1_d2p ! default ion to electron mass ratio
   cycleZFirst = 1
   PM = [0.0, 0.0, 10.0]
   ZM = [0.0, 0.0, 10.0]
   XOI = .1

! Check command line input parameters
   narg = command_argument_count()
   iArg = 1
   do
      if (iArg > narg) exit
      call get_command_argument(iArg, inputParameter)
      if (printDebugInfo) write (*, *) "Input parameter:", inputParameter
      select case (adjustl(inputParameter))
      case ("-help", "-h", "--help")
         write (*, *) "usage: whamp [-help] [-debug] [-maxiterations <number>] [-file <modelFilename>] [-outfile <outputFilename>]"
         stop
      case ("-debug")
         printDebugInfo = .true.
         if (printDebugInfo) write (*, *) "Enable debugging"
      case ("-file")
         if (iArg == narg) then
            write (*, *) "ERROR: File name not given"
            stop
         end if
         iArg = iArg + 1
         call get_command_argument(iArg, modelFilename)
         if (printDebugInfo) write (*, *) "Reading file: ", modelFilename
         call read_input_file(modelFilename)
      case ("-outfile")
         if (iArg == narg) then
            write (*, *) "ERROR: Output file name not given"
            stop
         end if
         iArg = iArg + 1
         call get_command_argument(iArg, output_filename)
         if (printDebugInfo) write (*, *) "Output file: ", trim(output_filename)
      case ("-maxiterations")
         if (iArg == narg) then
            write (*, *) "ERROR: maxiterations is not specified"
            stop
         end if
         iArg = iArg + 1
         call get_command_argument(iArg, inputParameter)
         read (inputParameter, *) maxIterations
         if (printDebugInfo) write (*, *) "Max iterations: ", maxIterations
      case default
         if (printDebugInfo) write (*, *) "Option '", trim(inputParameter), "' is unknown"
      end select
      iArg = iArg + 1
   end do

   !
   IERR = 0

   loop_plasma_update: do
      isChangedPlasmaModel = .false. ! changed to .true. in code when new plasma parameters are entered
      DEN = 0.d+0 ! total electron density?
      RED = 0.d+0
      loop_species: DO J = 1, 10
         REN(J) = mi_o_me*ASS(J) ! mass of the particle normalized to m_e, TODO: check for hardcoded m_i/m_e
         IF (REN(J) .EQ. 0.) REN(J) = 1.
         T(J) = TA(J)/TA(1) ! temperature of species normalized by species 1 temperature
         IF (DN(J) .EQ. 0.) cycle loop_species ! ignore 0 density species
         JMA = J ! count how many species there are until we hit 0 in DN
         RED = RED + DN(J)/REN(J) ! m_e*( n_e/m_e + n_i/m_i + ... )
         IF (ASS(J) .EQ. 0.) DEN = DEN + DN(J)
      end do loop_species

      
      !
      RN = REN(1) ! mass of the first species particle normalized to m_e
      REDN = DN(1)/RN ! density of the first species times m_e/m_s
      !                  ****  NORMALIZED TEMPERATURES AND VELOCITIES.  ****
      DO J = 1, JMA
         REN(J) = REN(J)/RN ! normalize mass of the particle to the first species particle
         T(J) = T(J)*REN(J) ! Temperature times normalized mass (by species 1)
         ST(J) = SQRT(T(J)) ! sqrt of normalized temperature (by species 1) times charge
      end do
      ! We are computing frequncy in kHz which is why we need 10^6 multiplier and (2pi)^2
      DEK = 12405._d2p ! 10^6 * 4 pi^2 m_e epsilon_0 / e^2 = 12404.4
      PFQ = RED/DEK
      PX = SQRT(PFQ) ! plasma frequency measured in Hz = sum_s (n_s e^2)/(m_s epsilon_0) (see Fitzpatrick Introduction to Plasmas, equation 5.33)
      PXN = SQRT(REDN/DEK) !  plasma frequency of species 1
      XA = XC/RN ! gyrofrequency of the first species 1 to which dispersion relation will be normalized
      TR = TA(1)/RN ! temperature of the first species divided by its mass over m_e ratio
      CV = TR*(1022._d2p + TR)/(511._d2p + TR)**2
      CV = 1./SQRT(CV)
      DEK = DEK*RN
      DET = 255.499_d2p ! 0.001 c^2 m_e/ 2 e (becuse T is measured in keV)
      VTH = SQRT(TR/DET) ! thermal velocity over speed of light of the first species
      BETA = (VTH*PXN/XA)**2 ! beta of the first species, ratio of thermal velocity to gyrofrequency
      !
      call print_plasma_parameters()
      !                  ****  ASK FOR INPUT!  ****
      loop_typin: do
         !for new plasma skip calling typin until convergence checked
         if (.not. isChangedPlasmaModel) call TYPIN(isChangedPlasmaModel, KFS)
         if (isChangedPlasmaModel) cycle loop_plasma_update

         ! Add quit check here
         if (KFS == 999) then  ! Use 999 as quit signal from TYPIN
            write(*, *) 'Quitting program...'
            stop
         end if

         isChangedPlasmaModel = .false.
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
            if (.not. rootFindingConverged) then
               PRINT 125, P, Z, X, I, IRK
125            FORMAT(2X, 'NO CONVERGENCE!'/'  KP=', F6.3, '  KZ=',&
                      & F6.4, '  X=', E12.2, E12.2/'  I=', I3, '  IRK=', I3/)
               IF (KFS .EQ. 1) PLG = 1.D99 ! end cycling in P
               IF (KFS .EQ. 2) ZLG = 1.D99 ! end cycling in Z
            end if
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
               CALL OUTPT
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
            end if
            if (solutionIsTooHeavilyDamped) then
               PRINT *, ' TOO HEAVILY DAMPED!'
               PRINT *, '   '
               IERR = 0
               CALL OUTPT
               IF (KFS .EQ. 1) PLG = 1.D99
               IF (KFS .EQ. 2) ZLG = 1.D99
            end if

            if (KFS == 1) then        ! cycle first P
               PLG = PLG + PM(3)
               !                   **** UPDATE P AND Z.  ****
               if (PLG .GE. PM(1) .AND. PLG .LE. PM(2)) then
                  P = PLG + PZL*(10.**PLG - PLG)
               else
                  ZLG = ZLG + ZM(3)
                  print *
                  if (ZLG .LT. ZM(1) .OR. ZLG .GT. ZM(2)) then
                     cycle loop_typin
                  end if
                  KV = 1
                  PLG = PLO
                  P = PVO
                  Z = ZLG + PZL*(10.**ZLG - ZLG)
               end if
            else if (KFS == 2) then ! cycle first Z
               ZLG = ZLG + ZM(3)
               IF (ZLG .GE. ZM(1) .AND. ZLG .LE. ZM(2)) then
                  Z = ZLG + PZL*(10.**ZLG - ZLG)
               else
                  PLG = PLG + PM(3)
                  print *
                  IF (PLG .LT. PM(1) .OR. PLG .GT. PM(2)) cycle loop_typin
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
      end do loop_typin
   end do loop_plasma_update
contains
   subroutine print_plasma_parameters
      !                  ****  PRINT PLASMA PARAMETERS.  ****
      PRINT 101, PX, XA, PXN
101   FORMAT('# TOTAL PLASMA FREQ.:', F11.5,&
             & ' KHZ; SPEC 1 GYRO FREQ.:', F10.5, ' KHZ; ',&
             & 'SPEC 1 PLASMA FREQ.:', F11.5, ' KHZ; ')
      PRINT 102, VTH, BETA
102   FORMAT('# SPECIES 1 parallel V_TH/C:', F11.6, &
             & ';  SPECIES 1 parallel BETA:', F11.6)
      DO J = 1, JMA
103      FORMAT('# ', A3, '  DN=', 1PE12.5, '  T=', 0PF9.5, '  D=', F4.2,&
                &'  A=', F4.2, '  B=', F4.2, ' VD=', F5.2)
         PRINT 103, species_symbol(ASS(J)), DN(J), TA(J), DD(J), AA(J, 1), AA(J, 2), VD(J)
      end do
      !
   end subroutine
   subroutine root_finding
      ! Newton's iteration method with small adjustments:
      ! 1) convergence criteria on relative and absolute size of CX
      ! 2) convergence criteria on relative change in D

      complex(kind=d2p)         :: CX                ! correction in Newton's iteration method

      CALL DIFU(2, JMA, IERR)
      IF (IERR .NE. 0) solutionIsTooHeavilyDamped = .true.
      !                  ****  START OF ITERATION.  ****
      if (printDebugInfo) write (*, *) 'START:', '. X=', X, 'D=', D, 'DX=', DX ! DEBUG
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
            if (printDebugInfo) &
                & write (*, '(I2,A ,I2 ,A,2E16.8,A,2E16.8 ,A,2E16.8,A,2E16.8)')&
                         & I, '.', IRK, '. X=', X, ' CX=', CX, ' D=', D, ' DX=', DX
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
                     PRINT *, ' Local minima!'
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
   pure function species_symbol(mass) result(symbol)
      real(kind=d2p), intent(in) :: mass
      character(5)   :: symbol
      if (mass == 0) then
         symbol = 'e-'
      else if (mass == 1) then
         symbol = 'H+'
      else if (mass == 2) then
         symbol = 'He++'
      else if (mass == 4) then
         symbol = 'He+'
      else if (mass == 16) then
         symbol = 'O+'
      else
         write (symbol, '(a,I3)') 'm=', mass
      end if
   end function
end program WHAMP
