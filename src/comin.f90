module comin
   implicit none
   integer, parameter, private :: d2p = 8  !d2p=kind(1.0d0)

   real(kind=d2p) :: ASS(10)   ! particle mass in unit of mp, if 0 then use e- mass
   real(kind=d2p) :: AA(10, 2) ! alpha parameters in distribution function
   real(kind=d2p) :: DD(10)    ! Delta_j in distribution function
   real(kind=d2p) :: VD(10)    ! drift velocity normalized against parallel thermal velocity
   real(kind=d2p) :: DN(10)    ! number density of each plasma component cm^-3
   real(kind=d2p) :: TA(10)    ! temperature in keV
   real(kind=d2p) :: PM(3)     ! array (pmin,pmax,pstep)
   real(kind=d2p) :: ZM(3)     ! array (zmin,zmax,zstep)
   real(kind=d2p) :: XOI       ! initialization frequency (real part)
   real(kind=d2p) :: XC        ! XC  - gyrofrequency (kHZ)
   real(kind=d2p) :: PZL       ! option L, =0 linear input, =1 logarythmic input of P and Z
   real(kind=d2p) :: BETA      ! beta of the first species, ratio of thermal velocity to gyrofrequency
   real(kind=d2p) :: mi_o_me   ! mass ratio between ions and electrons typically set to 1836.1

   logical        :: printDebugInfo = .false.
   character(len=100) :: output_filename = 'whamp_output.txt'  ! Default output filename

   integer        :: cycleZFirst   ! which direction vary first (KFS in WHAMP) ?? how to handle

!  ------     below are variables calculated by WHAMP or WHAMP_ENGINE -----------
!
   integer           :: JMA    ! the number of last plasma species with non-zero density
   real(kind=d2p)    :: PP(10) ! ??
   real(kind=d2p)    :: ZZ(10) ! ??
   real(kind=d2p)    :: CV     ! speed of light / thermal velocity of 1st species
   complex(kind=d2p) :: XX(10)
   complex(kind=d2p) :: XP(10) ! wpj^2/w^2
   real(kind=d2p)    :: PX     ! plasma frequency
   real(kind=d2p)    :: DEN    ! total electron density

!  ------     root finding parameters
!
   integer(kind=4)      :: maxIterations = 50

end module comin

