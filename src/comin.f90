module comin
implicit none
!integer, parameter,private :: d2p=kind(1.0d0)
integer, parameter,private :: d2p=8

integer :: JMA             ! the number of last plasma species with non-zero density

real(kind=d2p) :: PP(10)   ! ??
real(kind=d2p) :: ZZ(10)   ! ??
real(kind=d2p) :: AA(10,2) ! alpha parameters in distribution function
real(kind=d2p) :: DD(10)   !
real(kind=d2p) :: ASS(10)  ! mass of components
real(kind=d2p) :: VD(10)   ! drift velocity 
real(kind=d2p) :: DN(10)   ! density of each plasma component cm^-3
real(kind=d2p) :: TA(10)   ! temperature in keV
real(kind=d2p) :: CV       ! speed of light / thermal velocity of 1st species
real(kind=d2p) :: PM(3)    ! array (pmin,pmax,pstep)
real(kind=d2p) :: ZM(3)    ! array (zmin,zmax,zstep)
real(kind=d2p) :: XOI      ! initialization frequency (real part)
real(kind=d2p) :: XC       ! XC  - gyrofrequency (kHZ)
real(kind=d2p) :: PZL      ! option L, =0 linear input, =1 logarythmic input of P and Z

complex(kind=d2p) :: XX(10)
complex(kind=d2p) :: XP(10) 

end module comin




