!**************************************************************
!
! File comin.h 

 module comin
        implicit none
! find the kind of a high precision variable, by finding 
! the kind of 1.0d0
        integer, parameter,private :: d2p=kind(1.0d0)

!common block plasma
        integer :: JMA
        real(kind=d2p) :: BVEC(3),DBVDR(3,4),DBDR(4),DNDR(4,10),DTDR(4,10),&
        & VDRIFT(10),ZIGN 
        complex(kind=d2p) :: GE(6,4) 

!common block XPZ
        real(kind=d2p) :: PP(10),ZZ(10),AA(10,2),DD(10),ASS(10),VD(10),&
       & DN(10),TA(10),CV,PM(3),ZM(3),XOI,XC,PZL 
        complex(kind=d2p) :: XX(10),XP(10) 

 end module comin

! CV  - speed of light / thermal velocity of 1st species
! JMA - the number of last plasma species with non-zero density
! PM  - array (pmin,pmax,pstep)
! PZL - option L, =0 linear input, =1 logarythmic input of P and Z
! ZM  - arra (zmin,zmax,zstep)
! XC  - gyrofrequency (kHZ)
! XOI - initialization frequency (real part)
