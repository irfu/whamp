!**************************************************************
!
! File blockdata.f 

module default_values
      implicit none
! find the kind of a high precision variable, by finding 
! the kind of 1.0d0
        integer, parameter :: d2p=kind(1.0d0)
      
!      COMMON /XPZ/ XX,PP,ZZ,A,B,D,ASS,VD,
!     C   DN,TA,XP,CV,PM,ZM,XOI,XC,PZL
      COMPLEX(kind=d2p) :: XX(1:10), XP(1:10)
      REAL(kind=d2p) PP( 1 : 10 ), ZZ( 1 : 10 ) 
      REAL(kind=d2p) PM( 1 : 3 ), ZM( 1 : 3 )
      REAL(kind=d2p) CV, XOI
      real(kind=d2p),dimension(10) ::&
        &  DN=[37.5e6,112.5e6,99.9e6,.1e6,0.0,0.0,0.0,0.0,0.0,0.0],&
     &  TA=[.001,.001,.001,.125,.015,.030,.080,.2,.24,.001],&
     &   D=[1.0,1.0,1.0,1.0,1.0,1.0,1.,1.,1.,1.],&
     &   A=[1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],&
     &   B=[0.10,0.10,0.10,0.10,0.10,0.10,.1,.1,.1,.1],&
     & ASS=[1.,16.,0.,0.,1.,1.,1.,1.,1.,1.],&
     &  VD=[0.,0.0,0.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0]
        real(kind=d2p) :: XC=777
        integer :: PZL=1
! **** CORRESPONDENCE BETWEEN ARRAY AND NAMES IN WHAMP:  ****
!  ARRAY(IOF+ 0 12 18 24 30 36 42  48 54 60 66 78 79 82 85  86  87
!            XX PP ZZ A  B  D  ASS VD DN TA XP CV PM ZM XOI XC  PZL
!
end module default_values

