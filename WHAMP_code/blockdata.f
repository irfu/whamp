***************************************************************
*
* File blockdata.f 

      BLOCK DATA
      IMPLICIT REAL*8 ( A - H, O - Z )
      COMMON /XPZ/ XX,PP,ZZ,A,B,D,ASS,VD,
     C   DN,TA,XP,CV,PM,ZM,XOI,XC,PZL
      COMPLEX*16 XX( 1 : 10 ), XP( 1 : 10 )
      REAL*8 PP( 1 : 10 ), ZZ( 1 : 10 ), A( 1 : 10 ), B( 1 : 10 )
      REAL*8 D( 1 : 10 )
      REAL*8 ASS( 1 : 10 ), VD( 1 : 10 ), DN( 1 : 10 ), TA( 1 : 10 )
      REAL*8 PM( 1 : 3 ), ZM( 1 : 3 )
      REAL*8 CV, XOI, XC, PZL
       DATA DN/37.5e6,112.5e6,99.9e6,.1e6,0.0,0.0,0.0,0.0,0.0,0.0/,
     #  TA/.001,.001,.001,.125,.015,.030,.080,.2,.24,.001/,
     #   D/1.0,1.0,1.0,1.0,1.0,1.0,1.,1.,1.,1./,
     #   A/1.,1.,1.,1.,1.,1.,1.,1.,1.,1./,
     #   B/0.10,0.10,0.10,0.10,0.10,0.10,.1,.1,.1,.1/,
     # ASS/1.,16.,0.,0.,1.,1.,1.,1.,1.,1./,
     #  VD/0.,0.0,0.0,2,0.0,0.0,0.0,0.0,0.0,0.0/,
     #  XC/777/,
     # PZL/1/
C **** CORRESPONDENCE BETWEEN ARRAY AND NAMES IN WHAMP:  ****
C  ARRAY(IOF+ 0 12 18 24 30 36 42  48 54 60 66 78 79 82 85  86  87
C            XX PP ZZ A  B  D  ASS VD DN TA XP CV PM ZM XOI XC  PZL
C
      END

