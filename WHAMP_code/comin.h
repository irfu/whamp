***************************************************************
*
* File comin.h 

      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XX,XP,GE 

      COMMON /PLASMA/JMA,BVEC(3),DBVDR(3,4),DBDR(4),
     C DNDR(4,10),DTDR(4,10),GE(6,4),VDRIFT(10),ZIGN  

      COMMON /XPZ/ XX(10),PP(10),ZZ(10),AA(10,2),DD(10),ASS(10),VD(10),
     C  DN(10),TA(10),XP(10),CV,PM(3),ZM(3),XOI,XC,PZL 
   
