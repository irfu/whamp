***************************************************************
*
* File rtay.f 


      SUBROUTINE RTAY(Y,AL,RC)
      IMPLICIT REAL*8 ( A - H, O - Z )
      COMPLEX*16 Y,Y2,RC(2,2),PN,PYN,COT
C                  ******** TAYLOR SERIES ********
      Y2=Y*Y
   10 PN=Y/(Y2-1.)
      PYN=-Y*(Y2+1.)/(Y2-1.)**2
      RC(1,1)=PN
      RC(1,2)=PN
      RC(2,1)=PYN
      RC(2,2)=PYN
C
      DO 1 I=2,100
      COT=(2*I-1)/(Y2-I**2)*AL
      PYN=COT*(PYN-2.*Y2/(Y2-I**2)*PN)
      PN=COT*PN
      RC(1,1)=RC(1,1)+PN
      RC(2,1)=RC(2,1)+PYN
      RC(1,2)=RC(1,2)+I*PN
      RC(2,2)=RC(2,2)+I*PYN
      T=ABS(PN)*1.E8
      IF(T.LT.ABS(RC(1,1))) GOTO 2
    1 CONTINUE
    2 CONTINUE
      END
