! Assymptotic series estimate
SUBROUTINE RASY(Y, AL, RC)
   implicit none
   integer, parameter :: d2p = kind(1.0d0)
   integer :: M, N
   real(kind=d2p) :: A, AL, AY, C, T
   COMPLEX(kind=d2p) :: Y, Y2, COT, P, PY, PP, PPY, PN, PYN, QN, QYN, RC(2, 2)
   real(kind=d2p), parameter :: PI = 3.14159265358979_d2p

   Y2 = Y*Y
   COT = COS(PI*Y)/SIN(PI*Y)
   !                  1.E99 IS TOO BIG FOR S/370 HARDWARE. SET TO LARGEST
   !                  POSSIBLE FOR IBM MACHINES
   !     C=1.E99
   C = 7.2d35
   PN = -Y/AL
   PYN = PN
   A = 1./(AL*SQRT(2.*PI*AL))
   QN = PI*Y2*COT*A
   QYN = QN*(2.-Y*PI*COT) - Y*PI**2*Y2*A
   !
   P = PN + QN
   PY = PYN + QYN
   PP = -PN - 1.5*QN
   PPY = -PYN - 1.5*QYN
   AY = ABS(Y) + 2.
   !
   iteration: DO N = 1, 100
      M = N - 1
      PYN = (PYN*(M*M - Y2) - 2.*Y2*PN)/((2*M + 1)*AL)
      PN = PN*(M*M - Y2)/((2*M + 1)*AL)
      QYN = (QYN*((M + .5)**2 - Y2) - 2.*Y2*QN)/(2.*N*AL)
      QN = QN*((M + .5)**2 - Y2)/(2.*N*AL)
      IF (M .LT. AY) then
      else
         C = N*(ABS(PN) + ABS(QN))
         IF (C .LE. 1.E-7*ABS(PP)) exit iteration
         IF (C .GE. T) exit iteration
      end if
      P = P + PN + QN
      PY = PY + PYN + QYN
      PP = PP - (N + 1.)*PN - (N + 1.5)*QN
      PPY = PPY - (N + 1.)*PYN - (N + 1.5)*QYN
      T = C
   end DO iteration
   !
   RC(1, 1) = P + PN + QN
   RC(2, 1) = PY + PYN + QYN
   RC(1, 2) = PP + P
   RC(2, 2) = PPY + PY
end SUBROUTINE RASY
