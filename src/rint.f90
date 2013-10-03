!*************************************************************
!
! File rint.f 


SUBROUTINE RINT(YY,AL,RC)
  implicit none
  ! find the kind of a high precision variable, by finding 
  ! the kind of 1.0d0
  integer, parameter :: d2p=kind(1.0d0)

  !                  ******** NUMERICAL INTEGRATION ********
  integer :: I
  real(kind=d2p) :: AL,B,C,G,SIG,T,UL,Z,X,XO,XY,YR
  complex(kind=d2p) :: YA
  complex(kind=d2p) :: RC(2,2),Y,YY,COT,D,EXF,F,H,O,P,R,RY,RP,RPY,S
  !      real(kind=d2p) :: A(16), W(16)
  !       ABSCISSAS FOR GAUSSIAN INTEGRATION
  real(kind=d2p),dimension(16) :: A=(/  -.989400934991649_d2p,-.944575023073232_d2p,&
       & -.865631202387831_d2p, -.755404408355003_d2p,-.617876244402643_d2p,&
       & -.458016777657227_d2p, -.281603550779258_d2p,-.095012509837637_d2p,&
       &                       .989400934991649_d2p, .944575023073232_d2p,&
       &  .865631202387831_d2p,  .755404408355003_d2p, .617876244402643_d2p,&
       &  .458016777657227_d2p,  .281603550779258_d2p, .095012509837637_d2p/),&
       &   W =(/                 .027152459411754_d2p, .062253523938647_d2p,&
       &  .095158511682492_d2p,  .124628971255533_d2p, .149595988816576_d2p,&
       &  .169156519395002_d2p,  .182603415044923_d2p, .189450610455068_d2p,&
       &                       .027152459411754_d2p, .062253523938647_d2p,&
       &  .095158511682492_d2p,  .124628971255533_d2p, .149595988816576_d2p,&
       &  .169156519395002_d2p,  .182603415044923_d2p, .189450610455068_d2p/)
  real(kind=d2p), parameter :: PI = 3.14159265358979d0
  !
  RC(1:2,1:2)=(0.0d0,0.0d0)

  IF(REAL(YY).LT.0.) THEN
     Y=-YY
     SIG=-1.
  ELSE
     Y=YY
     SIG=1.
  END IF
  YA=DIMAG(Y)
  YR=real(Y)
  UL=PI-REAL(2.8d0*Y/(36.d0+Y))
  COT=COS(PI*Y)/SIN(PI*Y)
  D=PI*(1.d0+COT**2)
  C=YR/AL
  XO=LOG(C+SQRT(1.d0+C**2))
  !
  DO  I=1,16
     X=UL/2.d0*(1.+A(I))
     Z=SIN(X)
     C=COS(X)
     G=YR/AL*X/Z
     T=SQRT(1.d0+G**2)
     B=LOG(G+T)
     G=(1.d0/X-C/Z)*G/T
     T=AL*(T*C-1.d0)
     Z=REAL(EXP(X*YA))
     C=.5*(Z+1.d0/Z)
     S=(0.d0,.5d0)*(Z-1.d0/Z)
     F=COT+G
     H=1.d0-G*COT
     EXF=EXP(T-Y*B)
     O=B*C+X*S
     P=X*C-B*S
     XY=X*YR
     R=(F*C+H*S)*EXF
     RY=(F*O-H*P+D*(C-G*S))*EXF
     RP=((F*T-H*XY)*C+(H*T+F*XY)*S)*EXF
     RPY=(F*(T*O-XY*P)-H*(T*P+XY*O)+((T+XY*G)*C-(G*T-XY)*S)*D)*EXF
     !
     X=XO/2.d0*(1.d0+A(I))
     Z=EXP(X)
     C=(Z+1.d0/Z)/2.d0-1.d0
     P=EXP(AL*C-Y*X)
     RC(1,1)=RC(1,1)+W(I)*(UL*R+XO*P)
     RC(2,1)=RC(2,1)-W(I)*(UL*RY+XO*X*P)
     RC(1,2)=RC(1,2)+W(I)*(UL*RP+XO*AL*C*P)
     RC(2,2)=RC(2,2)-W(I)*(UL*RPY+XO*AL*X*C*P)
  end DO
  !
  O=Y/AL
  P=Y**2/2.d0
  RC(1,1)=O*(Y*RC(1,1)/2.d0-1.d0)*SIG
  RC(2,1)=2.d0*RC(1,1)+O*(P*RC(2,1)+1.d0)*SIG
  RC(1,2)=Y*O*RC(1,2)/2.d0*SIG
  RC(2,2)=2.d0*RC(1,2)+O*P*RC(2,2)*SIG

end SUBROUTINE RINT
