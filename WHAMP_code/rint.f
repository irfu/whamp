**************************************************************
*
* File rint.f 


      SUBROUTINE RINT(YY,AL,RC)
      IMPLICIT REAL*8 ( A - H, O - Z )
C                  ******** NUMERICAL INTEGRATION ********
      COMPLEX*16 RC(2,2),Y,YY,COT,D,EXF,F,H,O,P,R,RY,RP,RPY,S
      DIMENSION A(16), W(16)
C       ABSCISSAS FOR GAUSSIAN INTEGRATION
      DATA A/               -.98940 09349 91649,-.94457 50230 73232,
     D -.86563 12023 87831, -.75540 44083 55003,-.61787 62444 02643,
     D -.45801 67776 57227, -.28160 35507 79258,-.09501 25098 37637,
     D                       .98940 09349 91649, .94457 50230 73232,
     D  .86563 12023 87831,  .75540 44083 55003, .61787 62444 02643,
     D  .45801 67776 57227,  .28160 35507 79258, .09501 25098 37637/,
     D   W /                 .02715 24594 11754, .06225 35239 38647,
     D  .09515 85116 82492,  .12462 89712 55533, .14959 59888 16576,
     D  .16915 65193 95002,  .18260 34150 44923, .18945 06104 55068,
     D                       .02715 24594 11754, .06225 35239 38647,
     D  .09515 85116 82492,  .12462 89712 55533, .14959 59888 16576,
     D  .16915 65193 95002,  .18260 34150 44923, .18945 06104 55068/,
     D PI/3.14159265358979/
C
      CALL ZEROC2( RC, 1, 2, 1, 2 )
      IF(REAL(YY).LT.0.) THEN
         Y=-YY
         SIG=-1.
      ELSE
         Y=YY
         SIG=1.
      END IF
      YA=DIMAG(Y)
      YR=Y
      UL=PI-2.8*Y/(36.+Y)
      COT=COS(PI*Y)/SIN(PI*Y)
      D=PI*(1.+COT**2)
      C=YR/AL
      XO=LOG(C+SQRT(1.+C**2))
C
      DO 10 I=1,16
      X=UL/2.*(1.+A(I))
      Z=SIN(X)
      C=COS(X)
      G=YR/AL*X/Z
      T=SQRT(1.+G**2)
      B=LOG(G+T)
      G=(1./X-C/Z)*G/T
      T=AL*(T*C-1.)
      Z=EXP(X*YA)
      C=.5*(Z+1./Z)
      S=(0.d0,.5d0)*(Z-1./Z)
      F=COT+G
      H=1.-G*COT
      EXF=EXP(T-Y*B)
      O=B*C+X*S
      P=X*C-B*S
      XY=X*YR
      R=(F*C+H*S)*EXF
      RY=(F*O-H*P+D*(C-G*S))*EXF
      RP=((F*T-H*XY)*C+(H*T+F*XY)*S)*EXF
      RPY=(F*(T*O-XY*P)-H*(T*P+XY*O)+((T+XY*G)*C-(G*T-XY)*S)*D)*EXF
C
      X=XO/2.*(1.+A(I))
      Z=EXP(X)
      C=(Z+1./Z)/2.-1.
      P=EXP(AL*C-Y*X)
      RC(1,1)=RC(1,1)+W(I)*(UL*R+XO*P)
      RC(2,1)=RC(2,1)-W(I)*(UL*RY+XO*X*P)
      RC(1,2)=RC(1,2)+W(I)*(UL*RP+XO*AL*C*P)
      RC(2,2)=RC(2,2)-W(I)*(UL*RPY+XO*AL*X*C*P)
   10 CONTINUE
C
      O=Y/AL
      P=Y**2/2.
      RC(1,1)=O*(Y*RC(1,1)/2.-1.)*SIG
      RC(2,1)=2.*RC(1,1)+O*(P*RC(2,1)+1.)*SIG
      RC(1,2)=Y*O*RC(1,2)/2.*SIG
      RC(2,2)=2.*RC(1,2)+O*P*RC(2,2)*SIG
      END 

      SUBROUTINE ZEROC2( ARRAY, LBD1, UBD1, LBD2, UBD2 )
      IMPLICIT REAL*8 ( A - H, O - Z )
*     .... ARGUMENTS
      INTEGER LBD1, UBD1, LBD2, UBD2
      COMPLEX*16 ARRAY( LBD1 : UBD1,  LBD2 : UBD2 )
*     .... VARIABLES
      INTEGER I1, I2
      DO 20 I2 = LBD2, UBD2
         DO 10 I1 = LBD1, UBD1
            ARRAY( I1, I2 ) = ( 0.0, 0.0 )
   10    CONTINUE
   20 CONTINUE
      END 

