*********************************************************
*
* File av.f
* This subroutine calculates the individual contribution
* of plasma species to the energy density and the energy 
* flux. 
* Energy density units - electric field energy density is 1
* Energy flux units - part of the total energy flux 
*
* Loop through all species 
* For each species 
*
* Variables
*           EE(6,4) - dielectric tensor and its derivatives
*             ENDEN - energy density in particular species
*           ENFLUXZ - energy flux z-direction	in part. species
*           ENFLUXP - energy flux p-direction in part. species
*               JMA - number of species
*           ENFIELD - energy in fields
*              POYN - Poynting flux 

      SUBROUTINE AV
      include 'comin.h'
      include 'comcout.h'
c      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A,B
      COMPLEX*16 XSI(6,4), DF, U1, U2, U3, U12, U13, U32
      COMPLEX*16 ERG, DDER(3,3), EE(6,4)
      COMPLEX*16 POYN(3),ENFLZ,ENFLP,DEP(3,3),DEZ(3,3)

	PI=3.1415926535897
C
C
C *** estimate the number of species JMA
       DO 102 J=1,10
      IF(DN(J).EQ.0.) GOTO 103
      JMA=J
  102 CONTINUE
  103 CONTINUE
C *** calculate Poynting flux POYN in whamp units
      coef=10.0/4.0/PI/2.0
      coef=coef*4/8.8542e-12*CV/2.9979e8
      POYN(1)=real(efl(2)*conjg(bfl(3))-efl(3)*conjg(bfl(2)))*coef
      POYN(2)=real(efl(3)*conjg(bfl(1))-efl(1)*conjg(bfl(3)))*coef
      POYN(3)=real(efl(1)*conjg(bfl(2))-efl(2)*conjg(bfl(1)))*coef

C *** calculate energy in fields
      A=EFL(1)*CONJG(EFL(1))+EFL(2)*CONJG(EFL(2))+EFL(3)*CONJG(EFL(3))    
      B=BFL(1)*CONJG(BFL(1))+BFL(2)*CONJG(BFL(2))+BFL(3)*CONJG(BFL(3))
      ENFIELD=1.0+B/A*299.79*299.79
            
C *** Loop through all species

      DO 112 J=1,JMA
c
C              *********** FORM DIELECTRIC TENSOR ************
      DO 131 K=1,4
      DO 131 I=1,6
  131 E(I,K)=(0.,0.)
      E(1,1)=1.
      E(4,1)=1.
      E(6,1)=1.
C
cav      DO 135 J=1,JMA
      E(1,1)=E(1,1)-XP(J)
      E(4,1)=E(4,1)-XP(J)
      E(6,1)=E(6,1)-XP(J)
      IF(AA(J,1).NE.AA(J,2)) GOTO 132
      AA(J,2)=0.
      DD(J)=1.
  132 IB=1
      DF=XP(J)/(AA(J,1)*(AA(J,1)-AA(J,2)))
      Q=AA(J,1)-DD(J)*AA(J,2)
	IERR=0
  133 CALL CHI(XSI,J,IB,4,IERR)
      IF(IERR.NE.0) RETURN 
      DO 134 K=1,4
      DO 134 I=1,6
      E(I,K)=E(I,K)+DF*Q*XSI(I,K)
  134 CONTINUE
C
      IF(IB.EQ.2) GOTO 135
      Q=(DD(J)-1.)*AA(J,1)
      IF(Q.EQ.0.) GOTO 135
      IB=2
      GOTO 133
  135 CONTINUE
C                       *** DIELECTRIC TENSOR COMPUTED ***
C
C       ******* FORM REFRACTIVE INDEX, CV=SPEED OF LIGHT/THERM. SPEED. *
      U1=PP(1)*CV/XX(1)
      U3=ZZ(1)*CV/XX(1)
      U12=U1*U1
      U32=U3*U3
      U2=U12+U32
      U13=2.*U1*U3
C
C        ****** COMPLETE X-DERIVATIVE OF DIELECTRIC TENSOR ******
      E(1,2)=E(1,2)-2.*(E(1,1)-1.)
      E(2,2)=E(2,2)-2.* E(2,1)
      E(3,2)=E(3,2)-2.* E(3,1)
      E(4,2)=E(4,2)-2.*(E(4,1)-1.)
      E(5,2)=E(5,2)-2.* E(5,1)
      E(6,2)=E(6,2)-2.*(E(6,1)-1.)

	DO 151 I=1,6
      DO 151 K=1,4
  151  EE(I,K)=E(I,K)
C
      DDER(1,1)=EE(1,1)+EE(1,2)+U32
      DDER(1,2)=EE(2,1)+EE(2,2)
      DDER(1,3)=EE(3,1)+EE(3,2)-U1*U3
      DDER(2,2)=EE(4,1)+EE(4,2)+U2
      DDER(2,3)=EE(5,1)+EE(5,2)
      DDER(3,3)=EE(6,1)+EE(6,2)+U12
      DDER(2,1)=-DDER(1,2)
      DDER(3,1)=DDER(1,3)
      DDER(3,2)=-DDER(2,3)
C
      ERG=(0.,0.)
      DO 20 I=1,3
      DO 20 K=1,3
   20 ERG=ERG+ CONJG(EFL(I))*DDER(I,K)*EFL(K)
C
      ENDEN=ERG*CONJG(ERG)/(ERG+CONJG(ERG))*2.0-ENFIELD

C
      DEZ(1,1)=EE(1,3)
      DEZ(1,2)=EE(2,3)
      DEZ(1,3)=EE(3,3)
      DEZ(2,2)=EE(4,3)
      DEZ(2,3)=EE(5,3)
      DEZ(3,3)=EE(6,3)
      DEZ(2,1)=-DEZ(1,2)
      DEZ(3,1)=DEZ(1,3)
      DEZ(3,2)=-DEZ(2,3)

C
      DEP(1,1)=EE(1,4)
      DEP(1,2)=EE(2,4)
      DEP(1,3)=EE(3,4)
      DEP(2,2)=EE(4,4)
      DEP(2,3)=EE(5,4)
      DEP(3,3)=EE(6,4)
      DEP(2,1)=-DEP(1,2)
      DEP(3,1)=DEP(1,3)
      DEP(3,2)=-DEP(2,3)

      ENFLZ=(0.,0.)
      ENFLP=(0.,0.)
      DO 141 I=1,3
      DO 141 K=1,3
      ENFLZ=ENFLZ+ CONJG(EFL(I))*DEZ(I,K)*EFL(K)
      ENFLP=ENFLP+ CONJG(EFL(I))*DEP(I,K)*EFL(K)
  141 CONTINUE

      ENFLUXZ=ENFLZ*CONJG(ENFLZ)/(ENFLZ+CONJG(ENFLZ))*2.0
      ENFLUXP=ENFLP*CONJG(ENFLP)/(ENFLP+CONJG(ENFLP))*2.0

	A=2.9979e8/CV*8.8542e-12/4.0
	ENFLUXZ=-ENFLUXZ*real(X)/Z*A
	ENFLUXP=-ENFLUXP*real(X)/P*A
      write(*,152) ENDEN,ENFLUXP,ENFLUXZ
  152 FORMAT (' enden= ',e9.3,' enfl_p=',e9.3,' enfl_z= ',e9.3,' ',$) 


  112 CONTINUE

      END


