***************************************************************
*
* File output.f 

      SUBROUTINE OUTPT
      include 'comin.h'
      include 'comcout.h' 
      COMPLEX*16 andbz
      INTEGER IFILE  
      SAVE
      PARAMETER( IFILE = 5 )
      CHARACTER CC*1
      PARAMETER( CC = '$' )
      CHARACTER UPP_LTRS*26,LOW_LTRS*26
      PARAMETER( UPP_LTRS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' )
      PARAMETER( LOW_LTRS = 'abcdefghijklmnopqrstuvwxyz' )
      CHARACTER IOU*20
      CHARACTER OFIL*20
      CHARACTER IC*1
C
      PI=3.14159265358979
      K=0
    1 K=K+1
      IF(K.GT.KMX) PRINT 6
      IF(K.GT.KMX) RETURN
      IC = IOU( K : K )
      IF( IC .EQ. '/' )THEN
         PRINT 6
      ELSE
        LINDEX = INDEX(UPP_LTRS,IC)
        IF(LINDEX.EQ.0)LINDEX = INDEX(LOW_LTRS,IC)
*              A  B  C  D  E  F  G  H  I  J  K  L    M
         GOTO(10,18, 1,24,16,10,20,781,1, 1, 1, 773, 775,
     A       777, 771,12, 1,26,22,28,275,793,32, 791,795,14) LINDEX
*              N  O    P  Q  R  S  T  U   V  W    X  Y  Z
      END IF
    6 FORMAT('  ')
      GOTO 1
C
C ** added different stuff by Andris **
c     * empty *
 771  continue
      WRITE(*,772) andE
 772  FORMAT ( ' elip= ',e8.2,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
c ** bp/bz amplitude **
 773  WRITE(*,774) SQRT(dimag(BFL(1))*dimag(BFL(1))+real(BFL(2))
     A *real(BFL(2)))/ABS(BFL(3))
 774  FORMAT( ' bp/bz= ', e8.2, ' ', $)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
c ** bx/by amplitude **
 775  WRITE(*,776) ABS(dimag(BFL(1))/real(BFL(2)))
 776  FORMAT( ' bx/by= ', e8.2,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
c ** ellipticity min(bx,by)/max(bx.by) **
 777  WRITE(*,778) MIN(ABS(dimag(BFL(1))),ABS(real(BFL(2))))
     A /MAX(abs(dimag(BFL(1))),ABS(real(BFL(2))))*(real(BFL(1))*
     A dimag(BFL(2))-dimag(BFL(1))*real(BFL(2)))/ABS(real(BFL(1))*
     A dimag(BFL(2))-dimag(BFL(1))*real(BFL(2)))
 778  FORMAT ( ' e= ',e8.2,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
c ** bz-bx phase angle **      
 791  andbz=BFL(3)*conjg(BFL(1))/ABS(BFL(1))
      andphi=acos(real(andbz)/max(abs(real(andbz)),abs(andbz)))*180/PI
      andphi=sign(andphi,dimag(andbz))
      WRITE(*,792) andphi
 792  FORMAT ( ' phi(bz-bx)= ',f8.2,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
 793  continue
c     write Poynting vector uW/m^2
      coef_poynt = 10.0/4.0/PI/2.0
      Sx = real(efl(2)*conjg(bfl(3))-efl(3)*conjg(bfl(2)))*coef_poynt
      Sy = real(efl(3)*conjg(bfl(1))-efl(1)*conjg(bfl(3)))*coef_poynt
      Sz = real(efl(1)*conjg(bfl(2))-efl(2)*conjg(bfl(1)))*coef_poynt
      write(*,794) Sx,Sy,Sz
 794  FORMAT (' Sx= ',e8.2, ' Sy= ',e8.2, ' Sz= ',e9.3,' ',$) 
c     * abs(e)/abs(b) *
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
 795  continue
      CALL AV
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
 781  dle=EFL(1)*CONJG(EFL(1))+EFL(2)*CONJG(EFL(2))+EFL(3)*CONJG(EFL(3))    
      dlb=BFL(1)*CONJG(BFL(1))+BFL(2)*CONJG(BFL(2))+BFL(3)*CONJG(BFL(3))
      dla=sqrt(dle/dlb)
      WRITE(*,782)  dla
  782  FORMAT(' E/B= ',E9.3,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
 10   WRITE(*, 11) X
   11 FORMAT( ' ',1pE14.7,1PE10.2,'  ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
   12 WRITE(*,13)P
   13 FORMAT(' ',F12.7,'  ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
   14 WRITE(*,15)Z
   15 FORMAT(' ',F12.7,'  ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
   16 WRITE(*,17) EFL
   17 FORMAT(' EX=',F7.4,F8.4,'  EY=',F7.4,F8.4,
     A  '  EZ=',F7.4,F8.4,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
   18 WRITE(*,19) BFL
 19   FORMAT(' BX=',1PE10.2,1PE10.2,'  BY=',1PE10.2,1PE10.2,
     A  '  BZ=',1PE10.2,1PE10.2,' ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
   20 WRITE(*,21) VG
   21 FORMAT(' VGP= ',1PE9.2,'  VGZ= ',1PE9.2,'  ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
   22 WRITE(*,23) SG
c   22 sabs=DIMAG(X)/SQRT(VG(1)*VG(1)+VG(2)*VG(2))*1.e6
c      WRITE(*,231) sabs
   23 FORMAT(' SGP= ',1PE9.2,'  SGZ= ',E9.2,'  ',$)
c  231 FORMAT(' ',1PE9.2,'  ',$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
   24 WRITE(*,25) D, DX, DZ, DP
   25 FORMAT(1P,'D=',2E10.2,'  DX=',2E10.2,'  DZ=',2E10.2,
     F ' DP=',2E10.2,/,$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
   26 WRITE(*,27) RI
   27 FORMAT(' RI=',1P,2E10.2,$)
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 6
  275 WRITE(*,276) ENE*2.0
  276 FORMAT(' ene= ',1PE10.2,$)
C      PRINT 6
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
   28 DO 30 J=1,6
      N=1+J/4+J/6
      M=J-J/4*2-J/6
      WRITE(*,29)
     #(N,M,E(J,I), I=1,4)
   29 FORMAT(' E',2I1,'=',2(1PE10.2),'  EX',2I1,'=',2(1PE10.2),
     F  'EZ',2I1,'=',2(1PE10.2),'  EP',2I1,'=',2(1PE10.2),/)
   30 CONTINUE
      PRINT 6
      IF(IC.NE.'A' .and. ic .ne. 'a') GOTO 1
      PRINT 31, XX, PP ,ZZ
   31 FORMAT( 1P, ' XX=',12E12.3/' PP=',6E12.3/' ZZ=',6E12.3/)
      GOTO 1
C
   32 CALL WRFI(0)
      IF(IOU(1:5).EQ.'W    ')RETURN
      GOTO 1
C
      ENTRY INOUT
c  101 WRITE(*,'(A,$)')'#OUTPUT: '
  101 WRITE(*,'(A)')'#OUTPUT: ' 
      print*
      READ( *, '( A )', IOSTAT = IOS )  IOU
      IF( IOS .NE. 0 )THEN
         IOU = ' '
      END IF
      DO 103 K=1,20
      IC = IOU( K : K )
      IF(IC.EQ.' ') GOTO 103
      IF(IC.EQ.'H') GOTO 104
      KMX=K
      IF(IC.EQ.'W'.OR.IC.EQ.'w')CALL WRFI(2)
  103 CONTINUE
      RETURN
C
  104 PRINT 105
  105 FORMAT( ' THE OUTPUT IS DETERMINED BY A STRING OF LETTERS:'//
     F  ' A     ALL AVAILABLE OUTPUT.'/
     F  ' B     WAVE MAGNETIC FIELD COMPONENTS.'/
     F  ' D     DISPERSION FUNCTION AND DERIVATIVES.'/
     F  ' E     WAVE ELECTRIC FIELD COMPONENTS.'/
     F  ' F     FREQUENCY.'/
     F  ' G     GROUP VELOCITY COMPONENTS.'/
     F  ' P     PERPENDICULAR COMPONENT OF WAVE VECTOR.'/
     F  ' R     REFRACTIVE INDEX.'/
     F  ' S     SPATIAL GROWTH-RATES.'/
     F  ' T     DIELECTRIC TENSOR AND DERIVATIVES.'/
     F  ' Z     Z-COMPONENT OF WAVE VECTOR.'//
     F  ' THE RESULTS ARE  NORMALLY PRINTED ON ONE LINE IN THE ORDER'/
     F  ' THEY ARE SPECIFIED. A NEW LINE IS OBTAINED BY INSERTING'/
     F  ' A "/" IN THE STRING.'/
     F  ' EXAMPLE: OUTPUT: PZF/E'/
     F  ' THE WAVE NUMBERS AND THE FREQUENCY ARE PRINTED ON ONE LINE,'/
     F  ' AND THE ELECTRIC FIELD COMPONENTS ON THE NEXT.'//)
      GOTO 101
      END 




