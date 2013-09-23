!**************************************************************
!
! File output.f 

SUBROUTINE OUTPT
  use comin
  use comcout
  implicit none
  ! find the kind of a high precision variable, by finding 
  ! the kind of 1.0d0
  integer, parameter :: d2p=kind(1.0d0)

  real(kind=d2p) :: andE,andphi,coef_poynt,dla,dlb,dle,Sx,Sy,Sz
  COMPLEX(kind=d2p) :: andbz
  INTEGER :: I,IOS,J,K,KMX,LINDEX,M,N  
  SAVE
  integer,PARAMETER :: IFILE = 5 
  CHARACTER CC*1
  PARAMETER( CC = '$' )
  CHARACTER UPP_LTRS*26,LOW_LTRS*26
  PARAMETER( UPP_LTRS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' )
  PARAMETER( LOW_LTRS = 'abcdefghijklmnopqrstuvwxyz' )
  CHARACTER IOU*20
  CHARACTER OFIL*20
  CHARACTER IC*1
  !
  real(kind=d2p),parameter :: PI=3.14159265358979_d2p
  K=0
  output_loop: do 
     K=K+1
     IF(K.GT.KMX) PRINT 6
     IF(K.GT.KMX) RETURN
     IC = IOU( K : K )
     IF( IC .EQ. '/' )THEN
        PRINT 6
6       FORMAT('  ')
     ELSE
        LINDEX = INDEX(UPP_LTRS,IC)
        IF(LINDEX.EQ.0) LINDEX = INDEX(LOW_LTRS,IC)
        !              A  B  C  D  E  F  G  H  I  J  K  L    M
        !         GOTO(10,18, 1,24,16,10,20,781,1, 1, 1, 773, 775,
        !     A       777, 771,12, 1,26,22,28,275,793,32, 791,795,14) LINDEX
        !              N  O    P  Q  R  S  T  U   V  W    X  Y  Z
        !
        select case(LINDEX)
           ! ** added different stuff by Andris **
           !     * empty *
        case(15) ! O
           WRITE(*,772) andE
772        FORMAT ( ' elip= ',e8.2,' ',$)
           ! ** bp/bz amplitude **
        case(12) ! L
           WRITE(*,774) SQRT(dimag(BFL(1))*dimag(BFL(1))+real(BFL(2))&
                & *real(BFL(2)))/ABS(BFL(3))
774        FORMAT( ' bp/bz= ', e8.2, ' ', $)
        case(13) !M
           ! ** bx/by amplitude **
           write(*,776) abs(dimag(bfl(1))/real(bfl(2)))
776        format( ' bx/by= ', e8.2,' ',$)
           ! ** ellipticity min(bx,by)/max(bx.by) **
        case(14) !N
           write(*,778) min(abs(dimag(bfl(1))),abs(real(bfl(2))))&
                & /MAX(abs(dimag(BFL(1))),ABS(real(BFL(2))))*(real(BFL(1))*&
                & dimag(BFL(2))-dimag(BFL(1))*real(BFL(2)))/ABS(real(BFL(1))*&
                & dimag(BFL(2))-dimag(BFL(1))*real(BFL(2)))
778        FORMAT ( ' e= ',e8.2,' ',$)
           ! ** bz-bx phase angle **     
        case(24) ! X
           andbz=BFL(3)*conjg(BFL(1))/ABS(BFL(1))
           andphi=acos(real(andbz)/max(abs(real(andbz)),abs(andbz)))*180/PI
           andphi=sign(andphi,dimag(andbz))
           WRITE(*,792) andphi
792        FORMAT ( ' phi(bz-bx)= ',f8.2,' ',$)
        case(22) ! V
           !     write Poynting vector uW/m^2
           coef_poynt = 10.0/4.0/PI/2.0
           Sx = real(efl(2)*conjg(bfl(3))-efl(3)*conjg(bfl(2)))*coef_poynt
           Sy = real(efl(3)*conjg(bfl(1))-efl(1)*conjg(bfl(3)))*coef_poynt
           Sz = real(efl(1)*conjg(bfl(2))-efl(2)*conjg(bfl(1)))*coef_poynt
           write(*,794) Sx,Sy,Sz
794        FORMAT (' Sx= ',e8.2, ' Sy= ',e8.2, ' Sz= ',e9.3,' ',$) 
           !     * abs(e)/abs(b) *
        case(25) ! Y
           CALL AV
        case(8) !H
           dle=EFL(1)*CONJG(EFL(1))+EFL(2)*CONJG(EFL(2))+EFL(3)*CONJG(EFL(3))    
           dlb=BFL(1)*CONJG(BFL(1))+BFL(2)*CONJG(BFL(2))+BFL(3)*CONJG(BFL(3))
           dla=sqrt(dle/dlb)
           WRITE(*,782)  dla
782        FORMAT(' E/B= ',E9.3,' ',$)
           PRINT 6
        case(1,6) !A,F
           WRITE(*, 11) X
11         FORMAT( ' ',1pE14.7,1PE10.2,'  ',$)
        case(16) ! P
           WRITE(*,13)P
13         FORMAT(' ',F12.7,'  ',$)
        case(26) ! Z
           WRITE(*,15)Z
15         FORMAT(' ',F12.7,'  ',$)
        case(5) ! E
           WRITE(*,17) EFL
17         FORMAT(' EX=',F7.4,F8.4,'  EY=',F7.4,F8.4,&
                &  '  EZ=',F7.4,F8.4,' ',$)
        case(2) !B
           WRITE(*,19) BFL
19         FORMAT(' BX=',1PE10.2,1PE10.2,'  BY=',1PE10.2,1PE10.2,&
                &  '  BZ=',1PE10.2,1PE10.2,' ',$)
        case(7) ! G
           WRITE(*,21) VG
21         FORMAT(' VGP= ',1PE9.2,'  VGZ= ',1PE9.2,'  ',$)
        case(19) !s
           WRITE(*,23) SG
           !   22 sabs=DIMAG(X)/SQRT(VG(1)*VG(1)+VG(2)*VG(2))*1.e6
           !      WRITE(*,231) sabs
23         FORMAT(' SGP= ',1PE9.2,'  SGZ= ',E9.2,'  ',$)
           !  231 FORMAT(' ',1PE9.2,'  ',$)
        case(4) ! D
           WRITE(*,25) D, DX, DZ, DP
25         FORMAT(1P,'D=',2E10.2,'  DX=',2E10.2,'  DZ=',2E10.2,&
                & ' DP=',2E10.2,/,$)
        case(18) !R
           WRITE(*,27) RI
27         FORMAT(' RI=',1P,2E10.2,$)
           PRINT 6
        case(21) !U
           WRITE(*,276) ENE*2.0
276        FORMAT(' ene= ',1PE10.2,$)
           !      PRINT 6
        case(20) ! T
           DO 30 J=1,6
              N=1+J/4+J/6
              M=J-J/4*2-J/6
              WRITE(*,29) (N,M,E(J,I), I=1,4)
29            FORMAT(' E',2I1,'=',2(1PE10.2),'  EX',2I1,'=',2(1PE10.2),&
                   &  'EZ',2I1,'=',2(1PE10.2),'  EP',2I1,'=',2(1PE10.2),/)
30         end do
           PRINT 6
           PRINT 31, XX, PP ,ZZ
31         FORMAT( 1P, ' XX=',12E12.3/' PP=',6E12.3/' ZZ=',6E12.3/)
        case(3,9:11,17)
           ! do nothing
        case(23) ! W
           CALL WRFI(0)
           IF(IOU(1:5).EQ.'W    ')RETURN
        end select
     end if
  end do output_loop
  !
  ENTRY INOUT
  !  101 WRITE(*,'(A,$)')'#OUTPUT: '
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
103  CONTINUE
     RETURN
     !
104  PRINT 105
105  FORMAT( ' THE OUTPUT IS DETERMINED BY A STRING OF LETTERS:'//&
          &  ' A     ALL AVAILABLE OUTPUT.'/&
          &  ' B     WAVE MAGNETIC FIELD COMPONENTS.'/&
          &  ' D     DISPERSION FUNCTION AND DERIVATIVES.'/&
          &  ' E     WAVE ELECTRIC FIELD COMPONENTS.'/&
          &  ' F     FREQUENCY.'/&
          &  ' G     GROUP VELOCITY COMPONENTS.'/&
          &  ' P     PERPENDICULAR COMPONENT OF WAVE VECTOR.'/&
          &  ' R     REFRACTIVE INDEX.'/&
          &  ' S     SPATIAL GROWTH-RATES.'/&
          &  ' T     DIELECTRIC TENSOR AND DERIVATIVES.'/&
          &  ' Z     Z-COMPONENT OF WAVE VECTOR.'//&
          &  ' THE RESULTS ARE  NORMALLY PRINTED ON ONE LINE IN THE ORDER'/&
          &  ' THEY ARE SPECIFIED. A NEW LINE IS OBTAINED BY INSERTING'/&
          &  ' A "/" IN THE STRING.'/&
          &  ' EXAMPLE: OUTPUT: PZF/E'/&
          &  ' THE WAVE NUMBERS AND THE FREQUENCY ARE PRINTED ON ONE LINE,'/&
          &  ' AND THE ELECTRIC FIELD COMPONENTS ON THE NEXT.'//)
     GOTO 101
end SUBROUTINE OUTPT




