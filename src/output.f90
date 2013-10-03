SUBROUTINE OUTPT
use comin
use comcout
implicit none
integer, parameter :: d2p=kind(1.0d0)

real(kind=d2p) :: andE,andphi,coef_poynt,dla,dlb,dle,Sx,Sy,Sz
COMPLEX(kind=d2p) :: andbz
INTEGER :: I,IOS,J,K,KMX,M,N  
SAVE
integer,PARAMETER :: IFILE = 5 
CHARACTER CC*1
PARAMETER( CC = '$' )
CHARACTER IOU*20
CHARACTER IC*1
105  FORMAT( &
    &  ' The output is determined by a string of letters:'//&
    &  ' b     wave magnetic field components.'/&
    &  ' d     dispersion function and derivatives.'/&
    &  ' e     wave electric field components.'/&
    &  ' f     frequency <real,imaginery>.'/&
    &  ' g     group velocity components.'/&
    &  ' h     |e|/|b| [mV/nT].'/&
    &  ' l     |bp|/|bz|.'/&
    &  ' m     Im[bx]/Re[by].'/&
    &  ' n     ellipticity bx/by.'/&
    &  ' o     ellipticity general'/&
    &  ' p     perpendicular component of wave vector.'/&
    &  ' r     refractive index.'/&
    &  ' s     spatial growth-rates.'/&
    &  ' t     dielectric tensor and derivatives.'/&
    &  ' z     z-component of wave vector.'/&
    &  ' u     total wave energy / energy in electric field.'/&
    &  ' v     Poynting flux (in uW/m2 for <E^2>=0.5(mV/m)^2).'/&
    &  ' x     phase of bz against bx (/+/ means bz is in front of bx).'/&
    &  ' y     energy density and flux of each plasma component.'/& 
    &  ' The results are  normally printed on one line in the order'/&
    &  ' they are specified. A new line is obtained by inserting'/&
    &  ' a "/" in the string.'/&
    &  ' Example: output: pzf/e'/&
    &  ' The wave numbers and the frequency are printed on one line,'/&
    &  ' and the electric field components on the next.'//)

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
        select case (IC)
        case('O','o') ! O
            WRITE(*,772) andE
            772        FORMAT ( ' elip= ',e8.2,' ',$)
            ! ** bp/bz amplitude **
        case('L','l') ! L
            WRITE(*,774) SQRT(dimag(BFL(1))*dimag(BFL(1))+real(BFL(2))&
            & *real(BFL(2)))/ABS(BFL(3))
            774        FORMAT( ' bp/bz= ', e8.2, ' ', $)
        case('M','m') !M
            ! ** bx/by amplitude **
            write(*,776) abs(dimag(bfl(1))/real(bfl(2)))
            776        format( ' bx/by= ', e8.2,' ',$)
            ! ** ellipticity min(bx,by)/max(bx.by) **
        case('N','n') !N
            write(*,778) min(abs(dimag(bfl(1))),abs(real(bfl(2))))&
            & /MAX(abs(dimag(BFL(1))),ABS(real(BFL(2))))*(real(BFL(1))*&
            & dimag(BFL(2))-dimag(BFL(1))*real(BFL(2)))/ABS(real(BFL(1))*&
            & dimag(BFL(2))-dimag(BFL(1))*real(BFL(2)))
            778        FORMAT ( ' e= ',e8.2,' ',$)
            ! ** bz-bx phase angle **     
        case('X','x') ! X
            andbz=BFL(3)*conjg(BFL(1))/ABS(BFL(1))
            andphi=acos(real(andbz)/max(abs(real(andbz)),abs(andbz)))*180/PI
            andphi=sign(andphi,dimag(andbz))
            WRITE(*,792) andphi
            792        FORMAT ( ' phi(bz-bx)= ',f8.2,' ',$)
        case('V','v') ! V
            !     write Poynting vector uW/m^2
            coef_poynt = 10.0/4.0/PI/2.0
            Sx = real(efl(2)*conjg(bfl(3))-efl(3)*conjg(bfl(2)))*coef_poynt
            Sy = real(efl(3)*conjg(bfl(1))-efl(1)*conjg(bfl(3)))*coef_poynt
            Sz = real(efl(1)*conjg(bfl(2))-efl(2)*conjg(bfl(1)))*coef_poynt
            write(*,794) Sx,Sy,Sz
            794        FORMAT (' Sx= ',e8.2, ' Sy= ',e8.2, ' Sz= ',e9.3,' ',$) 
            !     * abs(e)/abs(b) *
        case('Y','y') ! Y
            CALL AV
        case('H','h') !H
            dle=real(EFL(1)*CONJG(EFL(1))+EFL(2)*CONJG(EFL(2))+EFL(3)*CONJG(EFL(3)))
            dlb=real(BFL(1)*CONJG(BFL(1))+BFL(2)*CONJG(BFL(2))+BFL(3)*CONJG(BFL(3)))
            dla=sqrt(dle/dlb)
            WRITE(*,782)  dla
            782        FORMAT(' E/B= ',E9.3,' ',$)
            PRINT 6
        case('F','f') !A,F
            WRITE(*, 11) X
            11         FORMAT( ' ',1pE14.7,1PE10.2,'  ',$)
        case('P','p') ! P
            WRITE(*,13)P
            13         FORMAT(' ',F12.7,'  ',$)
        case('Z','z') ! Z
            WRITE(*,15)Z
            15         FORMAT(' ',F12.7,'  ',$)
        case('E','e') ! E
            WRITE(*,17) EFL
            17         FORMAT(' EX=',F7.4,F8.4,'  EY=',F7.4,F8.4,&
            &  '  EZ=',F7.4,F8.4,' ',$)
        case('B','b') !B
            WRITE(*,19) BFL
            19         FORMAT(' BX=',1PE10.2,1PE10.2,'  BY=',1PE10.2,1PE10.2,&
            &  '  BZ=',1PE10.2,1PE10.2,' ',$)
        case('G','g') ! G
            WRITE(*,21) VG
            21         FORMAT(' VGP= ',1PE9.2,'  VGZ= ',1PE9.2,'  ',$)
        case('S','s') !s
            WRITE(*,23) SG
            !   22 sabs=DIMAG(X)/SQRT(VG(1)*VG(1)+VG(2)*VG(2))*1.e6
            !      WRITE(*,231) sabs
            23         FORMAT(' SGP= ',1PE9.2,'  SGZ= ',E9.2,'  ',$)
            !  231 FORMAT(' ',1PE9.2,'  ',$)
        case('D','d') ! D
            WRITE(*,25) D, DX, DZ, DP
            25         FORMAT(1P,'D=',2E10.2,'  DX=',2E10.2,'  DZ=',2E10.2,&
            & ' DP=',2E10.2,/,$)
        case('R','r') !R
            WRITE(*,27) RI
            27         FORMAT(' RI=',1P,2E10.2,$)
            PRINT 6
        case('U','u') !U
            WRITE(*,276) ENE*2.0
            276        FORMAT(' ene= ',1PE10.2,$)
            !      PRINT 6
        case('T','t') ! T
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
        case(' ') ! ignore spaces
        case default
            write(*,*) 'Unknown output request:',IC
        end select
    endif
end do output_loop
!
ENTRY INOUT
output_format_loop: do 
    WRITE(*,'(A)')'#OUTPUT: ' 
    print*
    READ( *, '( A )', IOSTAT = IOS )  IOU
    IF( IOS .NE. 0 )THEN
        IOU = ' '
    END IF
    DO K=1,20 ! loop through all input characters
        IC = IOU( K : K )
        select case (IC)
        case('H','h','?') ! help
            PRINT 105 ! print output options, defined in the beginning of file
            cycle output_format_loop  
        case(' ')
        case ('+') ! print debug info
            printDebugInfo = .true.
            IOU(K:K)=' '
        case ('-') ! remove printing debug info
            printDebugInfo = .false.
            IOU(K:K)=' '
        case default ! any character
            KMX=K
        end select
    end do
    exit output_format_loop
end do output_format_loop
end SUBROUTINE OUTPT
