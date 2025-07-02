SUBROUTINE WRFI(IAR)
   use comin
   use comcout
   implicit none
   integer, parameter :: d2p = kind(1.0d0)
   SAVE
   integer :: I, IAR, ISP, J, K, flag
   real(kind=d2p) :: EBU, ENEBU, GBU, gf, PBU, ptok, rn, SBU, VGBU, vth, XBU, ZBU
   DIMENSION PBU(50), ZBU(50), XBU(50), GBU(50), SBU(50), VGBU(2, 50),&
        & ENEBU(20), EBU(2, 3, 50)
   CHARACTER SPE(5)*3, OFIL*16, GR*3, WDFV*1, SWDF*1
   integer  :: IFS = 0, IB = 0
   DATA SPE/'E- ', 'H+ ', 'HE+', 'O+ ', '   '/
   !
   flag = IAR + 1
   if (flag == 1) then
      IF (IFS .EQ. 0) RETURN
      IB = IB + 1
      !      PBU(IB)=P
      !      ZBU(IB)=Z
      !      XBU(IB)=REAL(X)
      pbu(ib) = p*ptok
      zbu(ib) = z*ptok
      xbu(ib) = x*gf
      GBU(IB) = DIMAG(X)
      SBU(IB) = GBU(IB)/SQRT(VG(1)**2 + VG(2)**2 + 1.d-66)
      DO I = 1, 3
         EBU(1, I, IB) = REAL(EFL(I))
         EBU(2, I, IB) = DIMAG(EFL(I))
      end DO
      !      VGBU(1,IB)=VG(1)
      !      VGBU(2,IB)=VG(2)

      vgbu(1, ib) = vg(1)*vth
      vgbu(2, ib) = vg(2)*vth

      ENEBU(IB) = ENE
      IF (IB .LT. 10) RETURN
      flag = 4
   end if
   if (flag == 4) then
      IF (IB*IFS .LT. 1) RETURN
      IF (GR(1:1) .EQ. 'Y' .or. gr(1:1) .eq. 'y') THEN
         IF (WDFV .EQ. 'Y' .or. wdfv .eq. 'y') THEN
            DO I = 1, IB
               WRITE (10, 12, ERR=19) PBU(I), ZBU(I), XBU(I), GBU(I), SBU(I),&
                     &((EBU(K, J, I), K=1, 2), J=1, 3), ENEBU(I), (VGBU(J, I), J=1, 2)
               IF (SWDF .EQ. 'Y' .or. swdf .eq. 'y') THEN
                  WRITE (11, 16, ERR=19) PBU(I), ZBU(I), XBU(I), GBU(I), SBU(I)
               END IF
            end DO
         ELSE
            WRITE (10, 16, ERR=19) (PBU(I), ZBU(I), XBU(I), GBU(I), SBU(I), I=1, IB)
         END IF
      ELSE
         IF (WDFV .EQ. 'Y' .or. wdfv .eq. 'y') THEN
            DO I = 1, IB
               WRITE (10, 13, ERR=19) PBU(I), ZBU(I), XBU(I),&
                     &      ((EBU(K, J, I), K=1, 2), J=1, 3),&
                     &ENEBU(I), (VGBU(J, I), J=1, 2)
               IF (SWDF .EQ. 'Y' .or. swdf .eq. 'y') THEN
                  WRITE (11, 11, ERR=19) PBU(I), ZBU(I), XBU(I)
               END IF
            end DO
         ELSE
            WRITE (10, 11, ERR=19) (PBU(I), ZBU(I), XBU(I), I=1, IB)
         END IF
      END IF
12    FORMAT(' P=', 1pE13.7, ' Z=', E14.7, ' X=', E13.7, E14.7,&
             &' S=', E14.7, /, 'EX=', 2(E14.7), 'EY=', 2(E14.7), /,&
             &'EZ=', 2(E14.7), /, 'ENE=', E14.7, 'VGP=', E14.7, 'VGZ=', E14.7)
13    FORMAT(' P=', 1pE13.7, ' Z=', E14.7, ' X=', E13.7,&
             &/, 'EX=', 2(E14.7), 'EY=', 2(E14.7), /,&
             &'EZ=', 2(E14.7), /, 'ENE=', E14.7, 'VGP=', E14.7, 'VGZ=', E14.7)
11    FORMAT(' P=', 1pE13.7, ' Z=', E13.7, ' X=', E13.7)
16    FORMAT(' P=', 1pE13.7, ' Z=', E13.7, ' X=', E13.7, E14.7, ' S=', E14.7)
      IB = 0
      RETURN
      !
19    WRITE (*, *) 'WRFI: ERROR IN WRITE, FILE=', OFIL
      RETURN
   end if
   !
   if (flag == 2) then
      IF (IFS == 0) RETURN
      flag = 3
   end if
   !
   if (flag == 3) then
      do
         WRITE (*, '(A,$)') 'OUTPUT FILE:'
         READ (*, '(A)') OFIL
         IB = 0
         IF (OFIL .EQ. '        ') RETURN
         IF (IFS .EQ. 1) CLOSE (6)
         OPEN (10, FILE=OFIL, ERR=40)
         IFS = 1

         rn = mi_o_me*ass(1)
         if (rn .lt. 1.) rn = 1.
         gf = 1000.*xc/rn
         vth = sqrt(3.517388e14*ta(1)/rn)
         ptok = gf/vth

         WRITE (*, '(A,$)') 'GROWTH RATES?'
         READ (*, '(A)') GR
         WRITE (*, '(A,$)') 'VARIABLES FOR WDF RECONSTRUCTION?'
         READ (*, '(A)') WDFV
         SWDF = 'N'
         IF (WDFV .EQ. 'Y' .or. wdfv .eq. 'y') THEN
            WRITE (*, '(A,$)') 'WISH TO USE SEPARATE FILES (Y/N)'
            READ (*, '(A)') SWDF
            IF (SWDF .EQ. 'Y' .or. swdf .eq. 'y') THEN
               WRITE (*, '(A,$)') 'FILENAME?'
               READ (*, '(A)') OFIL
               OPEN (11, FILE=OFIL, ERR=40)
            END IF
         END IF
         !300 FORMAT(2X,A33,/)
         DO J = 1, 6
            ISP = SQRT(ASS(J))
            IF (ISP .LT. 4) ISP = ISP + 1
            WRITE (10, 33) SPE(ISP), DN(J), TA(J), DD(J), AA(J, 1), AA(J, 2), VD(J)
            IF (SWDF .EQ. 'Y') THEN
               WRITE (11, 33) SPE(ISP), DN(J), TA(J), DD(J), AA(J, 1), AA(J, 2), VD(J)
            END IF
         end DO
33       FORMAT(A3, ' N=', E9.3, ' T=', F8.4, ' D=', F4.2, ' A=', F4.2, ' B=',&
                &   F4.2, ' VD=', F5.2)
         RETURN
         !
40       WRITE (*, *) 'WRFI: ERROR IN OPEN, FILE=', OFIL
      end do
   end if
   RETURN
end SUBROUTINE WRFI
