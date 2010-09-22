***************************************************************
*
* File wrfi.f 

      SUBROUTINE WRFI(IAR)
      include 'comin.h'
      include 'comcout.h'
      SAVE
      DIMENSION PBU(50),ZBU(50),XBU(50),GBU(50),SBU(50),VGBU(2,50),
     # ENEBU(20), EBU(2,3,50)
      CHARACTER SPE(5)*3,MOD*5,OFIL*16,GR*3,WDFV*1,SWDF*1
      DATA IFS/0/,IB/0/,SPE/'E- ','H+ ','HE+','O+ ','   '/
C
      GOTO (1,20,30,10) IAR+1
    1 IF(IFS.EQ.0) RETURN
      IB=IB+1
c      PBU(IB)=P
c      ZBU(IB)=Z
c      XBU(IB)=REAL(X)
         pbu(ib) = p*ptok
         zbu(ib) = z*ptok
         xbu(ib) = x*gf
      GBU(IB)=DIMAG(X)
      SBU(IB)=GBU(IB)/SQRT(VG(1)**2+VG(2)**2+1.E-66)
      DO 2 I=1,3
        EBU(1,I,IB)=REAL(EFL(I))
        EBU(2,I,IB)=DIMAG(EFL(I))
   2  CONTINUE
c      VGBU(1,IB)=VG(1)
c      VGBU(2,IB)=VG(2)

          vgbu(1,ib) = vg(1)*vth
          vgbu(2,ib) = vg(2)*vth

      ENEBU(IB)=ENE
      IF(IB.LT.10) RETURN
   10 IF(IB*IFS.LT.1) RETURN
      IF(GR(1:1).EQ.'Y' .or. gr(1:1) .eq. 'y')THEN
        IF(WDFV.EQ.'Y' .or. wdfv .eq. 'y')THEN
          DO 21 I=1,IB
            WRITE(10,12,ERR=19)PBU(I),ZBU(I),XBU(I),GBU(I),SBU(I),
     #((EBU(K,J,I),K=1,2),J=1,3),ENEBU(I),(VGBU(J,I),J=1,2)
            IF(SWDF.EQ.'Y' .or. swdf .eq. 'y')THEN
            WRITE(11,16,ERR=19)PBU(I),ZBU(I),XBU(I),GBU(I),SBU(I)
            ENDIF
   21     CONTINUE
        ELSE
          WRITE(10,16,ERR=19)(PBU(I),ZBU(I),XBU(I),GBU(I),SBU(I),I=1,IB)
        ENDIF
      ELSE
        IF(WDFV.EQ.'Y' .or. wdfv .eq. 'y')THEN
          DO 22 I=1,IB
            WRITE(10,13,ERR=19)PBU(I),ZBU(I),XBU(I),
     3      ((EBU(K,J,I),K=1,2),J=1,3),
     #ENEBU(I),(VGBU(J,I),J=1,2)
            IF(SWDF.EQ.'Y' .or. swdf .eq. 'y')THEN
            WRITE(11,11,ERR=19)PBU(I),ZBU(I),XBU(I)
            ENDIF
   22     CONTINUE
        ELSE
          WRITE(10,11,ERR=19)(PBU(I),ZBU(I),XBU(I),I=1,IB)
        ENDIF
      ENDIF
   12 FORMAT(' P=',1pE13.7,' Z=',E14.7,' X=',E13.7,E14.7,
     #' S=',E14.7,/,'EX=',2(E14.7),'EY=',2(E14.7),/,
     #'EZ=',2(E14.7),/,'ENE=',E14.7,'VGP=',E14.7,'VGZ=',E14.7)
   13 FORMAT(' P=',1pE13.7,' Z=',E14.7,' X=',E13.7,
     #/,'EX=',2(E14.7),'EY=',2(E14.7),/,
     #'EZ=',2(E14.7),/,'ENE=',E14.7,'VGP=',E14.7,'VGZ=',E14.7)
   11 FORMAT(' P=',1pE13.7,' Z=',E13.7,' X=',E13.7)
   16 FORMAT(' P=',1pE13.7,' Z=',E13.7,' X=',E13.7,E14.7,' S=',E14.7)
   18 IB=0
      RETURN
C
   19 WRITE(*,*)'WRFI: ERROR IN WRITE, FILE=',OFIL
      RETURN
C
   20 IF(IFS.EQ.0) RETURN
C
   30 WRITE(*,'(A,$)')'OUTPUT FILE:'
      READ(*,'(A)')OFIL
      IB=0
      IF(OFIL.EQ.'        ') RETURN
      IF(IFS.EQ.1) CLOSE(6)
      OPEN(10,FILE=OFIL,ERR=40)
      IFS=1
 
             rn = 1836.1 * ass(1)
             if( rn .lt. 1.) rn = 1.
             gf = 1000.* xc/rn
             vth = sqrt(3.517388e14 * ta(1)/rn)
             ptok = gf/vth

      WRITE(*,'(A,$)')'GROWTH RATES?'
      READ(*,'(A)')GR
      WRITE(*,'(A,$)')'VARIABLES FOR WDF RECONSTRUCTION?'
      READ(*,'(A)')WDFV
      SWDF='N'
      IF(WDFV.EQ.'Y' .or. wdfv .eq. 'y')THEN
         WRITE(*,'(A,$)')'WISH TO USE SEPARATE FILES (Y/N)'
         READ(*,'(A)')SWDF
         IF(SWDF.EQ.'Y' .or. swdf .eq. 'y')THEN
           WRITE(*,'(A,$)')'FILENAME?'
           READ(*,'(A)')OFIL
           OPEN(11,FILE=OFIL,ERR=40)
         ENDIF
      ENDIF
  300 FORMAT(2X,A33,/)
      DO 32 J=1,6
      ISP=SQRT(ASS(J))
      IF(ISP.LT.4) ISP=ISP+1
      WRITE(10,33)SPE(ISP),DN(J),TA(J),DD(J),AA(J,1),AA(J,2),VD(J)
      IF(SWDF.EQ.'Y')THEN
        WRITE(11,33)SPE(ISP),DN(J),TA(J),DD(J),AA(J,1),AA(J,2),VD(J)
      ENDIF
   32 CONTINUE
   33 FORMAT(A3,' N=',E9.3,' T=',F8.4,' D=',F4.2,' A=',F4.2,' B=',
     F   F4.2,' VD=',F5.2)
      RETURN
C
   40 WRITE(*,*)'WRFI: ERROR IN OPEN, FILE=',OFIL
      GOTO 30
      RETURN
      END
