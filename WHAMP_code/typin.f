***************************************************************
*
* File typin.f 

      SUBROUTINE TYPIN(NPL,KFS)
C        ARGUMENTS: NPL    NEW PLASMA. WHEN A PLASMA PARAMETER IS
C                                      CHANGED, NPL IS SET TO 1.
C                   KFS    =1  P SPECIFIED LAST.
C                          =2  Z SPECIFIED LAST.
      IMPLICIT REAL*8 ( A - H, O - Z )
      INTEGER IFILE
      PARAMETER( IFILE = 5 )
      CHARACTER CC*1
*......... IF THE SYSTEM SUPPORTS PRINTER CONTROL CHARACTERS,
*......... "CC" CAN BE USED TO SUPPRESS CARRIAGE RETURN/LINE FEED.
      PARAMETER( CC = '$' ) 
      CHARACTER UPP_LTRS*26, DIGITS*10, LOW_LTRS*26
      PARAMETER( UPP_LTRS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' )
      PARAMETER( LOW_LTRS = 'abcdefghijklmnopqrstuvwxyz' )
      PARAMETER( DIGITS = '0123456789' )
      CHARACTER INP*80,IC*1
      COMMON /XPZ/ ARRAY(140)
      DIMENSION TV(2)
      DATA IOUT/2/,KP,KZ/1,1/
C **** CORRESPONDENCE BETWEEN ARRAY AND NAMES IN WHAMP:  ****
C  ARRAY(IOF+ 0 20 30 40 50 60 70  80 90 100 110 130 131 134 137  138  139
C            XX PP ZZ A  B  D  ASS VD DN TA  XP  CV  PM  ZM  XOI  XC   PZL
C
C **** CORRESPONDENCE BETWEEN ARRAY AND NAMES IN WHAMP:  ****
C  ARRAY(IOF+ 0 12 18 24 30 36 42  48 54 60 66 78 79 82 85  86  87
C            XX PP ZZ A  B  D  ASS VD DN TA XP CV PM ZM XOI XC  PZL
C
    1 IV=1000
      call wrfi(3)
      IOS = 1
* ... WHILE( IOS .NE. 0 )DO
   70 CONTINUE
      IF( IOS .EQ. 0 )GOTO 75             
c         write(*, '(A,$)') '#INPUT: ' 
         write(*, '(A )') '#INPUT: '  
         write(*,*) ' '
         READ( *, '( A )', IOSTAT = IOS )  INP
         IF( IOS .NE. 0 )THEN
            REWIND IFILE
         END IF
* ... END WHILE
      GOTO 70
   75 CONTINUE
      NC=0
    4 NC=NC+1
      IF(NC.GT.80) GOTO 6
      IC=INP(NC:NC)
      IF(IC.EQ.' ') GOTO 4
      IF(IV.GT.140) GOTO 7
      IF(INDEX(UPP_LTRS,IC).GT.0.OR.
     #   INDEX(LOW_LTRS,IC).GT.0)GOTO 6
      IF( INDEX( DIGITS, IC ) .EQ. 0 )GOTO 5
C
      TV(IE)=TV(IE)*DEK+( INDEX( DIGITS, IC ) - 1 )*DEC
      DEC=DEC*DEK/10.
      NV=1
      GOTO 4
C
    5 IF(IC.NE.',') GOTO 10
    6 IF(NV.NE.1) GOTO 7
      IF(IOF.GT.6) GOTO 25
      IF(IOF.GT.3.AND.IV.GT.110) GOTO 25
      IF(IOF.GT.1.AND.IV.GT.137) GOTO 25
      IF(IC.EQ.'E' .OR. IC .eq. 'e') GOTO 9
      ARRAY(IV+IOF)=TV(1)*10.**TV(2)
      IOF=IOF+1
      IF(IV.EQ.131) KP=IOF
      IF(IV.EQ.134) KZ=IOF
C
    7 IF(NC.GT.80) GOTO 50
      DEK=10.
      DEC=1.
      TV(1)=0.
      TV(2)=0.
      NV=0
      IE=1
      IF(IC.EQ.',') GOTO 4
      IV=1000
      IOF=1
      IF(IC.EQ.'A'.OR.IC.EQ.'a') IV=40
      IF(IC.EQ.'B'.OR.IC.EQ.'b') IV=50
      IF(IC.EQ.'C'.OR.IC.EQ.'c') IV=138
      IF(IC.EQ.'D'.OR.IC.EQ.'d') IV=60
      IF(IC.EQ.'F'.OR.IC.EQ.'f') IV=137
      IF(IC.EQ.'H'.OR.IC.EQ.'h') GOTO 30
      IF(IC.EQ.'L'.OR.IC.EQ.'l') IV=139
      IF(IC.EQ.'M'.OR.IC.EQ.'m') IV=70
      IF(IC.EQ.'N'.OR.IC.EQ.'n') IV=90
      IF(IC.EQ.'O'.OR.IC.EQ.'o') GOTO 8
      IF(IC.EQ.'P'.OR.IC.EQ.'p') IV=131
      IF(IC.EQ.'S'.OR.IC.EQ.'s') STOP
      IF(IC.EQ.'T'.OR.IC.EQ.'t') IV=100
      IF(IC.EQ.'V'.OR.IC.EQ.'v') IV=80
      IF(IC.EQ.'Z'.OR.IC.EQ.'z') IV=134
C
      IF(IV.GE.1000) GOTO 27
      IF(IV.LT.110.OR.IV.EQ.138) NPL=1
      IF(IV.EQ.131) KP=1
      IF(IV.EQ.131) KFS=1
      IF(IV.EQ.134) KZ=1
      IF(IV.EQ.134) KFS=2
      GOTO 4
C
    8 IOUT=0
      GOTO 4
C
    9 IE=2
      DEK=10.
      DEC=1.
      GOTO 4
C
   10 IF(IC.NE.'-') GOTO 11
      IF(TV(IE).NE.0.) GOTO 20
      DEC=-DEC
      GOTO 4
C
   11 IF(IC.NE.'.') GOTO 12
      IF(ABS(DEC).NE.1.) GOTO 20
      DEK=DEK/10.
      DEC=DEC/10.
      GOTO 4
C
   12 IF(IC.NE.')') GOTO 4
      IF(DEC.NE.1.) GOTO 20
      IF(IE.NE.1.OR.TV(1).LE.0.) GOTO 20
      IOF=TV(1)+.1
      TV(1)=0.
      NV=0
      GOTO 4
C
   20 PRINT 21,IC
   21 FORMAT(' AMBIGUITY CAUSED BY THE CHARACTER "',A1,'"')
   23 PRINT 24
   24 FORMAT(' THE REST OF THE LINE IS IGNORED. PLEASE TRY AGAIN!')
      GOTO 1
C
   25 TV(1)=TV(1)*10.**TV(2)
      PRINT 26,TV(1)
   26 FORMAT(' THE VALUE',E11.3,' WILL NOT FIT IN THE VARIABLE FIELD')
      GOTO 23
C
   27 CONTINUE
      IOS = 1
* ... WHILE( IOS .NE. 0 )DO
   80 CONTINUE
      IF( IOS .EQ. 0 )GOTO 85 
         WRITE(*,'(A,$)')'HELP, YES OR NO?'
         READ( *, '( A )', IOSTAT = IOS )  IC
         IF( IOS .NE. 0 )THEN
            REWIND IFILE
         END IF
* ... END WHILE
      GOTO 80
   85 CONTINUE
      IF(IC.EQ.'N' .or. IC .eq. 'n') GOTO 1
   30 PRINT 31
   31 FORMAT(' AN INPUT LINE MAY CONSIST OF UP TO 80 CHARACTERS.'/
     F  ' THE FORMAT IS:'/' NAME1=V11,V12,V13,...NAME2=V21,V22,...NAME'/
     F  ' THE NAMES ARE CHOSEN FROM THE LIST:'//
     F  ' NAME              PARAMETER'/
     F  ' A(I)              THE ALPHA1 PARAMETER IN THE DISTRIBUTION.'/
     F  '                   (I) IS THE COMPONENT NUMBER, I=1 - 6.'/
     F  ' B(I)              THE ALPHA2 PARAMETER IN THE DISTRIBUTION.'/
     F  ' C                 THE ELECTRON CYCLOTRON FREQ. IN KHZ.'/
     F  ' D(I)              THE DELTA PARAMETER IN THE DISTRIBUTION'/
     F  ' F                 FREQUENCY, START VALUE FOR ITERATION.'/
     F  ' L            L=1  THE P AND Z PARAMETERS ARE INTERPRETED'/
     F  '                   AS LOGARITHMS OF THE WAVE NUMBERS. THIS'/
     F  '                   OPTION ALLOWS FOR LOGARITHMIC STEPS.'/
     F  '              L=0  DEFAULT VALUE. LINEAR STEPS.'/
     F  ' M(I)              MASS IN UNITS OF PROTON MASS.'/
     F  ' N(I)              NUMBER DENSITY IN PART./CUBIC METER'/
     F  ' P(I)              PERPENDICULAR WAVE VECTOR COMPONENTS.'/
     F  '                   P(1) IS THE SMALLEST VALUE, P(2) THE'/
     F  '                   LARGEST VALUE, AND P(3) THE INCREMENT.'/
     F  ' S                 STOP! TERMINATES THE PROGRAM.')
      PRINT 32
   32 FORMAT(  ' T(I)              TEMPERATURE IN KEV'/
     F  ' V(I)              DRIFT VELOCITY / THERMAL VELOCITY.'/
     F  ' Z(I)              Z-COMPONENT OF WAVE VECTOR. I HAS THE'/
     F  '                   SAME MEANING AS FOR P(I).'/
     F  ' A NAME WITHOUT INDEX REFERS TO THE FIRST ELEMENT, "A" IS '/
     F  ' THUS EQUIVALENT TO "A(1)". THE VALUES V11,V12,.. MAY BE '/
     F  ' SPECIFIED IN I-, F-, OR E-FORMAT, SEPARATED BY COMMA(,).'/
     F  ' THE "=" IS OPTIONAL, BUT MAKES THE INPUT MORE READABLE.'/
     F  ' EXAMPLE: INPUT:A1.,2. B(3).5,P=.1,.2,1.E-2'/
     F  ' THIS SETS A(1)=1., A(2)=2., B(3)=.5, P(1)=.1, P(2)=.2,'/
     F  ' AND P(3)=.01. IF THE INCREMENT P(3)/Z(3) IS NEGATIVE, P/Z'/
     F  ' WILL FIRST BE SET TO P(2)/Z(2) AND THEN STEPPED DOWN TO'/
     F  ' P(1)/Z(1)'/
     F  ' THE LAST SPECIFIED OF P AND Z WILL VARY FIRST.'/
     F  ' IF THE LETTER "O" (WITHOUT VALUE) IS INCLUDED, YOU WILL'/
     F  '  BE ASKED TO SPECIFY A NEW OUTPUT FORMAT.'/)
      GOTO 1
c   50 IF(ARRAY(138).GT.0.) GOTO 51
   50 IF(ABS(ARRAY(138)).GT.0.) GOTO 51
      WRITE(*,'(A,$)')' START FREQUENCY'
      READ*,ARRAY(138)
   51 GOTO(52,53,54,55) KP
   52 WRITE(*,'(A,$)') ' PERP. WAVE VECTOR UNDEFINED!'
      GOTO 1
   53 ARRAY(133)=ARRAY(132)
   54 ARRAY(134)=ARRAY(133)-ARRAY(132)
   55 IF(ARRAY(134).EQ.0.) ARRAY(134)=10.
      GOTO (56,57,58,59) KZ
   56 WRITE(*,'(A,$)')' PARALLEL WAVE VECTOR UNDEFINED!'
      GOTO 1
   57 ARRAY(136)=ARRAY(135)
   58 ARRAY(137)=ARRAY(136)-ARRAY(135)
   59 IF(ARRAY(137).EQ.0.) ARRAY(137)=10. 
C     IF(IOUT*NEW.EQ.1)CALL WRFI(1)
      IF(IOUT.NE.1) CALL INOUT
      IOUT=1
      RETURN
      END

