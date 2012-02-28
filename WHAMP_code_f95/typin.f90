!**************************************************************
!
! File typin.f 

SUBROUTINE TYPIN(NPL,KFS)
  !        ARGUMENTS: NPL    NEW PLASMA. WHEN A PLASMA PARAMETER IS
  !                                      CHANGED, NPL IS SET TO 1.
  !                   KFS    =1  P SPECIFIED LAST.
  !                          =2  Z SPECIFIED LAST.i
  use comin
  implicit none
  ! find the kind of a high precision variable, by finding 
  ! the kind of 1.0d0
  integer, parameter :: d2p=kind(1.0d0)

  integer :: IE,IOF,IOS,KFS,NC,NPL,NV
  integer :: flag_success_reading,flag_ambiguous_character,flag_too_long_number,flag_next_number
  real :: DEC,DEK
  INTEGER,parameter :: IFILE=5
  CHARACTER CC*1
  !......... IF THE SYSTEM SUPPORTS PRINTER CONTROL CHARACTERS,
  !......... "CC" CAN BE USED TO SUPPRESS CARRIAGE RETURN/LINE FEED.
  PARAMETER( CC = '$' ) 
  CHARACTER UPP_LTRS*26, DIGITS*10, LOW_LTRS*26
  PARAMETER( UPP_LTRS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' )
  PARAMETER( LOW_LTRS = 'abcdefghijklmnopqrstuvwxyz' )
  PARAMETER( DIGITS = '0123456789' )
  CHARACTER INP*80,IC*1,input_variable*1
!  real :: ARRAY(140)
  real :: TV(2),input_number
  integer :: IOUT=2,KP=1,KZ=1
  ! in correspondence assume IOF=1
  ! **** CORRESPONDENCE BETWEEN ARRAY AND NAMES IN WHAMP:  ****
  !  ARRAY(IOF+ 0 20 30 40 50 60 70  80 90 100 110 130 131 134 137  138  139
  !            XX PP ZZ A  B  D  ASS VD DN TA  XP  CV  PM  ZM  XOI  XC   PZL
  !
  ! **** CORRESPONDENCE BETWEEN ARRAY AND NAMES IN WHAMP:  ****
  !  ARRAY(IOF+ 0 12 18 24 30 36 42  48 54 60 66 78 79 82 85  86  87
  !            XX PP ZZ A  B  D  ASS VD DN TA XP CV PM ZM XOI XC  PZL
  !
  ! 
  ! IC  - the character value
  ! INP - input line
  ! IOF - the field number in the ARRAY, usually  in code ARRAY(IV+IOF)=TV(1)*10.**TV(2)
  ! IV  - input_value, the number in array, see above
  ! NC  - number of character in line
  ! NV  - =1 number input, =0 number input ended
  ! TV  - number is constructed as TV(1)*10.**TV(2)

  NPL=0 ! default is that plasma parameters are not changed during input
  flag_success_reading=0 ! default value that did not succeed to read the line
  input_loop: do 
     input_variable=' ' ! default emppty
     call wrfi(3)
     IOS = 1

     do
        IF( IOS .EQ. 0 ) exit             
        !         write(*, '(A,$)') '#INPUT: ' 
        write(*, '(A )') '#INPUT: '  
        write(*,*) ' '
        READ( *, '( A )', IOSTAT = IOS )  INP
        IF( IOS .NE. 0 )THEN
           REWIND IFILE
        END IF
     end do

     NC=0
     input_line_loop: do
        flag_ambiguous_character=0 ! default assume non-alpha character
        flag_too_long_number=0     ! default assume input numbers are not too long
        flag_next_number=0         ! default that we enter the first number of variable
        NC=NC+1
        IF(NC.LE.80) then          ! still within input line 80 character limit
           IC=INP(NC:NC)
           IF(IC.EQ.' ') cycle input_line_loop
           IF(input_variable==' ') then 
              flag_next_number=1 
           else
              IF(INDEX(UPP_LTRS,IC).GT.0.OR.&
                   &   INDEX(LOW_LTRS,IC).GT.0) then          ! IC is alphabetic character
              else IF( INDEX( DIGITS, IC ) .GT. 0 ) then      ! IC is digit
                 !
                 TV(IE)=TV(IE)*DEK+( INDEX( DIGITS, IC ) - 1 )*DEC
                 DEC=DEC*DEK/10.
                 NV=1
                 cycle input_line_loop
              else IF(IC.NE.',') then ! strange character
                 IF(IC.EQ.'-') then
                    IF(TV(IE).NE.0.) then
                       flag_ambiguous_character=1
                       exit input_line_loop
                    end if
                    DEC=-DEC
                    cycle input_line_loop
                 end if
                 !
                 IF(IC.EQ.'.') then
                    IF(ABS(DEC).NE.1.) then 
                       flag_ambiguous_character=1
                       exit input_line_loop
                    end if
                    DEK=DEK/10.
                    DEC=DEC/10.
                    cycle input_line_loop
                 end if
                 !
                 IF(IC.NE.')') cycle input_line_loop
                 IF(DEC.NE.1.) then 
                    flag_ambiguous_character=1
                    exit input_line_loop
                 end if
                 IF(IE.NE.1.OR.TV(1).LE.0.) then 
                    flag_ambiguous_character=1
                    exit input_line_loop
                 end if
                 IOF=TV(1)+.1
                 TV(1)=0.
                 NV=0
                 cycle input_line_loop
              end if
           end IF
        end if

        if (flag_next_number==0) then ! still working on the same number input
           IF(NV.EQ.1) then
              IF(IOF.GT.10)  then ! more than 10 plasma components 
                 flag_too_long_number=1              
                 exit input_line_loop
              end if
              IF (IOF.GT.3.AND.(index('cflpz',input_variable)>0))  then  ! too many components for variable
                 flag_too_long_number=1              
                 exit input_line_loop
              end if
              IF (IOF.GT.1.AND.(index('cfl',input_variable)>0)) then ! scalars have more than 1 component
                 flag_too_long_number=1              
                 exit input_line_loop
              end if
              IF((IC.EQ.'E') .OR. (IC .eq. 'e')) then
                 IE=2
                 DEK=10.
                 DEC=1.
                 cycle input_line_loop
              end if
              input_number=TV(1)*10.**TV(2)
             
              select case(input_variable)
              case('A','a')
                 AA(IOF,1)=input_number
              case('B','b')
                 AA(IOF,2)=input_number
              case('C','c')
                 XC=input_number
              case('D','d')
                 DD(IOF)=input_number
              case('F','f')
                 XOI=input_number
              case('L','l')
                 PZL=input_number
              case('M','m')
                 ASS(IOF)=input_number
              case('N','n')
                 DN(IOF)=input_number
              case('P','p')
                 PM(IOF)=input_number
              case('T','t')
                 TA(IOF)=input_number
              case('V','v')
                 VD(IOF)=input_number
              case('Z','z')
                 ZM(IOF)=input_number
              end select
              IOF=IOF+1
              IF(input_variable=='p') KP=IOF
              IF(input_variable=='z') KZ=IOF
           end if
        end if
        !
        IF(NC.GT.80) then
           flag_success_reading=1
           exit input_line_loop
        end if
        DEK=10.
        DEC=1.
        TV(1)=0.
        TV(2)=0.
        NV=0
        IE=1
        IF(IC.EQ.',') cycle input_line_loop
        input_variable=' '
        IOF=1
        IF(IC.EQ.'A'.OR.IC.EQ.'a') input_variable='a'
        IF(IC.EQ.'B'.OR.IC.EQ.'b') input_variable='b'
        IF(IC.EQ.'C'.OR.IC.EQ.'c') input_variable='c'
        IF(IC.EQ.'D'.OR.IC.EQ.'d') input_variable='d'
        IF(IC.EQ.'F'.OR.IC.EQ.'f') input_variable='f'
        IF(IC.EQ.'H'.OR.IC.EQ.'h') then
           call print_help
           cycle input_loop
        end if
        IF(IC.EQ.'L'.OR.IC.EQ.'l') input_variable='l'
        IF(IC.EQ.'M'.OR.IC.EQ.'m') input_variable='m'
        IF(IC.EQ.'N'.OR.IC.EQ.'n') input_variable='n'
        IF(IC.EQ.'O'.OR.IC.EQ.'o') then
           IOUT=0
           cycle input_line_loop
        end if
        IF(IC.EQ.'P'.OR.IC.EQ.'p') input_variable='p'
        IF(IC.EQ.'S'.OR.IC.EQ.'s') STOP
        IF(IC.EQ.'T'.OR.IC.EQ.'t') input_variable='t'
        IF(IC.EQ.'V'.OR.IC.EQ.'v') input_variable='v'
        IF(IC.EQ.'Z'.OR.IC.EQ.'z') input_variable='z'
        !
        IF(input_variable==' ') exit input_line_loop 
        if (index('abcdmntv',input_variable)>0) NPL=1 ! plasma parameters changed
        IF(input_variable=='p') KP=1
        IF(input_variable=='p') KFS=1
        IF(input_variable=='z') KZ=1
        IF(input_variable=='z') KFS=2
     end do input_line_loop
     if (flag_ambiguous_character==1) then
        PRINT 21,IC
21      FORMAT(' AMBIGUITY CAUSED BY THE CHARACTER "',A1,'"')
        PRINT 24
24      FORMAT(' THE REST OF THE LINE IS IGNORED. PLEASE TRY AGAIN!')
        cycle input_loop
     end if
     !
     if (flag_too_long_number==1) then  
        TV(1)=TV(1)*10.**TV(2)
        PRINT 26,TV(1)
26      FORMAT(' THE VALUE',E11.3,' WILL NOT FIT IN THE VARIABLE FIELD')
        PRINT 24
        cycle input_loop
     end if
     if (flag_success_reading==0) then
        IOS = 1
        do 
           IF( IOS .EQ. 0 ) exit 
           WRITE(*,'(A,$)')'HELP, YES OR NO?'
           READ( *, '( A )', IOSTAT = IOS )  IC
           IF( IOS .NE. 0 )THEN
              REWIND IFILE
           END IF
        end do
     end if
     IF(IC.EQ.'N' .or. IC .eq. 'n') cycle input_loop
     if (flag_success_reading==0) then
        call print_help
        cycle input_loop
     end if
     IF(ABS(XOI)<=0.) then
        WRITE(*,'(A,$)')' START FREQUENCY'
        READ*,XOI
     end if
     if (KP==1) then
        WRITE(*,'(A,$)') ' PERP. WAVE VECTOR UNDEFINED!'
        cycle input_loop
     end if

     if (KP==2) then 
        PM(2)=PM(1)
        KP=3
     end if
     if (KP==3) then
        PM(3)=PM(2)-PM(1)
        KP=4
     end if
     if (KP==4) then
        if (PM(3)==0) PM(3)=10
     end if

     if (KZ==1) then
        WRITE(*,'(A,$)')' PARALLEL WAVE VECTOR UNDEFINED!'
        cycle input_loop
     end if

     if (KZ==2) then
        ZM(2)=ZM(1)
        KZ=3
     end if
     if (KZ==3) then
        ZM(3)=ZM(2)-ZM(1)
        KZ=4
     end if
     if (KZ==4) then
        if (ZM(3)==0) ZM(3)=10
     end if
     exit input_loop
  end do input_loop

  !     IF(IOUT*NEW.EQ.1)CALL WRFI(1)
  IF(IOUT.NE.1) CALL INOUT
  IOUT=1
  RETURN

contains
  subroutine print_help
    PRINT 131
131 FORMAT(' AN INPUT LINE MAY CONSIST OF UP TO 80 CHARACTERS.'/&
         &  ' THE FORMAT IS:'/' NAME1=V11,V12,V13,...NAME2=V21,V22,...NAME'/ &
         &  ' THE NAMES ARE CHOSEN FROM THE LIST:'//&
         &  ' NAME              PARAMETER'/&
         &  ' A(I)              THE ALPHA1 PARAMETER IN THE DISTRIBUTION.'/&
         &  '                   (I) IS THE COMPONENT NUMBER, I=1 - 6.'/&
         &  ' B(I)              THE ALPHA2 PARAMETER IN THE DISTRIBUTION.'/&
         &  ' C                 THE ELECTRON CYCLOTRON FREQ. IN KHZ.'/&
         &  ' D(I)              THE DELTA PARAMETER IN THE DISTRIBUTION'/&
         &  ' F                 FREQUENCY, START VALUE FOR ITERATION.'/&
         &  ' L            L=1  THE P AND Z PARAMETERS ARE INTERPRETED'/&
         &  '                   AS LOGARITHMS OF THE WAVE NUMBERS. THIS'/&
         &  '                   OPTION ALLOWS FOR LOGARITHMIC STEPS.'/&
         &  '              L=0  DEFAULT VALUE. LINEAR STEPS.'/&
         &  ' M(I)              MASS IN UNITS OF PROTON MASS.'/&
         &  ' N(I)              NUMBER DENSITY IN PART./CUBIC METER'/&
         &  ' P(I)              PERPENDICULAR WAVE VECTOR COMPONENTS.'/&
         &  '                   P(1) IS THE SMALLEST VALUE, P(2) THE'/&
         &  '                   LARGEST VALUE, AND P(3) THE INCREMENT.'/&
         &  ' S                 STOP! TERMINATES THE PROGRAM.')
    PRINT 132
132 FORMAT(  ' T(I)              TEMPERATURE IN KEV'/&
         &  ' V(I)              DRIFT VELOCITY / THERMAL VELOCITY.'/&
         &  ' Z(I)              Z-COMPONENT OF WAVE VECTOR. I HAS THE'/&
         &  '                   SAME MEANING AS FOR P(I).'/&
         &  ' A NAME WITHOUT INDEX REFERS TO THE FIRST ELEMENT, "A" IS '/&
         &  ' THUS EQUIVALENT TO "A(1)". THE VALUES V11,V12,.. MAY BE '/&
         &  ' SPECIFIED IN I-, F-, OR E-FORMAT, SEPARATED BY COMMA(,).'/&
         &  ' THE "=" IS OPTIONAL, BUT MAKES THE INPUT MORE READABLE.'/&
         &  ' EXAMPLE: INPUT:A1.,2. B(3).5,P=.1,.2,1.E-2'/&
         &  ' THIS SETS A(1)=1., A(2)=2., B(3)=.5, P(1)=.1, P(2)=.2,'/&
         &  ' AND P(3)=.01. IF THE INCREMENT P(3)/Z(3) IS NEGATIVE, P/Z'/&
         &  ' WILL FIRST BE SET TO P(2)/Z(2) AND THEN STEPPED DOWN TO'/&
         &  ' P(1)/Z(1)'/&
         &  ' THE LAST SPECIFIED OF P AND Z WILL VARY FIRST.'/&
         &  ' IF THE LETTER "O" (WITHOUT VALUE) IS INCLUDED, YOU WILL'/&
         &  '  BE ASKED TO SPECIFY A NEW OUTPUT FORMAT.'/)

  end subroutine print_help

END SUBROUTINE TYPIN

