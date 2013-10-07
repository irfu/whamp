PROGRAM WHAMP 
  use comin
  use comcout
  use comoutput
  implicit none
  !integer, parameter :: d2p=kind(1.0d0)
  integer, parameter :: d2p=8
  logical :: isChangedPlasmaModel
  integer :: KFS,J
  CHARACTER FILENAME*(80)

  printDebugInfo = .true.

  CALL READ_INPUT_FILE(FILENAME)
  call whamp_engine ! call to get plasma parameters only
  !
  plasma_update_loop: do
     call print_plasma_parameters()
     isChangedPlasmaModel = .false.
     !                  ****  ASK FOR INPUT!  **** 
     typin_loop: do 
        !for new plasma skip calling typin until convergence checked
        if(.not.isChangedPlasmaModel) call  TYPIN(isChangedPlasmaModel,KFS) 
        if (KFS==1) cycleZFirst = 0
        if (KFS==2) cycleZFirst = 1
        if(isChangedPlasmaModel)      cycle plasma_update_loop
        isChangedPlasmaModel = .false.
        call whamp_engine
        write(*,*) 'PM=',PM
        write(*,*) 'cycleZfirst=',cycleZfirst
     end do typin_loop
  end do plasma_update_loop
  contains
  subroutine print_plasma_parameters
 !                  ****  PRINT PLASMA PARAMETERS.  ****
     PRINT 101,PX,XC,DEN
101  FORMAT('# PLASMA FREQ.:',  F11.4,&
          &       'KHZ GYRO FREQ.:',  F10.4,  'KHZ   ',&
          &       'ELECTRON DENSITY:',1PE11.5,  'M-3'   )
     DO  J=1,JMA
102     FORMAT('# ',  A3,  '  DN=',1PE12.5,  '  T=',0PF9.5,  '  D=',  F4.2,&
             &'  A=',  F4.2,  '  B=',  F4.2,  ' VD=',  F5.2)
        PRINT 102,species_symbol(ASS(J)),DN(J),TA(J),DD(J),AA(J,1),AA(J,2),VD(J)
     end do
     !
  end subroutine
  pure function species_symbol(mass) result(symbol)
          real(kind=d2p),intent(in) :: mass
          character(5)   :: symbol
          if (mass==0) then
              symbol = 'e-'
          else if (mass==1) then
              symbol = 'H+'
          else if (mass==2) then
              symbol = 'He++'
          else if (mass==4) then
              symbol = 'He+'
          else if (mass==16) then
              symbol = 'O+'
          else
             write(symbol,'(a,I3)') 'm=',mass 
          endif
  end function
end program WHAMP
