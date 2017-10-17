!========================================================            
      module ScriptInput
      use Input_Format
!      integer, parameter :: maxLineLen = 100   
      contains
!========================================================            
      subroutine Script_ReadParameters
      use Constants
      use ParallelVar
      use Units
      use Input_Forcefield
      implicit none
!      integer(kind=8), intent(OUT) :: ncycle,nmoves
      integer :: i, ii, j, nArgs
      integer :: iLine, lineStat, AllocateStat
      integer :: nLines, nForceLines, lineBuffer
      integer, allocatable :: lineNumber(:), ffLineNumber(:)
      character(len=maxLineLen), allocatable :: lineStore(:)
      character(len=maxLineLen), allocatable :: forcefieldStore(:)
      character(len=30) :: command, command2, dummy
      character(len=50) :: fileName
      character(len=50) :: forcefieldFile
      

    
!      Get the filename from the command line. 
      nArgs = command_argument_count()
      if(nArgs > 1) then
        stop "This program only takes one argument"
      elseif(nArgs == 1) then
        call get_command_argument(1, fileName)
        call LoadFile(lineStore, nLines, lineNumber, fileName)
        write(*,*) "File successfully loaded!"
      elseif(nArgs == 0) then
        write(*,*) "ERROR! No Input File has been given!"
        stop
      endif
!      Read in the script file


!      fileName = "ScriptTest.dat"

!      This block counts the number of lines in the input file to determine how large the lineStorage array needs to be.

!      call setDefaults(seed, screenEcho)

      lineBuffer = 0
      do iLine = 1, nLines
        if(lineBuffer .gt. 0) then
          lineBuffer = lineBuffer - 1
          cycle
        endif
        lineStat = 0        
!        write(*,*) trim(adjustl(lineStore(iLine)))
        call getCommand(lineStore(iLine), command, lineStat)
!        call GetXCommand(lineStore(iLine), dummy, 2, lineStat)
!        write(*,*) dummy, command
        call LowerCaseLine(command)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 1) then
          cycle
        endif 

        select case(trim(adjustl( command )))
        case("set")
           call setCommand( lineStore(iLine), lineStat )
        case("create")
           call createCommand( lineStore(iLine), lineStat )
        case("modify")
           call modifyCommand( lineStore(iLine), lineStat )
        case("forcefield")
          call Script_Forcefield( lineStore(iLine), lineStat )
        case default
          write(*,"(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
          write(*,*) trim(adjustl(lineStore(iLine)))
          stop 
        end select

        ! Ensure that the called processes exited properly.
        if(lineStat .eq. -1) then
          write(*,"(A,1x,I10)") "ERROR! Parameters for command on line:", lineNumber(iLine)
          write(*, "(A)") "could not be understood. Please check command for accuracy and try again."
          write(*,*) trim(adjustl(lineStore(iLine)))
          stop 
        endif
        
      enddo
!      write(*,*) "Finished Reading Input Script."      
      deallocate(lineStore)

 
!      call ReadInitialConfiguration
!      write(*,*) "Finished Reading Initial Configuration."
!      call RecenterCoordinates
!      call ReadInitialGasPhase     
!      write(*,*) "Finished Setting up Gas Phase Configuration." 
!      if(useWHAM) then
!        if(.not. useUmbrella) then
!          write(*,*) "ERROR! The WHAM method can not be used if no Umbrella sampling variables are given!"
!          stop
!        endif
!        nWhamItter = ceiling(dble(ncycle)/dble(intervalWHAM))
!        call WHAM_Initialize
!      endif



      end subroutine
!========================================================            
      subroutine setCommand(line, lineStat)
      use Units
      use VarPrecision
      use Common_MolDef
      use ParallelVar
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: dummy, command!, stringValue
!      character(len=15) :: fileName      
      logical :: logicValue
      integer :: j
      integer :: intValue
      real(dp) :: realValue
      

      lineStat  = 0

      read(line,*) dummy, command
      call LowerCaseLine(command)
      select case(trim(adjustl(command)))
!        case("cycles")
!          read(line,*) dummy, command, realValue  
!          nCycle = nint(realValue)
!        case("moves")
!          read(line,*) dummy, command, realValue        
!          ncycle2 = nint(realValue)   
!        case("screenecho")
!          read(line,*) dummy, command, logicValue
!          screenEcho = logicValue
        case("atomtypes")
          read(line,*) dummy, command, intValue
          nAtomTypes = intValue
        case("rng_seed")
          read(line,*) dummy, command, intValue
          if(intValue < 0) then
            call system_clock(intValue)
            seed = mod(intValue,10000)
          else
            seed = intValue
          endif
          seed = p_size*seed + myid
          write(*,*) seed
!        case("out_energyunits")
!          read(line,*) dummy, command, outputEngUnits
!          outputEConv = FindEngUnit(outputEngUnits)
!        case("out_distunits")
!          read(line,*) dummy, command, outputLenUnits   
!          outputLenConv = FindLengthUnit(outputLenUnits)
        case default
          lineStat = -1
      end select

     
      end subroutine
!========================================================            
      subroutine createCommand(line, lineStat)
      use BoxData, only: BoxArray
      use CubicBoxDef, only: CubeBox
      use ForcefieldData, only: EnergyCalculator, nForceFields
      use VarPrecision
      use Units
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: dummy, command!, stringValue
      logical :: logicValue
      integer :: i
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      lineStat  = 0

      read(line,*) dummy, command
      call LowerCaseLine(command)
      select case(trim(adjustl(command)))
        case("boxes")
           read(line,*) dummy, command, intValue
           if( .not. allocated(BoxArray) ) then
             allocate(BoxArray(1:intValue), stat = AllocateStat)

             allocate(CubeBox::BoxArray(1)%box, stat = AllocateStat)
           else
             write(*,*) "ERROR! The create box command has already been used and can not be called twice"
             stop
           endif

        case("energycalculators")
           read(line,*) dummy, command, intValue
           nForceFields = intValue
           if( .not. allocated(EnergyCalculator) ) then
             allocate(EnergyCalculator(1:intValue), stat = AllocateStat)
           else
             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
             stop
           endif
 
        case default
          lineStat = -1
      end select

      IF (AllocateStat /= 0) STOP "*** Not enough memory ***"
     
      end subroutine   
!========================================================            
      subroutine modifyCommand(line, lineStat)
      use BoxData, only: BoxArray
      use ForcefieldData, only: EnergyCalculator, nForceFields
      use VarPrecision
      use Units
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: dummy, command, command2
      logical :: logicValue
      integer :: i
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      lineStat  = 0

!      read(line,*) dummy, command
      call GetXCommand(line, command, 2, lineStat)
      call LowerCaseLine(command)
      select case(adjustl(trim(command)))
        case("box")
           call GetXCommand(line, command2, 3, lineStat)
           read(command2, *) intValue

           call BoxArray(intValue)%box%IOProcess(line, lineStat)
        case default
           lineStat = -1
      end select

      IF (AllocateStat /= 0) STOP "*** Not enough memory ***"
     
      end subroutine   
!========================================================            
!     The purpose of this subroutine is to lower case a given character string. 
      subroutine IO_ErrorControl(iLine, lineNumber, lineStore, lineStat)
      implicit none
      integer, intent(in) :: iLine, lineStat
      integer, intent(in) :: lineNumber(:)
      character(len=maxLineLen), intent(in) :: lineStore(:)

      if(lineStat .eq. -1) then
        write(*,"(A,2x,I10)") "ERROR! Unknown Variable Name on Line", lineNumber(iLine)
        write(*,*) lineStore(iLine)
        stop 
      endif

      end subroutine
!========================================================            
      end module
