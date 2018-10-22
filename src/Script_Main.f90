!========================================================            
! 
!========================================================            
      module ScriptInput
      use Input_Format
!      integer, parameter :: maxLineLen = 100   
      contains
!========================================================            
      subroutine Script_ReadParameters
      use Constants
      use Input_Forcefield
      use Input_AnalysisType, only: Script_AnalysisType
      use Input_Sampling, only: Script_SamplingType
      use Input_NeighType, only: Script_NeighType
      use Input_Initialize, only: Script_Initialize
      use SimControl, only: TimeStart, TimeEnd
      use SimMonteCarlo, only: RunMonteCarlo
      use ParallelVar
      use Units
      implicit none
      integer :: i, ii, j, nArgs
      integer :: iLine, lineStat, AllocateStat
      integer :: nLines, nForceLines, lineBuffer
      integer, allocatable :: lineNumber(:)

      character(len=maxLineLen), allocatable :: lineStore(:)
      character(len=30) :: command, command2, dummy
      character(len=50) :: fileName
      character(len=50) :: forcefieldFile

    
!      Get the filename from the command line. 
      nArgs = command_argument_count()
      if(nArgs > 0) then
!        stop "This program only takes one argument"
           
!      elseif(nArgs == 1) then
        call get_command_argument(1, fileName)
        nLines = 0
        call LoadFile(lineStore, nLines, lineNumber, fileName)
        if(nLines == 0) then
          write(*,*) "ERROR! Input file is empty or could not be read!"
          stop
        else
          write(nout, *) "File successfully loaded!"
        endif
      elseif(nArgs == 0) then
        write(*,*) "ERROR! No Input File has been given!"
        stop
      endif

!      This block counts the number of lines in the input file to determine how large the lineStorage array needs to be.

      do iLine = 1, nLines
        call LowerCaseLine(lineStore(iLine))
      enddo
!      call setDefaults(seed, screenEcho)

      lineBuffer = 0
      do iLine = 1, nLines
!        write(*,*) trim(adjustl(lineStore(iLine)))
        if(lineBuffer > 0) then
          lineBuffer = lineBuffer - 1
          cycle
        endif
        lineStat = 0        
        call getCommand(lineStore(iLine), command, lineStat)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 1) then
          cycle
        endif 

        select case(trim(adjustl( command )))
          case("create")
            call createCommand(iLine, linestore, lineBuffer, lineStat)

          case("forcefield")
            call GetXCommand(lineStore(iLine), filename, 2, lineStat)  
            call Script_ReadFieldFile(filename, lineStat)

          case("modify")
            call modifyCommand( lineStore(iLine), lineStat )

          case("neighlist")
            call Script_NeighType( lineStore(iLine), lineStat)

          case("samplingtype")
            call Script_SamplingType(iLine, lineStore, lineStat)

          case("set")
            call setCommand( lineStore(iLine), lineStat )

          case("run")
            call CPU_TIME(TimeStart)
            call RunMonteCarlo
            call CPU_TIME(TimeEnd)


            
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
      if(allocated(lineStore)) then
        deallocate(lineStore)
      endif

      call Script_Initialize
!      call Script_SafetyCheck

      end subroutine
!========================================================            
      subroutine setCommand(line, lineStat)
      use Common_MolInfo
      use Common_NeighData
      use ParallelVar
      use RandomGen, only: initSeed
      use SimControl, only: nMoves, nCycles, screenFreq
      use Units, only: outEngUnit, outLenUnit, outAngUnit,  &
                       FindEngUnit, FindLengthUnit, FindAngularUnit
      use VarPrecision
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: command, command2
      logical :: logicValue
      integer :: intValue
      real(dp) :: realValue
      

      lineStat  = 0

      call GetXCommand(line, command, 2, lineStat)
      select case(trim(adjustl(command)))
        case("cycles")
          call GetXCommand(line, command, 3, lineStat)
          read(command, *) realValue  
          nCycles = nint(realValue, kind=8)

        case("moves")
          call GetXCommand(line, command, 3, lineStat)
          read(command, *) realValue  
          nMoves = nint(realValue, kind=8)

!        case("screenecho")
!          read(line,*) dummy, command, logicValue
!          screenEcho = logicValue

        case("rng_seed")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2,*) intValue
          initSeed = intValue

        case("angleunits")
          call GetXCommand(line, command2, 3, lineStat)
          outAngUnit = FindAngularUnit(command2)

        case("distunits")
          call GetXCommand(line, command2, 3, lineStat)
          outLenUnit = FindLengthUnit(command2)


        case("energyunits")
          call GetXCommand(line, command2, 3, lineStat)
          outEngUnit = FindEngUnit(command2)

        case("neighskin")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          neighSkin = realValue

        case("screenfrequency")
          call GetXCommand(line, command2, 3, lineStat)
          read(command2, *) realValue
          screenFreq = nint(realValue)

        case default
          lineStat = -1
      end select

     
      end subroutine
!========================================================            
      subroutine createCommand(iLine, linestore, lineBuffer, lineStat)
      use AnalysisData, only: AnalysisArray, analyCommon
      use BoxData, only: BoxArray
      use BoxPresets, only: Preset
      use TrajData, only: TrajArray
      use MCMoveData, only: Moves, MoveProb
      use ForcefieldData, only: EnergyCalculator, nForceFields

      use Input_SimBoxes, only: Script_BoxType
      use Input_Forcefield, only: Script_FieldType
      use Input_AnalysisType, only: Script_AnalysisType
      use Input_Constraint, only: Script_Constraint
      use Input_TrajType, only: Script_TrajType
      use Input_Moves, only: Script_MCMoves
      use Input_LoadCoords, only: Script_ReadCoordFile

      use VarPrecision
      use Units
      implicit none
      integer, intent(in) :: iLine
      character(len=maxLineLen), intent(in) :: linestore(:) 
      integer, intent(out) :: lineStat, lineBuffer

      character(len=30) :: dummy, command, command2
      character(len=50) :: fileName
      logical :: logicValue
      integer :: i, curLine
      integer :: intValue, AllocateStat, nItems
      real(dp) :: realValue
      

      lineStat  = 0
      AllocateStat = 0

      read(linestore(iLine),*) dummy, command
      call FindCommandBlock(iLine, lineStore, "end_create", lineBuffer)
      nItems = lineBuffer - 1
!      write(*,*) nItems

      select case(trim(adjustl(command)))
        case("analysis") 
           if( .not. allocated(AnalysisArray) ) then
             allocate(AnalysisArray(1:nItems), stat = AllocateStat)
             allocate(analyCommon(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_AnalysisType(linestore(curLine), i, lineStat)
             enddo   
             do i = 1, nItems
               call AnalysisArray(i) % func % CastCommonType(analyCommon(i)%val)
             enddo   

           else
             write(*,*) "ERROR! The create analysis command has already been used and can not be called twice"
             stop
           endif

         !-------------------------------------------------------------------------------------
        case("boxes")
           if( .not. allocated(BoxArray) ) then
             allocate(BoxArray(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call GetXCommand(lineStore(curLine), command2, 1, lineStat)
               if( trim(adjustl(command2)) == "fromfile") then
                   call GetXCommand(lineStore(curLine), command2, 2, lineStat)
                   fileName = ""
                   fileName(1:30) = command2(1:30)
                   call Script_ReadCoordFile(fileName, i, lineStat)

               elseif( trim(adjustl(command2)) == "preset") then
                   call Preset(curLine, lineStore, i, lineStat)

               else
                   call Script_BoxType(linestore(curLine), i, lineStat)
               endif
               BoxArray(i)%box%boxID = i 
             enddo             
           else
             write(*,*) "ERROR! The create box command has already been used and can not be called twice"
             stop
           endif

         !-------------------------------------------------------------------------------------
        case("constraint")
           if( .not. allocated(BoxArray) ) then
             stop "Box array not allocated!"
           endif
           call GetXCommand(lineStore(iLine), command2, 3, lineStat)
           read(command2,*) intValue
           call Script_Constraint(lineStore, iLine, intValue, lineBuffer, lineStat)


!        case("energyfunctions")
!           if( .not. allocated(EnergyCalculator) ) then
!             nForceFields = nItems
!             allocate(EnergyCalculator(1:nItems), stat = AllocateStat)
!             do i = 1, nItems
!               curLine = iLine + i
!               call Script_FieldType(linestore(curLine), i, lineStat)
!               call EnergyCalculator(i)%Method%Constructor
!             enddo   
!           else
!             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
!             stop
!           endif
!
         !-------------------------------------------------------------------------------------
        case("moves") 
           if( .not. allocated(Moves) ) then
             nForceFields = nItems
             allocate(Moves(1:nItems), stat = AllocateStat)
             allocate(MoveProb(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_MCMoves(linestore(curLine), i, lineStat)
             enddo   
           else
             write(*,*) "ERROR! The create moves command has already been used and can not be called twice"
             stop
           endif

         !-------------------------------------------------------------------------------------
        case("trajectory") 
           if( .not. allocated(TrajArray) ) then
             allocate(TrajArray(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_TrajType(linestore(curLine), i, lineStat)
             enddo   
           else
             write(*,*) "ERROR! The create trajectory command has already been used and can not be called twice"
             stop
           endif
 
        case default
          write(*,*) command
          lineStat = -1
      end select

      IF (AllocateStat /= 0) then
        write(*,*) AllocateStat
        STOP "Allocation Error in Create Command"
      endif
     
      end subroutine   
!========================================================            
      subroutine modifyCommand(line, lineStat)
      use BoxData, only: BoxArray
      use Template_SimBox
      use ForcefieldData, only: EnergyCalculator, nForceFields
      use MCMoveData, only: Moves
      use CommonSampling, only: sampling
      use VarPrecision
      use Units
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: dummy, command, command2
      logical :: logicValue
      integer :: intValue, AllocateStat
      real(dp) :: realValue
      

      lineStat  = 0
      AllocateStat = 0
!      read(line,*) dummy, command
      call GetXCommand(line, command, 2, lineStat)
      call LowerCaseLine(command)
      select case(adjustl(trim(command)))
        case("box")
           call GetXCommand(line, command2, 3, lineStat)
           read(command2, *) intValue
           call BoxArray(intValue) % box % IOProcess(line, lineStat)

        case("move")
           call GetXCommand(line, command2, 3, lineStat)
           read(command2, *) intValue
           call Moves(intValue) % Move % ProcessIO(line, lineStat)

        case("sampling")
           call Sampling % ProcessIO(line, lineStat)


        case default
           lineStat = -1
      end select

      IF (AllocateStat /= 0) STOP "Allocation Error in the Modify Command"
     
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
