!========================================================            
      module ScriptInput
      integer, parameter :: maxLineLen = 500   
      contains
!========================================================            
      subroutine Script_ReadParameters
      use Constants
      use ParallelVar
      use Units
      implicit none
!      integer(kind=8), intent(OUT) :: ncycle,nmoves
      integer :: i, ii, j, nArgs
      integer :: iLine, lineStat, AllocateStat
      integer :: nLines, nForceLines, lineBuffer
      integer, allocatable :: lineNumber(:), ffLineNumber(:)
      character(len=maxLineLen), allocatable :: lineStore(:)
      character(len=maxLineLen), allocatable :: forcefieldStore(:)
      character(len=25) :: command, command2, dummy
      character(len=50) :: fileName
      character(len=50) :: forcefieldFile
      

    
!      Get the filename from the command line. 
      nArgs = command_argument_count()
      if(nArgs > 1) then
        stop "This program only takes one argument"
      elseif(nArgs == 1) then
        call get_command_argument(1, fileName)
        call LoadFile(lineStore, nLines, lineNumber, fileName)
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
        call LowerCaseLine(command)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 1) then
          cycle
        endif 

        select case(trim(adjustl( command )))
        case("set")
           call setCommand( lineStore(iLine), lineStat )
           if(lineStat .eq. -1) then
             write(*,"(A,2x,I10)") "ERROR! Unknown Variable Name on Line", lineNumber(iLine)
             write(*,*) lineStore(iLine)
             stop 
           endif
        case("create")
           call createCommand( lineStore(iLine), lineStat )
           if(lineStat .eq. -1) then
             write(*,"(A,2x,I10)") "ERROR! Unknown Variable Name on Line", lineNumber(iLine)
             write(*,*) lineStore(iLine)
             stop 
           endif
!        case("forcefieldfile")
!          read(lineStore(iLine),*) dummy, command2
!          forcefieldFile =  trim( adjustl( command2 ) )
!          call LoadFile(forcefieldStore, nForceLines, ffLineNumber, forcefieldFile)
!          call ScriptForcefield(forcefieldStore)
!          if(allocated(forcefieldStore)) then
!            deallocate(forcefieldStore)
!          endif
        case default
          write(*,"(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
          write(*,*) lineStore(iLine)
          stop 
        end select
        
      enddo
!      write(*,*) "Finished Reading Input Script."      
!      deallocate(lineStore)

 
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
      use ParallelVar
      implicit none
      character(len=maxLineLen), intent(in) :: line      
      integer, intent(out) :: lineStat

      character(len=30) :: dummy, command!, stringValue
!      character(len=15) :: fileName      
      logical :: logicValue
      integer :: j
      integer :: intValue, AllocateStat
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
        case("rng_seed")
          read(line,*) dummy, command, intValue
          if(intValue < 0) then
            call system_clock(intValue)
            seed = mod(intValue,10000)
          else
            seed = intValue
          endif
          seed = p_size*seed + myid
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
           else
             write(*,*) "ERROR! The create box command has already been used and can not be called twice"
             stop
           endif
!           do i = 1, intValue
!             BoxArray(i)%
!           enddo
        case default
          lineStat = -1
      end select

     
      end subroutine   
!========================================================            
      subroutine LoadFile(lineArray, nLines, lineNumber, fileName)
      use ParallelVar, only: myid
!      use SimParameters, only: echoInput
      implicit none
      character(len=maxLineLen),allocatable,intent(inout) :: lineArray(:)
      character(len=50), intent(in) :: fileName
      integer, allocatable, intent(inout) :: lineNumber(:)

      character(len=maxLineLen),allocatable :: rawLines(:)
      character(len=25) :: command
      integer, intent(out) :: nLines
      integer :: i,iLine, nRawLines, lineStat, AllocateStat, InOutStat

      open(unit=54, file=trim(adjustl(fileName)), status='OLD', iostat=InOutStat)    

      if(InOutStat .gt. 0) then
        if(myid .eq. 0) then
          write(*,*) "The file specified in either the input script or on the command line could not be opened."
          write(*,*) "Please double check to ensure the file exists and the name in the input is accurate."
          write(*,*) "Attempted to open file:", trim(adjustl(fileName))
        endif
        stop
      endif

!      This block counts the number of lines in the input file to determine how large the lineStorage array needs to be.
      nRawLines = 0
      do iLine = 1, nint(1d7)
        read(54,*,iostat=lineStat)
        if(lineStat .lt. 0) then
          exit
        endif
        nRawLines = nRawLines + 1
      enddo
      rewind(54)


!      Read in the file line by line
      allocate(rawLines(1:nRawLines), stat = AllocateStat)

      do iLine = 1, nRawLines
        read(54,"(A)") rawLines(iLine)
!        if(echoInput) then
!          write(35,*) rawLines(iLine)        
!        endif
!        call LowerCaseLine(rawLines(iLine))
!        write(*,*) rawLines(iLine)
      enddo
      close(54) 

      nLines = 0
      do iLine = 1, nRawLines
        lineStat = 0
        call getCommand(rawLines(iLine), command, lineStat)
!         If line is empty or commented, move to the next line.         
        if(lineStat .eq. 0) then
          nLines = nLines + 1
        endif
      enddo

      allocate( lineArray(1:nLines) )
      allocate( lineNumber(1:nLines) )
      i = 0
      do iLine = 1, nRawLines
        lineStat = 0        
        call getCommand(rawLines(iLine), command, lineStat)
        if(lineStat .eq. 0) then
          i = i + 1
          call CleanLine(rawLines(iLine), lineArray(i))
!          lineArray(i) = rawLines(iLine)
          lineNumber(i) = iLine
!          write(*,*) lineNumber(i), lineArray(i)
        endif 
      enddo

      deallocate(rawLines)
    
      end subroutine

!========================================================            
!     This subrotuine searches a given input line for the first command. 
      subroutine GetCommand(line, command, lineStat)
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: line
      character(len=25), intent(out) :: command
      integer, intent(out) :: lineStat
      integer :: i, sizeLine, lowerLim, upperLim

      sizeLine = len( line )
      lineStat = 0
      i = 1
!      Find the first non-blank character in the string
      do while(i .le. sizeLine)
        if(ichar(line(i:i)) .ne. ichar(' ')) then
           !If a non-blank character is found, check first to see if it is the comment character.
          if(ichar(line(i:i)) .eq. ichar('#')) then
            lineStat = 1
            return
          else
            exit
          endif
        endif
        i = i + 1
      enddo
!      If no characters are found the line is empty, 
      if(i .ge. sizeLine) then
        lineStat = 1
        return
      endif
      lowerLim = i

!      
      do while(i .le. sizeLine)
        if(line(i:i) .eq. " ") then
          exit
        endif
        i = i + 1
      enddo
      upperLim = i

      command = line(lowerLim:upperLim)
     
      end subroutine
!========================================================            
!     This subrotuine searches a given input line for comments 
      subroutine CleanLine(inputline, cleanedLine)
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: inputline
      character(len=25), intent(out) :: cleanedLine
      integer :: i, sizeLine

      sizeLine = len(inputline)
      i = 1
!      Find the first non-blank character in the string
      do while(i .le. sizeLine)
        if(ichar(inputline(i:i)) .eq. ichar('#')) then
          exit
        endif
      enddo

      cleanedLine = inputline(1:i)
     
      end subroutine
!========================================================
      subroutine FindCommandBlock(iLine, lineStore, endCommand, lineBuffer)
      implicit none
      integer, intent(in) :: iLine
      character(len=maxLineLen), intent(in) :: lineStore(:)   
      character(len=*), intent(in) :: endCommand   
      integer, intent(out) :: lineBuffer
      logical :: found
      integer :: i, lineStat, nLines
      character(len=35) :: command 


      command = " "
      nLines = size(lineStore)
      found = .false.
      do i = iLine + 1, nLines
        call GetCommand(lineStore(i), command, lineStat)
        call LowerCaseLine(command)
!        write(*,*)  dummy
        if( trim(adjustl(command)) .eq. trim(adjustl(endCommand)) ) then
          lineBuffer = i - iLine
          found = .true.
          exit
        endif
      enddo

      if(.not. found) then
        write(*,*) "ERROR! A command block was opened in the input script, but no closing END statement found!"
        write(*,*) lineStore(iLine)
        stop
      endif

      end subroutine
!========================================================            
!     The purpose of this subroutine is to lower case a given character string. 
      subroutine LowerCaseLine(line)
      implicit none
      character(len=*),intent(inout) :: line
      integer, parameter :: offset = ichar("a") - ichar("A")
      integer :: i,sizeLine
      integer :: curVal, newVal

      sizeLine = len(line)

      do i = 1, sizeLine
        curVal = ichar(line(i:i))
        if(curVal .le. ichar("Z")) then
          if(curVal .ge. ichar("A")) then
            newVal = curVal + offSet
            line(i:i) = char(newVal)
          endif
        endif
      enddo

      end subroutine
!========================================================            
      end module
