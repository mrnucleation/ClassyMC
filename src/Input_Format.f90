!========================================================            
      module Input_Format
      integer, parameter :: maxLineLen = 100   
      contains
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
        if(lineStat < 0) then
          exit
        endif
        nRawLines = nRawLines + 1
      enddo
      rewind(54)


!      Read in the file line by line
      allocate(rawLines(1:nRawLines), stat = AllocateStat)
      IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

!      write(*,*) nRawLines

      do iLine = 1, nRawLines
        read(54,"(A)") rawLines(iLine)
!        if(echoInput) then
!          write(35,*) rawLines(iLine)        
!        endif
!        call LowerCaseLine(rawLines(iLine))
!        write(*,*) trim(rawLines(iLine))
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

      IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

      allocate( lineArray(1:nLines), stat = AllocateStat )
      allocate( lineNumber(1:nLines), stat = AllocateStat )
      i = 0
      do iLine = 1, nRawLines
        lineStat = 0        
        call getCommand(rawLines(iLine), command, lineStat)
        if(lineStat .eq. 0) then
          i = i + 1
          call CleanLine(rawLines(iLine), lineArray(i))
!          lineArray(i) = rawLines(iLine)
          lineNumber(i) = iLine
        endif 
      enddo

      IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

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
!     This subrotuine searches a given input line for the first command. 
      subroutine GetXCommand(line, command, comNum, lineStat)
      use VarPrecision
      implicit none
      integer, intent(in) :: comNum
      character(len=*), intent(in) :: line
      character(len=25), intent(out) :: command
      
      integer, intent(out) :: lineStat
      integer :: i, sizeLine, lowerLim, upperLim
      integer :: curNum


      sizeLine = len( line )
      lineStat = 0
      i = 1
      curNum = 0
!      Find the first non-blank character in the string
      do while(curNum < comNum)
        do while(i .le. sizeLine)
          if(ichar(line(i:i)) .ne. ichar(' ')) then
            exit
          endif
          i = i + 1
        enddo
!        If no characters are found the line is empty, 
        if(i .ge. sizeLine) then
          lineStat = 1
          return
        endif
        lowerLim = i
      
        do while(i .le. sizeLine)
          if(line(i:i) .eq. " ") then
            exit
          endif
          i = i + 1
        enddo
        if(i .ge. sizeLine) then
          lineStat = 1
          return
        endif
        upperLim = i
        curNum = curNum + 1
      enddo

      command = line(lowerLim:upperLim)
     
      end subroutine
!========================================================            
!     This subrotuine searches a given input line for comments 
      subroutine CleanLine(inputline, cleanedLine)
      use VarPrecision
      implicit none
      character(len=*), intent(in) :: inputline
      character(len=maxLineLen), intent(out) :: cleanedLine
      integer :: i, sizeLine

      sizeLine = len(inputline)
      i = 1
!      Find the first non-blank character in the string
      do while(i .le. sizeLine)
        if(ichar(inputline(i:i)) .eq. ichar('#')) then
          exit
        endif
        i = i + 1
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
