!================================================================================
module Input_LoadCoords
  use VarPrecision
  integer, parameter :: coordChunkLen = 300
  integer, parameter :: maxCoordRecordLen = 4096
  contains
!================================================================================
  subroutine Script_ReadCoordFile(filename, boxNum, lineStat)
    use BoxData, only: BoxArray
    use Input_Format, only: maxLineLen, GetXCommand, ReplaceText
    use Input_SimBoxes, only: Script_BoxType
    use APeriodicOrthoDef, only: APeriodicOrthoBox
    use ParallelVar, only: nout, myid
    implicit none
    character(len=50), intent(inout) :: fileName      
    integer, intent(in) :: boxNum
    integer, intent(out) :: lineStat

    integer :: iCharacter
    integer :: i, nLines,AllocateStat, iLine, lineBuffer
    integer, allocatable :: lineNumber(:)
    character(len=coordChunkLen), allocatable :: lineStore(:)

    character(len=30) :: command, val
    character(len=maxLineLen) :: fullLine, newline
    character(len=30) :: idString
    integer :: nAtoms



    lineStat  = 0
    AllocateStat = 0

    !Replace the & character with the thread ID in the file name
    do iCharacter = 1, len(filename)
      if(filename(iCharacter:iCharacter) == "&") then
        write(idString, *) myid
        filename = ReplaceText(filename, "&", trim(adjustl(idString)))
        exit
      endif
    enddo


    write(nout,*) "Loading file ", trim(adjustl(fileName))
    call LoadCoordFile(lineStore, nLines, lineNumber, fileName)
    lineBuffer = 0
    nAtoms = 0
    do iLine = 1, nLines
      if(lineBuffer > 0) then
        lineBuffer = lineBuffer - 1
        cycle
      endif
      lineStat = 0
      fullLine = lineStore(iLine)
      call GetXCommand(fullLine, command, 1, lineStat)
      if(lineStat .eq. 1) then
        cycle
      endif

!      write(nout,*) trim(adjustl(lineStore(iLine)))
      select case(trim(adjustl( command )))
        case("boxtype")
          call GetXCommand(fullLine, val, 2, lineStat)
          newline = ""
          newline(1:30) = val(1:30)
          call Script_BoxType(newline, boxNum, lineStat)
        case("dimension")
          call BoxArray(boxNum)%box%LoadDimension(fullLine, lineStat)
        case("boundary")
          select type(box => BoxArray(boxNum)%box)
            type is(APeriodicOrthoBox)
              call box%LoadBoundaryFlags(fullLine, lineStat)
            class default
              write(nout,*) "WARNING! 'boundary' keyword is only valid for 'aperiodic_ortho' box type. Ignoring."
          end select
        case("mol")
          call ReadMolCounts(lineStore, iLine, nLines, BoxArray(boxNum)%box%NMol, lineBuffer, lineStat)
        case("molmin")
          call ReadMolCounts(lineStore, iLine, nLines, BoxArray(boxNum)%box%NMolMin, lineBuffer, lineStat)
        case("molmax")
          call ReadMolCounts(lineStore, iLine, nLines, BoxArray(boxNum)%box%NMolMax, lineBuffer, lineStat)
        case default
          if(allocated(BoxArray(boxNum)%box) ) then
            call BoxArray(boxNum)%box%LoadAtomCoord(fullLine, lineStat)
            nAtoms = nAtoms + 1
          else
            write(nout,*) "ERROR! Boxed must be defined before loading atomic coorindates!"
            error stop
          endif
      end select

      IF (AllocateStat /= 0) error STOP "*** Unable to read coordinate file ***"

      ! Ensure that the called processes exited properly.
      if(lineStat .eq. -1) then
        write(0,"(A,1x,I10)") "ERROR! Parameters for command on line:", lineNumber(iLine)
        write(0, "(A)") "could not be understood. Please check command for accuracy and try again."
        write(0,*) trim(adjustl(lineStore(iLine)))
        error stop
      endif

    enddo

    deallocate(lineNumber)
    deallocate(lineStore)

    BoxArray(boxNum)%box%nAtoms = nAtoms

  end subroutine
!================================================================================
  subroutine LoadCoordFile(lineArray, nLines, lineNumber, fileName)
    use Input_Format, only: LowerCaseLine
    implicit none
    character(len=coordChunkLen), allocatable, intent(out) :: lineArray(:)
    character(len=50), intent(in) :: fileName
    integer, allocatable, intent(out) :: lineNumber(:)
    integer, intent(out) :: nLines

    character(len=maxCoordRecordLen) :: rawLine, workingLine
    character(len=50) :: modFileName
    integer :: i, iLine, lineStat, allocateStat, inOutStat, commentPos, splitPos

    modFileName = fileName
    do i = 1, len(modFileName)
      if(modFileName(i:i) == '"') modFileName(i:i) = ' '
    enddo

    open(newunit=inOutStat, file=trim(adjustl(modFileName)), status='OLD', iostat=lineStat)
    if(lineStat /= 0) then
      write(*,*) "The coordinate file could not be opened: ", trim(adjustl(modFileName))
      error stop
    endif

    nLines = 0
    do
      read(inOutStat, "(A)", iostat=lineStat) rawLine
      if(lineStat < 0) exit
      if(lineStat > 0) error stop "*** Unable to read coordinate file ***"

      commentPos = index(rawLine, "#")
      if(commentPos > 0) rawLine(commentPos:) = " "
      workingLine = adjustl(rawLine)
      do while(len_trim(workingLine) > 0)
        nLines = nLines + 1
        if(len_trim(workingLine) <= coordChunkLen) exit
        splitPos = scan(workingLine(:coordChunkLen), " ", back=.true.)
        if(splitPos <= 1) error stop "ERROR! Coordinate file contains a value longer than the input buffer."
        workingLine = adjustl(workingLine(splitPos + 1:))
      enddo
    enddo
    rewind(inOutStat)

    if(nLines == 0) error stop "ERROR! Coordinate file is empty!"

    allocate(lineArray(nLines), lineNumber(nLines), stat=allocateStat)
    if(allocateStat /= 0) error stop "*** Unable to load coordinate file ***"

    i = 0
    iLine = 0
    do
      read(inOutStat, "(A)", iostat=lineStat) rawLine
      if(lineStat < 0) exit
      if(lineStat /= 0) error stop "*** Unable to read coordinate file ***"
      iLine = iLine + 1

      commentPos = index(rawLine, "#")
      if(commentPos > 0) rawLine(commentPos:) = " "
      workingLine = adjustl(rawLine)

      do while(len_trim(workingLine) > 0)
        i = i + 1
        if(len_trim(workingLine) <= coordChunkLen) then
          lineArray(i) = workingLine
          call LowerCaseLine(lineArray(i))
          lineNumber(i) = iLine
          exit
        endif

        splitPos = scan(workingLine(:coordChunkLen), " ", back=.true.)
        if(splitPos <= 1) error stop "ERROR! Coordinate file contains a value longer than the input buffer."
        lineArray(i) = workingLine(:splitPos - 1)
        call LowerCaseLine(lineArray(i))
        lineNumber(i) = iLine
        workingLine = adjustl(workingLine(splitPos + 1:))
      enddo
    enddo
    close(inOutStat)

  end subroutine
!================================================================================
  subroutine ReadMolCounts(lineStore, firstLine, nLines, counts, lineBuffer, lineStat)
    use Input_Format, only: GetXCommand
    implicit none
    character(len=*), intent(in) :: lineStore(:)
    integer, intent(in) :: firstLine, nLines
    integer, intent(out) :: counts(:)
    integer, intent(out) :: lineBuffer, lineStat

    character(len=30) :: value
    integer :: iLine, iValue, nValues, readStat, firstValue

    lineBuffer = 0
    lineStat = 0
    nValues = 0

    do iLine = firstLine, nLines
      if(iLine == firstLine) then
        firstValue = 2
      else
        firstValue = 1
      endif

      iValue = firstValue
      do
        call GetXCommand(lineStore(iLine), value, iValue, lineStat)
        if(lineStat /= 0) exit

        nValues = nValues + 1
        read(value, *, iostat=readStat) counts(nValues)
        if(readStat /= 0) then
          lineStat = -1
          return
        endif

        if(nValues == size(counts)) then
          lineBuffer = iLine - firstLine
          lineStat = 0
          return
        endif
        iValue = iValue + 1
      enddo
    enddo

    lineStat = -1
  end subroutine
!================================================================================
end module
!================================================================================
