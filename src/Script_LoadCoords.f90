!================================================================================
module Input_LoadCoords
  use VarPrecision
  contains
!================================================================================
  subroutine Script_ReadCoordFile(filename, boxNum, lineStat)
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nMolTypes
    use ForcefieldData, only: nForceFields
    use Input_Format, only: maxLineLen, LoadFile, GetXCommand, ReplaceText
    use Input_SimBoxes, only: Script_BoxType
    use ParallelVar, only: nout, myid
    implicit none
    character(len=50), intent(inout) :: fileName      
    integer, intent(in) :: boxNum
    integer, intent(out) :: lineStat

    integer :: iCharacter
    integer :: i, j, nLines,AllocateStat, iLine, lineBuffer
    integer, allocatable :: lineNumber(:)
    character(len=maxLineLen), allocatable :: lineStore(:)

    character(len=30) :: dummy, command, val
    character(len=maxLineLen) :: newline
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
    call LoadFile(lineStore, nLines, lineNumber, fileName)
    lineBuffer = 0
    nAtoms = 0
    do iLine = 1, nLines
      if(lineBuffer > 0) then
        lineBuffer = lineBuffer - 1
        cycle
      endif
      lineStat = 0
      call GetXCommand(lineStore(iLine), command, 1, lineStat)
      if(lineStat .eq. 1) then
        cycle
      endif

!      write(nout,*) trim(adjustl(lineStore(iLine)))
      select case(trim(adjustl( command )))
        case("boxtype")
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          newline = ""
          newline(1:30) = val(1:30)
          call Script_BoxType(newline, boxNum, lineStat)
        case("dimension")
          call BoxArray(boxNum)%box%LoadDimension(lineStore(iLine), lineStat)
        case("mol")
          read(lineStore(iLine), *) dummy, (BoxArray(boxNum)%box%NMol(j), j=1,nMolTypes) 
        case("molmin")
          read(lineStore(iLine), *) dummy, (BoxArray(boxNum)%box%NMolMin(j), j=1,nMolTypes) 
        case("molmax")
          read(lineStore(iLine), *) dummy, (BoxArray(boxNum)%box%NMolMax(j), j=1,nMolTypes) 
        case default
          if(allocated(BoxArray(boxNum)%box) ) then
            call BoxArray(boxNum)%box%LoadAtomCoord(lineStore(iLine), lineStat)
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
end module
!================================================================================
