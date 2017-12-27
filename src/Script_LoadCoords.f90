!================================================================================
module Input_LoadCoords
  use VarPrecision
  contains
!================================================================================
  subroutine Script_ReadCoordFile(filename, boxNum, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Input_Format, only: maxLineLen, LoadFile, GetXCommand
    use Input_SimBoxes, only: Script_BoxType
    use Common_MolInfo, only: nMolTypes
    use BoxData, only: BoxArray
    implicit none
    character(len=50), intent(in) :: fileName      
    integer, intent(in) :: boxNum
    integer, intent(out) :: lineStat

    integer :: i, j, nLines, AllocateStat, iLine, lineBuffer
    integer, allocatable :: lineNumber(:)
    character(len=maxLineLen), allocatable :: lineStore(:)

    character(len=30) :: dummy, command, val
    character(len=maxLineLen) :: newline
    integer :: nAtoms

    lineStat  = 0
    write(*,*) "Loading file ", trim(adjustl(fileName))
    call LoadFile(lineStore, nLines, lineNumber, fileName)
    lineBuffer = 0
    nAtoms = 0
    do iLine = 1, nLines
      if(lineBuffer .gt. 0) then
        lineBuffer = lineBuffer - 1
        cycle
      endif
      lineStat = 0
      call GetXCommand(lineStore(iLine), command, 1, lineStat)
      if(lineStat .eq. 1) then
        cycle
      endif

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
          call BoxArray(boxNum)%box%LoadAtomCoord(lineStore(iLine), lineStat)
          nAtoms = nAtoms + 1
      end select

      IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

      ! Ensure that the called processes exited properly.
      if(lineStat .eq. -1) then
        write(*,"(A,1x,I10)") "ERROR! Parameters for command on line:", lineNumber(iLine)
        write(*, "(A)") "could not be understood. Please check command for accuracy and try again."
        write(*,*) trim(adjustl(lineStore(iLine)))
        stop
      endif

    enddo

    deallocate(lineNumber)
    deallocate(lineStore)

    BoxArray(boxNum)%box%nAtoms = nAtoms

  end subroutine
!================================================================================
end module
!================================================================================
