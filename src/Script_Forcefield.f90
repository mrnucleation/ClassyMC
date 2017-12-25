!================================================================================
module Input_Forcefield
  use VarPrecision
  use Input_Format
  use ForcefieldData
  use Input_FieldType

  real(dp) :: engUnit = 1E0_dp 
  real(dp) :: lenUnit = 1E0_dp
  real(dp) :: angUnit = 1E0_dp
  contains
!================================================================================
  subroutine Script_Forcefield(line, lineStat)
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command 
    character(len=30) :: fileName      
    logical :: logicValue
    integer :: j
    integer :: intValue
    integer :: FFNum
    real(dp) :: realValue

    lineStat  = 0
    read(line, *) dummy, FFNum, command
    call LowerCaseLine(command)
    select case(trim(adjustl(command)))
      case("readfile")
        read(line,*) dummy, FFNum, command, fileName
!          call EnergyCalculator(FFNum) % Method % ReadParFile(fileName)
        case default
          lineStat = -1
      end select


  end subroutine

!================================================================================
  subroutine Script_ReadFieldFile(filename, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Input_Format, only: maxLineLen, LoadFile
    use Common_MolInfo, only: MolData, nMolTypes
    implicit none
    character(len=50), intent(in) :: fileName      
    integer, intent(out) :: lineStat

    integer :: nLines, AllocateStat, iLine, lineBuffer
    integer, allocatable :: lineNumber(:)
    character(len=maxLineLen), allocatable :: lineStore(:)

    character(len=30) :: command, val
    logical :: logicValue
    integer :: intValue
    real(dp) :: realValue

    lineStat  = 0
    call LoadFile(lineStore, nLines, lineNumber, fileName)
    lineBuffer = 0
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
      call LowerCaseLine(command)

      select case(trim(adjustl( command )))
        case("forcefield")
          call FindCommandBlock(iLine, lineStore, "end_forcefield", lineBuffer)
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue

        case("forcefieldtype")
          call Script_Forcefield(lineStore(iLine), lineStat)

        case("molecule")
          call FindCommandBlock(iLine, lineStore, "end_molecule", lineBuffer)
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue
          call Script_ReadMolDef( lineStore(iLine+1:iLine+lineBuffer-1), intValue, lineStat )

        case("moleculetypes")
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue
          if(allocated(MolData)) then
            lineStat = -1
            return
          endif
          nMolTypes = intValue
          allocate(MolData(1:nMolTypes), stat = AllocateStat)

        case("units")
          call Script_SetUnits(lineStore(iLine), lineStat)

        case default
          write(*,"(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
          write(*,*) trim(adjustl(lineStore(iLine)))
          stop
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

  end subroutine
!================================================================================
  subroutine Script_ReadMolDef(cmdBlock, molType, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Input_Format, only: maxLineLen, LoadFile
    use Common_MolInfo, only: MolData
    implicit none
    character(len=maxLineLen), intent(in) :: cmdBlock(:)
    integer, intent(in) :: molType
    integer, intent(out) :: lineStat
    integer :: i, j, intValue(1:5), iLine, lineBuffer, nLines, curLine
    integer :: AllocateStat, nItems
    character(len=30) :: command

    lineStat  = 0
    lineBuffer = 0
    nLines = size(cmdBlock)
    do iLine = 1, nLines
      if(lineBuffer .gt. 0) then
        lineBuffer = lineBuffer - 1
        cycle
      endif
      lineStat = 0
      call GetXCommand(cmdBlock(iLine), command, 1, lineStat)
      if(lineStat .eq. 1) then
        cycle
      endif
      call LowerCaseLine(command)

      select case(trim(adjustl( command )))
        case("atoms")
           if( .not. allocated(MolData(molType)%atomType) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_atoms", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nAtoms = nItems
             allocate(MolData(molType)%atomType(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) intValue(1)
               MolData(molType)%atomType(i) = intValue(1)
             enddo   
           else
             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
             stop
           endif

        case("bonds")
           if( .not. allocated(MolData(molType)%bond) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_bonds", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nBonds = nItems
             allocate(MolData(molType)%bond(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,3)
               MolData(molType)%bond(i)%bondType = intValue(1)
               MolData(molType)%bond(i)%mem1 = intValue(2)
               MolData(molType)%bond(i)%mem2 = intValue(3)
             enddo   
           else
             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
             stop
           endif 

        case("angles")
           if( .not. allocated(MolData(molType)%angle) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_angles", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nAngles = nItems
             allocate(MolData(molType)%bond(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,4)
               MolData(molType)%angle(i)%angType = intValue(1)
               MolData(molType)%angle(i)%mem1 = intValue(2)
               MolData(molType)%angle(i)%mem2 = intValue(3)
               MolData(molType)%angle(i)%mem3 = intValue(4)
             enddo   
           else
             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
             stop
           endif
        case("torsion")
           if( .not. allocated(MolData(molType)%torsion) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_torsion", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nTors = nItems
             allocate(MolData(molType)%torsion(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,5)
               MolData(molType)%torsion(i)%torsType = intValue(1)
               MolData(molType)%torsion(i)%mem1 = intValue(2)
               MolData(molType)%torsion(i)%mem2 = intValue(3)
               MolData(molType)%torsion(i)%mem3 = intValue(4)
               MolData(molType)%torsion(i)%mem4 = intValue(5)
             enddo   
           else
             write(*,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
             stop
           endif

        case default
          write(*,"(A,2x,I10)") "ERROR! Unknown Command on Line"
          write(*,*) trim(adjustl(cmdBlock(iLine)))
          stop
      end select

    enddo

  end subroutine
!================================================================================
  subroutine Script_SetUnits(line, lineStat)
    use Units
    implicit none
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: command, unitType
    character(len=30) :: fileName      
    integer :: FFNum

    lineStat  = 0
    call GetXCommand(line, command, 2, lineStat)
    select case(trim(adjustl(command)))
      case("energy")
        call GetXCommand(line, unitType, 3, lineStat)
        engUnit = FindEngUnit(unitType)
      case("angle")
        call GetXCommand(line, unitType, 3, lineStat)
        angUnit = FindAngularUnit(unitType)
      case("length")
        call GetXCommand(line, unitType, 3, lineStat)
        lenUnit = FindLengthUnit(unitType)
      case default
        lineStat = -1
    end select


  end subroutine


!================================================================================
end module
!================================================================================
