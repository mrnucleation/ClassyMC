!================================================================================
module Input_Forcefield
  use VarPrecision
  use Input_Format
  use ForcefieldData

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
  subroutine Script_FieldType(line, FFNum, lineStat)
    use ForcefieldData, only: nForceFields
    use FF_Pair_LJ_Cut, only: Pair_LJ_Cut
    use FF_Pair_LJ_Cut_NoNei, only: Pair_LJ_Cut_NoNei
    use FF_Pair_Tersoff, only: Pair_Tersoff
    use ParallelVar, only: nout
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in) :: FFNum
    integer, intent(out) :: lineStat

    character(len=30) :: dummy, command, FF_Type
    character(len=30) :: fileName      
    logical :: logicValue
    integer :: j
    real(dp) :: realValue

    lineStat  = 0
    read(line, *) FF_Type

    !Safety check to ensure that the index number is within proper bounds
    call LowerCaseLine(FF_Type)
    select case(trim(adjustl(FF_Type)))
      case("lj_cut")
        allocate(Pair_LJ_Cut::EnergyCalculator(FFNum) % Method)
        write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as LJ_Cut style"

      case("lj_cut_nonei")
        allocate(Pair_LJ_Cut_NoNei::EnergyCalculator(FFNum) % Method)
        write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as LJ_Cut (No Neighbor List) style"

      case("tersoff")
        allocate(Pair_Tersoff::EnergyCalculator(FFNum) % Method)
        write(nout,"(A,I2,A)") "Forcefield", FFNum, " allocated as Tersoff style"

      case default
        write(*,*) "Here"
        lineStat = -1

      end select

  end subroutine
!================================================================================
  subroutine Script_ReadFieldFile(filename, lineStat)
    use ForcefieldData, only: nForceFields
    use ParallelVar, only: nout
    use Input_Format, only: maxLineLen, LoadFile
    implicit none
    character(len=30), intent(in) :: fileName      
    integer, intent(out) :: lineStat

    integer :: nLines
    integer, allocatable :: lineNumber(:)
    character(len=maxLineLen), allocatable :: lineStore(:)

    character(len=30) :: command, val
    logical :: logicValue
    integer :: intVal
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

        case("forcefieldtype")]
          call Script_Forcefield(lineStore(iLine), lineStat)

        case("molecule")
          call FindCommandBlock(iLine, lineStore, "end_molecule", lineBuffer)
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue

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
    implicit none
    character(len=maxLineLen), intent(in) :: cmdBlock(:)
    integer, intent(in) :: molType
    integer, intent(out) :: lineStat
    integer :: j, intValue(1:5)

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
           if( .not. allocated(EnergyCalculator) ) then
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

  end subroutine
!================================================================================
  subroutine Script_SetUnits(line, lineStat)
    use Units, only: FindEngUnit, FindLengthUnit, FindAngularUnit
    implicit none
    character(len=maxLineLen), intent(in) :: line
    integer, intent(out) :: lineStat

    character(len=30) :: command, unitType
    character(len=30) :: fileName      
    integer :: FFNum

    lineStat  = 0
    call LowerCaseLine(line)
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
