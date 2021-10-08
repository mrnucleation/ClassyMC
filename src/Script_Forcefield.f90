!================================================================================
module Input_Forcefield
  use VarPrecision
  use ForcefieldData
  use Input_FieldType
  use Input_AngleType
  use Input_BondType
  use Input_TorsionType
  use Input_Format
  use Input_MiscType

  real(dp) :: engUnit = 1E0_dp 
  real(dp) :: lenUnit = 1E0_dp
  real(dp) :: angUnit = 1E0_dp
  contains
!================================================================================
  subroutine Script_ReadFieldFile(filename, lineStat)
    use ForcefieldData, only: nForceFields, EnergyCalculator
    use ParallelVar, only: nout
    use Input_Format, only: maxLineLen, LoadFile
    use Common_MolInfo, only: AtomData, MolData, BondData, AngleData, TorsionData, &
                             nMolTypes, nAtomTypes, nBondTypes, nAngleTypes, &
                            nTorsionTypes
    implicit none
    character(len=50), intent(in) :: fileName      
    integer, intent(out) :: lineStat

    integer :: i, nLines, nItems, curLine, AllocateStat, iLine, lineBuffer
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
!      write(*,*) trim(adjustl( lineStore(iLine) ))
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
!        -----------------------------------------------------------------------------
        case("forcefield")
          call FindCommandBlock(iLine, lineStore, "end_forcefield", lineBuffer)
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue
          if( .not. allocated(EnergyCalculator) ) then
            write(0,*) "ERROR! The forcefield type must be defined before parameters can be modifed."
            error stop
          endif
          nItems = lineBuffer - 1
          do i = 1, nItems
            curLine = iLine + i
            call EnergyCalculator(intValue)%Method%ProcessIO(lineStore(curLine))
          enddo   
!        -----------------------------------------------------------------------------
        case("forcefieldtype")
          call FindCommandBlock(iLine, lineStore, "end_forcefieldtype", lineBuffer)
          nItems = lineBuffer - 1
          if( .not. allocated(EnergyCalculator) ) then
            nForceFields = nItems
            allocate(EnergyCalculator(1:nItems), stat = AllocateStat)
            do i = 1, nItems
              curLine = iLine + i
              call Script_FieldType(linestore(curLine), i, lineStat)
              if(lineStat < 0) then
                write(0,*) "ERROR! Unknown forcefield type."
                error stop
              endif
              call EnergyCalculator(i)%Method%Constructor
            enddo   
          else
            write(0,*) "ERROR! The forcefieldtype has already been used and can not be called twice"
            error stop
          endif

!        -----------------------------------------------------------------------------
        case("atomdef")
          call FindCommandBlock(iLine, lineStore, "end_atomdef", lineBuffer)
          nItems = lineBuffer - 1
          if( .not. allocated(AtomData) ) then
            nAtomTypes = nItems
            allocate(AtomData(1:nItems), stat = AllocateStat)
            do i = 1, nItems
              curLine = iLine + i
              read(lineStore(curLine), *) AtomData(i)%symb, AtomData(i)%mass
            enddo   
          else
            write(0,*) "ERROR! The atomdef has already been used and can not be called twice"
            error stop
          endif

!        -----------------------------------------------------------------------------
        case("bonddef")
          call FindCommandBlock(iLine, lineStore, "end_bonddef", lineBuffer)
          nItems = lineBuffer - 1
          if( .not. allocated(BondData) ) then
            nBondTypes = nItems
            allocate(BondData(1:nItems), stat = AllocateStat)
            do i = 1, nItems
              curLine = iLine + i
!              read(lineStore(curLine), *) BondData(i)%rEq
              call Script_BondType(lineStore(curLine), i, lineStat)
!              write(*,*) lineStat
            enddo   
          else
            write(0,*) "ERROR! The BondDef has already been used and can not be called twice"
            error stop
          endif
!        -----------------------------------------------------------------------------
        case("angledef")
          call FindCommandBlock(iLine, lineStore, "end_angledef", lineBuffer)
          nItems = lineBuffer - 1
          if( .not. allocated(AngleData) ) then
            nAngleTypes = nItems
            allocate(AngleData(1:nItems), stat = AllocateStat)
            do i = 1, nItems
              curLine = iLine + i
!              read(lineStore(curLine), *) AngleData(i)%rEq
              call Script_AngleType(lineStore(curLine), i, lineStat)
!              write(*,*) lineStat
            enddo   
          else
            write(0,*) "ERROR! The AngleDef has already been used and can not be called twice"
            error stop
          endif

!        -----------------------------------------------------------------------------
        case("torsiondef")
          call FindCommandBlock(iLine, lineStore, "end_torsiondef", lineBuffer)
          nItems = lineBuffer - 1
          if( .not. allocated(TorsionData) ) then
            nTorsionTypes = nItems
            allocate(TorsionData(1:nItems), stat = AllocateStat)
            do i = 1, nItems
              curLine = iLine + i
              call Script_TorsionType(lineStore(curLine), i, lineStat)
            enddo   
          else
            write(0,*) "ERROR! The TorsionDef has already been used and can not be called twice"
            error stop
          endif


!        -----------------------------------------------------------------------------
        case("molecule")
          call FindCommandBlock(iLine, lineStore, "end_molecule", lineBuffer)
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue
          if(intValue > nMolTypes) then
            write(0,*) "ERROR! Molecule index out of bounds in the forcefield file!"
            write(0,*) "Index Called: ", intValue
            error stop
          endif
          call Script_ReadMolDef( lineStore(iLine+1:iLine+lineBuffer-1), intValue, lineStat )

!        -----------------------------------------------------------------------------
        case("moleculetypes")
          call GetXCommand(lineStore(iLine), val, 2, lineStat)
          read(val, *) intValue
          if(allocated(MolData)) then
            lineStat = -1
            return
          endif
          nMolTypes = intValue
          allocate(MolData(1:nMolTypes), stat = AllocateStat)

!        -----------------------------------------------------------------------------
        case("units")
          call Script_SetUnits(lineStore(iLine), lineStat)

!        -----------------------------------------------------------------------------
        case default
          write(0,"(A,2x,I10)") "ERROR! Unknown Command on Line", lineNumber(iLine)
          write(0,*) trim(adjustl(lineStore(iLine)))
          error stop
      end select

      IF (AllocateStat /= 0) error STOP "*** Unable to read forcefield file ***"

      ! Ensure that the called processes exited properly.
      if(lineStat .eq. -1) then
        write(0,"(A,1x,I10)") "ERROR! Parameters for command on line:", lineNumber(iLine)
        write(0, "(A)") " could not be understood. Please check command for accuracy and try again."
        write(0,*) trim(adjustl(lineStore(iLine)))
        error stop
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
    use Input_RegrowType, only: Script_RegrowType
    use Common_MolInfo, only: MolData, mostAtoms
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
!        -----------------------------------------------------------------------------
        case("regrowthtype")
             call Script_RegrowType(cmdBlock(iLine), MolType, lineStat)
!        -----------------------------------------------------------------------------
        case("atoms")
           if( .not. allocated(MolData(molType)%atomType) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_atoms", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nAtoms = nItems
             if(nItems > mostAtoms) then
               mostAtoms = nItems
             endif
             allocate(MolData(molType)%atomType(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) intValue(1)
               MolData(molType)%atomType(i) = intValue(1)
             enddo   
           else
             write(0,*) "ERROR! The create energycalculators command has already been used and can not be called twice"
             error stop
           endif

!        -----------------------------------------------------------------------------
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
             write(0,*) "ERROR! The bonds for this molecule has already been defined."
             error stop
           endif 

!        -----------------------------------------------------------------------------
        case("angles")
           if( .not. allocated(MolData(molType)%angle) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_angles", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nAngles = nItems
             allocate(MolData(molType)%angle(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,4)
               MolData(molType)%angle(i)%angleType = intValue(1)
               MolData(molType)%angle(i)%mem1 = intValue(2)
               MolData(molType)%angle(i)%mem2 = intValue(3)
               MolData(molType)%angle(i)%mem3 = intValue(4)
             enddo   
           else
             write(0,*) "ERROR! The create angle command has already been used and can not be called twice"
             error stop
           endif
!        -----------------------------------------------------------------------------
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
             write(0,*) "ERROR! The create torsion command has already been used and can not be called twice"
             error stop
           endif
!        -----------------------------------------------------------------------------
        case("misc")
           if( .not. allocated(MolData(molType)%Miscdata) ) then
             call FindCommandBlock(iLine, cmdBlock, "end_misc", lineBuffer)
             nItems = lineBuffer - 1
             MolData(molType)%nMisc = nItems
             allocate(MolData(molType)%miscdata(1:nItems), stat = AllocateStat)
             do i = 1, nItems
               curLine = iLine + i
               call Script_MiscType(cmdBlock(curLine), molType, i, lineStat)
             enddo   
           else
             write(0,*) "ERROR! The create misc command has already been used and can not be called twice"
             error stop
           endif
!        -----------------------------------------------------------------------------
        case default
          write(0,"(A,2x,I10)") "ERROR! Unknown Command on Line"
          write(0,*) trim(adjustl(cmdBlock(iLine)))
          error stop
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
        inEngUnit = FindEngUnit(unitType)
      case("angle")
        call GetXCommand(line, unitType, 3, lineStat)
        inAngUnit = FindAngularUnit(unitType)
      case("length")
        call GetXCommand(line, unitType, 3, lineStat)
        inLenUnit = FindLengthUnit(unitType)
      case default
        lineStat = -1
    end select


  end subroutine


!================================================================================
end module
!================================================================================
