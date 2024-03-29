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
          call Script_BuildGraph(intValue)

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
    integer :: atm1, atm2, atm3, atm4
    integer :: nIntra

    lineStat  = 0
    lineBuffer = 0
    nIntra = 0
    atm1 = 0
    atm2 = 0
    atm3 = 0
    atm4 = 0
    nItems = 0
    AllocateStat = 0
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
             allocate(MolData(molType)%nAtmBonds(1:MolData(molType)%nAtoms), stat = AllocateStat)
             allocate(MolData(molType)%atmBonds(1:nItems, 1:MolData(molType)%nAtoms), stat = AllocateStat)
             MolData(molType)%nAtmBonds = 0
             MolData(molType)%atmBonds = 0
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,3)
               MolData(molType)%bond(i)%bondType = intValue(1)
               atm1 = intValue(2)
               atm2 = intValue(3)
               MolData(molType)%bond(i)%mem1 = atm1
               MolData(molType)%bond(i)%mem2 = atm2

               nIntra = MolData(molType)%nAtmBonds(atm1) + 1
               MolData(molType)%nAtmBonds(atm1) = nIntra
               MolData(molType)%atmBonds(nIntra, atm1) = i

               nIntra = MolData(molType)%nAtmBonds(atm2) + 1
               MolData(molType)%nAtmBonds(atm2) = nIntra
               MolData(molType)%atmBonds(nIntra, atm2) = i
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
             allocate(MolData(molType)%nAtmAngles(1:MolData(molType)%nAtoms), stat = AllocateStat)
             allocate(MolData(molType)%atmAngles(1:nItems, 1:MolData(molType)%nAtoms), stat = AllocateStat)
             MolData(molType)%nAtmAngles = 0
             MolData(molType)%atmAngles = 0
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,4)
               atm1 = intValue(2)
               atm2 = intValue(3)
               atm3 = intValue(4)

               MolData(molType)%angle(i)%angleType = intValue(1)
               MolData(molType)%angle(i)%mem1 = atm1
               MolData(molType)%angle(i)%mem2 = atm2
               MolData(molType)%angle(i)%mem3 = atm3

               nIntra = MolData(molType)%nAtmAngles(atm1) + 1
               MolData(molType)%nAtmAngles(atm1) = nIntra
               MolData(molType)%atmAngles(nIntra, atm1) = i

               nIntra = MolData(molType)%nAtmAngles(atm2) + 1
               MolData(molType)%nAtmAngles(atm2) = nIntra
               MolData(molType)%atmAngles(nIntra, atm2) = i

               nIntra = MolData(molType)%nAtmAngles(atm3) + 1
               MolData(molType)%nAtmAngles(atm3) = nIntra
               MolData(molType)%atmAngles(nIntra, atm3) = i
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
             allocate(MolData(molType)%nAtmTorsions(1:MolData(molType)%nAtoms), stat = AllocateStat)
             allocate(MolData(molType)%atmTorsions(1:nItems, 1:MolData(molType)%nAtoms), stat = AllocateStat)
             MolData(molType)%nAtmTorsions = 0
             MolData(molType)%atmTorsions = 0 
             do i = 1, nItems
               curLine = iLine + i
               read(cmdBlock(curLine),*) (intValue(j), j=1,5)
               MolData(molType)%torsion(i)%torsType = intValue(1)
               atm1 = intValue(2)
               atm2 = intValue(3)
               atm3 = intValue(4)
               atm4 = intValue(5)
               MolData(molType)%torsion(i)%mem1 = atm1
               MolData(molType)%torsion(i)%mem2 = atm2
               MolData(molType)%torsion(i)%mem3 = atm3 
               MolData(molType)%torsion(i)%mem4 = atm4
               nIntra = MolData(molType)%nAtmTorsions(atm1) + 1
               MolData(molType)%nAtmTorsions(atm1) = nIntra
               MolData(molType)%atmTorsions(nIntra, atm1) = i

               nIntra = MolData(molType)%nAtmTorsions(atm2) + 1
               MolData(molType)%nAtmTorsions(atm2) = nIntra
               MolData(molType)%atmTorsions(nIntra, atm2) = i

               nIntra = MolData(molType)%nAtmTorsions(atm3) + 1
               MolData(molType)%nAtmTorsions(atm3) = nIntra
               MolData(molType)%atmTorsions(nIntra, atm3) = i

               nIntra = MolData(molType)%nAtmTorsions(atm4) + 1
               MolData(molType)%nAtmTorsions(atm4) = nIntra
               MolData(molType)%atmTorsions(nIntra, atm4) = i
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
  subroutine Script_BuildGraph(molnum)
    use Common_MolInfo, only: MolData
    implicit none
    integer, intent(in) :: molnum
    integer :: iGraph, iEdge
    integer :: iBond, iAngle, iTorsion
    integer :: nAtoms, intraType
    integer :: mem1, mem2, mem3, mem4

    nAtoms = MolData(molnum)%nAtoms

    call MolData(molnum)%molgraph%Constructor(nAtoms)

    do iBond = 1, MolData(molnum)%nBonds
      intraType = MolData(molnum)%bond(iBond)%bondType
      mem1 = MolData(molnum)%bond(iBond)%mem1
      mem2 = MolData(molnum)%bond(iBond)%mem2
      call MolData(molnum)%molgraph%AddEdge(mem1, mem2, intraType)
    enddo

    do iAngle = 1, MolData(molnum)%nAngles
      intraType = MolData(molnum)%angle(iAngle)%angleType
      mem1 = MolData(molnum)%angle(iAngle)%mem1
      mem2 = MolData(molnum)%angle(iAngle)%mem2
      mem3 = MolData(molnum)%angle(iAngle)%mem3
      call MolData(molnum)%molgraph%AddAngle(mem1, mem2, mem3, intraType)
    enddo

    do iTorsion = 1, MolData(molnum)%nTors
      intraType = MolData(molnum)%torsion(iTorsion)%torsType
      mem1 = MolData(molnum)%torsion(iTorsion)%mem1
      mem2 = MolData(molnum)%torsion(iTorsion)%mem2
      mem3 = MolData(molnum)%torsion(iTorsion)%mem3
      mem4 = MolData(molnum)%torsion(iTorsion)%mem4
      call MolData(molnum)%molgraph%AddTorsion(mem1, mem2, mem3, mem4, intraType)
    enddo



!    write(*,*) "Edges"
!    do iGraph = 1, size(MolData(molnum)%molgraph%nodes)
!      if(.not. allocated(MolData(molnum)%molgraph%nodes(iGraph)%edge) ) cycle
!      do iEdge = 1, size(MolData(molnum)%molgraph%nodes(iGraph)%edge)
!        write(*,*) iGraph, iEdge, MolData(molnum)%molgraph%nodes(iGraph)%edge(iEdge)
!      enddo
!    enddo
!    write(*,*)

!    write(*,*) "EndAngle"
!    do iGraph = 1, size(MolData(molnum)%molgraph%nodes)
!      if(.not. allocated(MolData(molnum)%molgraph%nodes(iGraph)%endangle) ) cycle
!      do iEdge = 1, size(MolData(molnum)%molgraph%nodes(iGraph)%endangle, 2)
!        write(*,*) iGraph, iEdge, MolData(molnum)%molgraph%nodes(iGraph)%endangle(1:2,iEdge) 
!      enddo
!    enddo

!    write(*,*) "MidAngle"
!    do iGraph = 1, size(MolData(molnum)%molgraph%nodes)
!      if(.not. allocated(MolData(molnum)%molgraph%nodes(iGraph)%midangle) ) cycle
!      do iEdge = 1, size(MolData(molnum)%molgraph%nodes(iGraph)%midangle, 2)
!        write(*,*) iGraph, iEdge, MolData(molnum)%molgraph%nodes(iGraph)%midangle(1:2,iEdge) 
!      enddo
!    enddo

!    write(*,*) "EndTorsion"
!    do iGraph = 1, size(MolData(molnum)%molgraph%nodes)
!      if(.not. allocated(MolData(molnum)%molgraph%nodes(iGraph)%endtorsion) ) cycle
!      do iEdge = 1, size(MolData(molnum)%molgraph%nodes(iGraph)%endtorsion, 2)
!        write(*,*) iGraph, iEdge, MolData(molnum)%molgraph%nodes(iGraph)%endtorsion(1:3,iEdge) 
!      enddo
!    enddo

!    write(*,*) MolData(molnum)%molgraph%IsConnected()



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
      case("angle")
        call GetXCommand(line, unitType, 3, lineStat)
        inAngUnit = FindAngularUnit(unitType)
      case("energy")
        call GetXCommand(line, unitType, 3, lineStat)
        inEngUnit = FindEngUnit(unitType)
      case("distance")
        call GetXCommand(line, unitType, 3, lineStat)
        inLenUnit = FindLengthUnit(unitType)
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
