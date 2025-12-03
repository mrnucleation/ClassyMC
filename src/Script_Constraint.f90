!================================================================================
module Input_Constraint
  use Input_Format
  use VarPrecision
!================================================================================
contains
!================================================================================
  subroutine Script_Constraint(lineStore, iLine, BoxNum, lineBuffer, lineStat)
    use Input_Format, only: GetXCommand
    use BoxData, only: BoxArray
    use ParallelVar, only: nout
    use Constrain_DistanceCriteria, only: DistCriteria
    use Constrain_EnergyCeiling, only: EnergyCeiling
    use Constrain_EnergyFloor, only: EnergyFloor
    use Constrain_FreezeType, only: FreezeType
!    use Constrain_HardWall, only: HardWall
    use Constrain_MolTotal, only: MolTotal
    use Constrain_MultiAtomDistanceCriteria, only: MultiAtomDistCrit

    implicit none
    character(len=maxLineLen), intent(in) :: linestore(:) 
    integer, intent(in) :: BoxNum, iLine
    integer, intent(out) :: lineStat, lineBuffer

    character(len=30) :: command
    integer :: i, intVal, curLine
    integer :: nItems, AllocateStat
    real(dp) :: realVal

    lineStat  = 0
    call FindCommandBlock(iLine, lineStore, "end_create", lineBuffer)
    nItems = lineBuffer - 1

    allocate(BoxArray(BoxNum)%box%Constrain(1:nItems), stat = AllocateStat)
    !Safety check to ensure that the index number is within proper bounds
    do i = 1, nItems
      curLine = iLine + i
      call GetXCommand(lineStore(curLine), command, 1, lineStat)
      select case(trim(adjustl(command)))
        case("distancecriteria")
          allocate( DistCriteria::BoxArray(BoxNum)%box%Constrain(i)%method )

        case("multidistancecriteria")
          allocate( MultiAtomDistCrit::BoxArray(BoxNum)%box%Constrain(i)%method )

        case("energyceiling")
          allocate( EnergyCeiling ::BoxArray(BoxNum)%box%Constrain(i)%method )


        case("energyfloor")
          allocate( EnergyFloor ::BoxArray(BoxNum)%box%Constrain(i)%method )

!        case("hardwall")
!          allocate( HardWall ::BoxArray(BoxNum)%box%Constrain(i)%method )

        case("freezetype")
          allocate( FreezeType ::BoxArray(BoxNum)%box%Constrain(i)%method )

        case("moltotal")
          allocate( MolTotal::BoxArray(BoxNum)%box%Constrain(i)%method )


        case default
          lineStat = -1
          return
      end select
      call BoxArray(BoxNum)%box%Constrain(i)%method%ProcessIO(lineStore(curLine), lineStat)
    enddo   

  end subroutine
!================================================================================
end module
!================================================================================
