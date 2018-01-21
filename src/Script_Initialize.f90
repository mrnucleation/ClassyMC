!=================================================
module Input_Initialize
use VarPrecision
contains
!=================================================
  subroutine Script_Initialize
    use Common_MolInfo, only: MolData, nMolTypes, mostAtoms
    use Exceptions, only: HardError
    use ParallelVar, only: p_size, myid, nout
    use BoxData, only: BoxArray
    use TrajData, only: TrajArray
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: initSeed, sgrnd
    use MCMoveData, only: MoveProb
    implicit none
    integer :: i, j
    integer :: AllocateStat
    real(dp) :: norm

    if(initSeed < 0) then
      call system_clock(initSeed)
      initSeed = mod(initSeed, 10000)
    endif
    initSeed = p_size*initSeed + myid
    write(nout,*) "Random Generator Seed:", initSeed
    call sgrnd(initSeed)

    if(.not. allocated(BoxArray) ) then
      call HardError(-1, "simulation box")
    endif

    if(.not. allocated(EnergyCalculator) ) then
      call HardError(-1, "forcefield function(s)")
    endif

    !Normalize the move probabilities
    norm = 0E0_dp
    do i = 1, size(MoveProb)
      norm = norm + MoveProb(i)
    enddo

    do i = 1, size(MoveProb)
      MoveProb(i) = MoveProb(i)/norm
    enddo
    write(nout, *) "Move Probabilities", MoveProb(:)

    
    !Box Initialization and Safety Checks
    do i = 1, size(BoxArray)
      ! Check to see if an energy calculator has been assigned, if not default to 1 
      if( .not. associated(BoxArray(i)%box%EFunc) ) then
        BoxArray(i)%box%EFunc => EnergyCalculator(1)
      endif
      ! Initialize the neighbor lists for each box
      do j = 1, size(BoxArray(i)%box%NeighList)
        call BoxArray(i)%box%NeighList(j)%constructor(i)
      enddo

      if( allocated(BoxArray(i)%box%Constrain) ) then
        do j = 1, size(BoxArray(i)%box%Constrain)
          call BoxArray(i)%box%Constrain(j)%method%Constructor(i)
        enddo
      endif

    enddo


  end subroutine
!=================================================
end module
!=================================================
