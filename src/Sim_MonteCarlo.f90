!===========================================================================
module SimMonteCarlo
!===========================================================================
contains
!===========================================================================
  subroutine RunMonteCarlo
    use VarPrecision
    use MPI

    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use Debug, only: Debug_DumpNeiList
    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use Output_DumpCoords, only: Output_DumpData
    use ParallelVar, only: myid, ierror, nout
    use RandomGen, only: sgrnd, ListRNG
    use SimControl, only: nMoves, nCycles

    implicit none
 
    logical :: accept
    integer :: i, j
    integer :: iAtom, moveNum
    integer(kind=8) :: iCycle, iMove
    real(dp) :: E_T, E_Final
    character(len=50) :: fileName

    do i = 1, size(BoxArray)
      call BoxArray(i) % box % ComputeEnergy
      call BoxArray(i) % box % NeighList(1) % BuildList
    enddo


    call Trajectory( int(0,kind=8), int(0,kind=8) )
    write(nout, *) "============================================"
    write(nout, *) "       Simulation Start!"
    write(nout, *) "============================================"

    !-------Main Monte Carlo Simulation Loop-------
    do iCycle = 1, nCycles

      !-----Start Move Loop
      do iMove = 1, nMoves
        moveNum = ListRNG(MoveProb)
        call Moves(moveNum) % Move % FullMove(BoxArray(1)%box, accept)
        call Analyze(iCycle, iMove, accept, .true.)
!        call Debug_DumpNeiList(1, 1, 1)
      enddo 
      !------End Move Loop
      if(mod(iCycle, 1000) == 0) then
        write(*,*) iCycle, BoxArray(1)%box%ETotal, (Moves(j)%Move%GetAcceptRate(), j=1, size(Moves))
      endif

      if(mod(iCycle, 100) == 0) then
        call BoxArray(1) % box % NeighList(1) % BuildList
      endif

!      if(mod(iCycle, 10) == 0) then
!        call Moves(1) % Move % Maintenance
!      endif

      call Analyze(iCycle, iMove, accept, .false.)
      call Trajectory(iCycle, iMove)
    enddo
    !-------End of Main Monte Carlo Simulation Loop-------
    
    call Debug_DumpNeiList(1, 1, 1)
    E_Final = BoxArray(1)%box%ETotal
    do i = 1, size(BoxArray)
      call BoxArray(i) % box % ComputeEnergy
      call BoxArray(i) % box % NeighList(1) % BuildList
    enddo
    call Debug_DumpNeiList(1, 1, 1)

    write(nout, *) "Culmative Energy:", E_Final
    write(nout, *) "Final Energy:",  BoxArray(1)%box%ETotal
    write(nout, *) "Difference:",  E_Final - BoxArray(1)%box%ETotal

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Finalize
      enddo
    endif

    call Output_DumpData
      
  end subroutine
!===========================================================================
  subroutine Analyze(iCycle, iMove, accept, moveloop)
    use AnalysisData, only: AnalysisArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    logical, intent(in) :: accept
    logical, intent(in) :: moveloop
    integer :: iAn

    if( allocated(AnalysisArray) ) then
      do iAn = 1, size(AnalysisArray)
        if( AnalysisArray(iAn)%func%perMove .eqv. moveloop) then
          if(mod(iCycle, AnalysisArray(iAn)%func%UpdateFreq) == 0) then
            call AnalysisArray(iAn) % func % Compute(accept)
          endif
        endif
      enddo
    endif
 
  end subroutine
!===========================================================================
  subroutine Maintenance(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i

    if(mod(iCycle, Sampling%maintFreq) == 0 ) then
      call Sampling % Maintenance
    endif

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        if(mod(iCycle, AnalysisArray(i)%func%maintFreq) == 0) then
          call AnalysisArray(i) % func % Maintenance
        endif
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        if(mod(iCycle, TrajArray(i)%traj%maintFreq) == 0) then
          call TrajArray(i) % traj % Maintenance
        endif
      enddo
    endif

    do i = 1, size(BoxArray)
      if(mod(iCycle, BoxArray(i)%box%maintFreq) == 0) then
        call BoxArray(i) % box % Maintenance
      endif
    enddo


    do i = 1, size(Moves)
      if(mod(iCycle, Moves(i)%move%maintFreq) == 0) then
        call Moves(i) % move % Maintenance
      endif
    enddo

  end subroutine
!===========================================================================
  subroutine Trajectory(iCycle, iMove)
    use TrajData, only: TrajArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: iTraj

    if( allocated(TrajArray) ) then
      do iTraj = 1, size(TrajArray)
        if(mod(iCycle, TrajArray(iTraj)%traj%outfreq) == 0) then
          call TrajArray(iTraj) % traj % WriteFrame
        endif
      enddo
    endif

  end subroutine
!===========================================================================
end module
!===========================================================================
