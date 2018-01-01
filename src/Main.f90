!===========================================================================
  program Classy
    use VarPrecision
    use MPI
    use SimControl, only: nMoves, nCycles
    use ParallelVar, only: myid, p_size, ierror, nout
    use ScriptInput, only: Script_ReadParameters
    use BoxData, only: BoxArray
    use TrajData, only: TrajArray
!    use AnalysisData, only: AnalysisArray
    use MCMoveData, only: Moves
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: sgrnd
    use Output_DumpCoords, only: Output_DumpData
    implicit none
 
    logical :: accept
    integer :: i, j
    integer :: iMoves, iAtom
    integer(kind=8) :: iCycle, iMove
    real(dp) :: E_T, E_Final
    real(dp) :: avgE, cnt
    character(len=50) :: fileName

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  

    call Script_ReadParameters

    do i = 1, size(BoxArray)
      call BoxArray(i) % box % ComputeEnergy
      call BoxArray(i) % box % BuildNeighList
    enddo

    write(nout, *) "============================================"
    write(nout, *) "       Simulation Start!"
    write(nout, *) "============================================"

    avgE = 0E0_dp
    cnt = 0E0_dp
    !-------Main Monte Carlo Simulation Loop-------
    do iCycle = 1, nCycles
      do iMoves = 1, nMoves
        call Moves(1) % Move % FullMove(BoxArray(1)%box, accept)
        avgE = avgE + BoxArray(1)%box%ETotal
        cnt = cnt + 1E0_dp
      enddo
      if(mod(iCycle, 100) == 0) then
        write(*,*) iCycle, BoxArray(1)%box%ETotal, Moves(1)%Move%GetAcceptRate()
      endif
      if(mod(iCycle, 1000) == 0) then
        call BoxArray(1) % box % NeighList(1) % BuildList
      endif
      if(mod(iCycle, 10) == 0) then
        call Moves(1) % Move % Maintenance
      endif

      if( allocated(TrajArray) ) then
        do i = 1, size(TrajArray)
          if(mod(iCycle, TrajArray(i)%traj%outfreq) == 0) then
            call TrajArray(i) % traj % WriteFrame
          endif
        enddo
      endif
    enddo
    !-------End of Main Monte Carlo Simulation Loop-------
    
    E_Final = BoxArray(1)%box%ETotal
    do i = 1, size(BoxArray)
      call BoxArray(i) % box % ComputeEnergy
    enddo

    write(nout, *) "Culmative Energy:", E_Final
    write(nout, *) "Final Energy:",  BoxArray(1)%box%ETotal
    write(nout, *) "Difference:",  E_Final - BoxArray(1)%box%ETotal
    write(nout, *) "Average Energy:", avgE/cnt

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call Output_DumpData
    write(nout, *) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
