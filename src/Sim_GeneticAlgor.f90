!===========================================================================
module SimGeneticAlgorithm
  integer :: poolSize 
  integer :: maxSize
!===========================================================================
contains
!===========================================================================
  subroutine RunGeneticAlgorithm
    use VarPrecision
    use MPI

    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use MCMoveData, only: Moves, MoveProb
    use Output_DumpCoords, only: Output_DumpData
    use ParallelVar, only: myid, ierror, nout
    use RandomGen, only: sgrnd, ListRNG
    use SimControl, only: nMoves, nCycles

    use MultiBoxMoveDef, only: MCMultiBoxMove
    implicit none
 
    logical :: accept
    integer :: iAtom, moveNum
    integer(kind=8) :: iCycle, iMove


    write(nout, *) "============================================"
    write(nout, *) "       Simulation Start!"
    write(nout, *) "============================================"

    !-------Main GA Simulation Loop-------
    do iCycle = 1, nCycles
      moveNum = ListRNG(MoveProb)
!      select type( move => Moves(moveNum) % Move )
!        type is (MCMultiBoxMove)
!          call Moves(moveNum) % Move % MultiBox(accept)
!      end select
    enddo
    !-------End of Main GA Simulation Loop-------
    

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       

    call Output_DumpData
      
  end subroutine
!===========================================================================
end module
!===========================================================================
