!===========================================================================
  program Classy
    use VarPrecision
    use MPI
    use ParallelVar, only: myid, p_size, ierror, nout
    use ScriptInput, only: Script_ReadParameters
    use BoxData, only: BoxArray
    use TrajData, only: TrajArray
    use MCMoveData, only: Moves
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: sgrnd
    use RSqListDef, only: RSqList
    implicit none
 
    integer :: i, j
    integer :: nMoves, iAtom
    real(dp) :: E_T, E_Final
    real(dp) :: avgE, cnt
    character(len=50) :: fileName

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  


    call Script_ReadParameters
    call sgrnd(1) 

!    allocate( RSqList::BoxArray(1)%box%NeighList(1:1) )
    call EnergyCalculator(1)%Method%Constructor
    call BoxArray(1)%box%NeighList(1)%constructor(1)

    call BoxArray(1) % box % ComputeEnergy
    call BoxArray(1) % box % BuildNeighList
    write(nout,*) "Simulation Start!"
    avgE = 0E0_dp
    cnt = 0E0_dp
    do nMoves = 1, nint(1d5)
       call Moves(1) % Move % FullMove(BoxArray(1)%box)
       avgE = avgE + BoxArray(1)%box%ETotal
       cnt = cnt + 1E0_dp
       if(mod(nMoves, 1000) == 0) then
         write(*,*) nMoves, BoxArray(1)%box%ETotal, Moves(1)%Move%GetAcceptRate()
       endif
       if(mod(nMoves, 1000) == 0) then
         call BoxArray(1) % box % BuildNeighList
       endif
       if(mod(nMoves, 100) == 0) then
         call Moves(1) % Move % Maintenance
       endif

       if( allocated(TrajArray) ) then
         do i = 1, size(TrajArray)
           if(mod(nMoves, TrajArray(i)%traj%outfreq) == 0) then
             call TrajArray(i) % traj % WriteFrame
           endif
         enddo
       endif

    enddo
    
    E_Final = BoxArray(1)%box%ETotal

    call BoxArray(1) % box % ComputeEnergy
    write(nout, *) "Culmative Energy:", E_Final
    write(nout, *) "Final Energy:",  BoxArray(1)%box%ETotal
    write(nout, *) "Difference:",  E_Final - BoxArray(1)%box%ETotal
    write(nout, *) "Average Energy:", avgE/cnt

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    write(nout,*) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
