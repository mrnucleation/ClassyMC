!===========================================================================
  program Classy
    use VarPrecision
    use MPI
    use ParallelVar, only: myid, p_size, ierror, nout
    use ScriptInput, only: Script_ReadParameters
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use ForcefieldData, only: EnergyCalculator
    use RandomGen, only: sgrnd
    use RSqListDef, only: RSqList
    use TrajectoryTemplate, only: trajectory
    use Traj_XYZ, only: trajXYZ
    implicit none
 
    integer :: i, j
    integer :: nMoves, iAtom
    real(dp) :: E_T, E_Final
    real(dp) :: avgE, cnt
    character(len=50) :: fileName

    class(trajectory), allocatable :: TrajFile

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  

    call sgrnd(1) 

    call Script_ReadParameters

    allocate( RSqList::BoxArray(1)%box%NeighList(1:1) )
    allocate( trajXYZ::TrajFile )
!    call BoxArray(1)%box%LoadCoordinates("Dummy.xyz")
    call EnergyCalculator(1)%Method%Constructor
    call BoxArray(1)%box%NeighList(1)%constructor(1)

    BoxArray(1)%box%temperature = 0.8E0_dp
    BoxArray(1)%box%beta = 1E0_dp/BoxArray(1)%box%temperature

    fileName = "Traj.xyz"
    call TrajFile%SetUnit(100)
    call TrajFile%SetBox(1)
    call TrajFile%SetFreq(1000)
    call TrajFile%SetFileName(fileName)
    call TrajFile%OpenFile

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
         call TrajFile % WriteFrame
       endif
       if(mod(nMoves, 1000) == 0) then
         call BoxArray(1) % box % BuildNeighList
       endif
       if(mod(nMoves, 100) == 0) then
         call Moves(1) % Move % Maintenance
       endif
    enddo
    call TrajFile % CloseFile
    
    E_Final = BoxArray(1)%box%ETotal

!    call BoxArray(1) % box % EFunc % Method % DetailedECalc( BoxArray(1)%box, BoxArray(1)%box%ETotal )
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
