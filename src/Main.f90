!===========================================================================
  program Classy
    use VarPrecision
    use MPI
    use ParallelVar, only: myid, p_size, ierror, nout
    use ScriptInput, only: Script_ReadParameters
    use BoxData, only: BoxArray, Constrain
    use SimpleSimBox, only: SimpleBox
    use CubicBoxDef, only: CubeBox
    use Common_MolDef, only: nAtomTypes
    use ForcefieldData, only: EnergyCalculator
    use AtomTranslation, only: AtomMolTranslate
    use RandomGen, only: sgrnd
    use DistanceCriteria, only: distcriteria
    use CommonSampling, only: sampling, metropolis
    implicit none
 
    integer :: nMoves, iAtom
    real(dp) :: E_T
    real(dp) :: avgE, cnt
    class(AtomMolTranslate), allocatable :: MCMover

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  

    call sgrnd(1) 

    call Script_ReadParameters

!    allocate( Constrain(1:1) )
!    allocate( distcriteria::Constrain(1)%Method )
    allocate( BoxArray(1)%box%NeighList(1:1) )
    call BoxArray(1)%box%LoadCoordinates("Dummy.xyz")
    allocate( metropolis::sampling )

    allocate( AtomMolTranslate::MCMover )
    call EnergyCalculator(1)%Method%Constructor

    BoxArray(1)%box%AtomType = 1
    BoxArray(1)%box%temperature = 0.7E0_dp
    BoxArray(1)%box%beta = 1E0_dp/BoxArray(1)%box%temperature
!    call BoxArray(1)%box%DummyCoords
    call BoxArray(1)%box % EFunc % Method % DetailedECalc( BoxArray(1)%box, BoxArray(1)%box%ETotal )

    open(unit=2, file="Traj.xyz")
    write(2,*) BoxArray(1)%box%nAtoms
    write(2,*) 
    do iAtom = 1, BoxArray(1)%box%nAtoms
      write(2,*) "Ar",BoxArray(1)%box%atoms(1, iAtom), BoxArray(1)%box%atoms(2, iAtom), BoxArray(1)%box%atoms(3, iAtom)
    enddo


    avgE = 0E0_dp
    cnt = 0E0_dp
    do nMoves = 1, nint(1d7)
       call MCMover % FullMove(BoxArray(1)%box)
       avgE = avgE + BoxArray(1)%box%ETotal
       cnt = cnt + 1E0_dp
       if(mod(nMoves, 1000) == 0) then
         write(*,*) nMoves, BoxArray(1)%box%ETotal
         write(2,*) BoxArray(1)%box%nAtoms
         write(2,*) 
         do iAtom = 1, BoxArray(1)%box%nAtoms
           write(2,*) "Ar",BoxArray(1)%box%atoms(1, iAtom), BoxArray(1)%box%atoms(2, iAtom), BoxArray(1)%box%atoms(3, iAtom)
         enddo
       endif
    enddo
    call BoxArray(1)%box % EFunc % Method % DetailedECalc( BoxArray(1)%box, BoxArray(1)%box%ETotal )

    write(*,*) "Average Energy:", avgE/cnt

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    write(nout,*) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
