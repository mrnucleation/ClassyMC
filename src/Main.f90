!===========================================================================
  program Classy
    use VarPrecision
    use MPI
    use ParallelVar, only: myid, p_size, ierror, nout
    use ScriptInput, only: Script_ReadParameters
    use BoxData, only: BoxArray, Constrain
    use Common_MolDef, only: nAtomTypes
    use ForcefieldData, only: EnergyCalculator
    use AtomTranslation, only: AtomMolTranslate
    use RandomGen, only: sgrnd
    use DistanceCriteria, only: distcriteria
    implicit none
 
    integer :: nMoves
    real(dp) :: E_T
    class(AtomMolTranslate), allocatable :: MCMover

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  

    call sgrnd(1) 

    call Script_ReadParameters
    allocate( BoxArray(1)%atoms(1:3, 1:2) )
    allocate( BoxArray(1)%AtomType(1:2) )
    allocate( BoxArray(1)%ETable(1:2) )
    allocate( BoxArray(1)%dETable(1:2) )
    allocate( Constrain(1:1) )
    allocate( distcriteria::Constrain(1)%Method )
    allocate( BoxArray(1)%NeighList(1:1) )

    allocate( AtomMolTranslate::MCMover )
!    allocate( BoxArray(1)%AtomType(1:2) )
    call EnergyCalculator(1)%Method%Constructor

    BoxArray(1)%AtomType = 1
    call BoxArray(1)%DummyCoords
    call EnergyCalculator(BoxArray(1)%ECalcer) % Method % DetailedECalc( BoxArray(1), BoxArray(1)%ETotal )

    open(unit=2, file="Traj.xyz")
         write(2,*) 2
         write(2,*) 
         write(2,*) "Ar",BoxArray(1)%atoms(1, 1), BoxArray(1)%atoms(2, 1), BoxArray(1)%atoms(3, 1)
         write(2,*) "Ar",BoxArray(1)%atoms(1, 2), BoxArray(1)%atoms(2, 2), BoxArray(1)%atoms(3, 2)
    do nMoves = 1, nint(1d6)
       call MCMover % FullMove(BoxArray(1))
       if(mod(nMoves, 1000) == 0) then
         write(*,*) nMoves
         write(2,*) 2
         write(2,*) 
         write(2,*) "Ar",BoxArray(1)%atoms(1, 1), BoxArray(1)%atoms(2, 1), BoxArray(1)%atoms(3, 1)
         write(2,*) "Ar",BoxArray(1)%atoms(1, 2), BoxArray(1)%atoms(2, 2), BoxArray(1)%atoms(3, 2)
       endif
    enddo
    call EnergyCalculator(BoxArray(1)%ECalcer) % Method % DetailedECalc( BoxArray(1), BoxArray(1)%ETotal )

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    write(nout,*) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
