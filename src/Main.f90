!===========================================================================
  program Classy
    use VarPrecision
    use MPI
    use ParallelVar, only: myid, p_size, ierror, nout
    use ScriptInput, only: Script_ReadParameters
    use BoxData, only: BoxArray
    use Common_MolDef, only: nAtomTypes
    use ForcefieldData, only: EnergyCalculator
    use AtomTranslation, only: AtomMolTranslate
    use RandomGen, only: sgrnd
    implicit none
 
    
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

    allocate( AtomMolTranslate::MCMover )
!    allocate( BoxArray(1)%AtomType(1:2) )
    call EnergyCalculator(1)%Method%Constructor

    BoxArray(1)%AtomType = 1
    call BoxArray(1)%DummyCoords
    call EnergyCalculator(BoxArray(1)%ECalcer)%Method%DetailedECalc( BoxArray(1), E_T)

    call MCMover % FullMove(BoxArray(1))
    call EnergyCalculator(BoxArray(1)%ECalcer)%Method%DetailedECalc( BoxArray(1), E_T)

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    write(nout,*) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
