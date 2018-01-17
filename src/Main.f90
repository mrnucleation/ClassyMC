!===========================================================================
  program Classy
    use MPI
    use ParallelVar, only: myid, p_size, ierror, nout
    use SimMonteCarlo, only: RunMonteCarlo
    use ScriptInput, only: Script_ReadParameters
    use VarPrecision
    implicit none

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  

    call Script_ReadParameters

    call RunMonteCarlo
    write(nout, *) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
