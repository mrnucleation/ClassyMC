!===========================================================================
  program Classy
    use MPI
    use SimControl, only: simType
    use ParallelVar, only: myid, p_size, ierror, nout
    use SimMonteCarlo, only: RunMonteCarlo
    use ScriptInput, only: Script_ReadParameters
    use VarPrecision
    implicit none
    real(dp) :: TimeStart, TimeEnd

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  


    if(myid == 0) then
      nout = 6
    else
      nout = 100 + myId
    endif

    call Script_ReadParameters

    call CPU_TIME(TimeStart)
    call RunMonteCarlo
    call CPU_TIME(TimeEnd)

    write(nout, *) "Total Simulation Time: ", TimeEnd - TimeStart
    write(nout, *) "Finished!"
    close(nout)
      
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   

  end program
!===========================================================================
