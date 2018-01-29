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
    character(len=100) :: format_string, fl_name
    character(len=1500) :: outFormat1
 
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  


    if(myid == 0) then
      nout = 6
    else
      nout = 100 
    endif

    if (myid .lt. 10) then
      format_string = "(A,I1,A)"
    elseif(myid .lt. 100) then
      format_string = "(A,I2,A)"
    elseif(myid .lt. 1000) then
      format_string = "(A,I3,A)"
    elseif(myid .lt. 10000) then
      format_string = "(A,I4,A)"          
    else
      format_string = "(A,I5,A)"      
    endif      

    write(fl_name,format_string) "ScreenOutput_", myid,".txt"      
    open( unit=100, file=trim(adjustl(fl_name)) ) 
 

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
