!===========================================================================
  program Classy
  !This is a test
#ifdef MPIPARALLEL
    use MPI
#endif

    use SimControl, only: simType, TimeStart, TimeEnd
    use ParallelVar, only: myid, p_size, ierror, nout
    use SimMonteCarlo, only: RunMonteCarlo
    use ScriptInput, only: Script_ReadParameters
#ifdef EMBPYTHON
    use forpy_mod, only: forpy_initialize, list, get_sys_path, forpy_finalize
#endif

    use VarPrecision
    implicit none

#ifdef EMBPYTHON
    type(list) :: paths
#endif

    character(len=100) :: format_string, fl_name
    character(len=1500) :: outFormat1
 
#ifdef MPIPARALLEL
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)  
#else
    myid = 0
    p_size = 1
#endif

#ifdef EMBPYTHON
    ierror = forpy_initialize()
    ierror = get_sys_path(paths)
    ierror = paths%append(".")
#endif

    if(myid == 0) then
      nout = 6
    else
      nout = 100 
    endif

    if (myid < 10) then
      format_string = "(A,I1,A)"
    elseif(myid < 100) then
      format_string = "(A,I2,A)"
    elseif(myid < 1000) then
      format_string = "(A,I3,A)"
    elseif(myid < 10000) then
      format_string = "(A,I4,A)"          
    elseif(myid < 100000) then
      format_string = "(A,I5,A)"          
    else
      format_string = "(A,I6,A)"      
    endif      

    write(nout,*) "============================================================"
    write(nout,*) "             ****  *    *****    ***    ***  * * "
    write(nout,*) "             *     *    * * *    *      *     *  "
    write(nout,*) "             ****  **** *   *  ***    ***     *  "
    write(nout,*) "============================================================"

    write(fl_name,format_string) "ScreenOutput_", myid,".txt"      
    open( unit=100, file=trim(adjustl(fl_name)) ) 
 

    call Script_ReadParameters

    write(nout, *) "Total Simulation Time: ", TimeEnd - TimeStart
    write(nout, *) "Finished!"
    close(nout)
      
#ifdef MPIPARALLEL
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       
    call MPI_FINALIZE(ierror)   
#endif

#ifdef EMBPYTHON
  call forpy_finalize
#endif

  end program
!===========================================================================
