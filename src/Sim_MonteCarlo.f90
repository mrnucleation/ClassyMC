!===========================================================================
module SimMonteCarlo
!===========================================================================
contains
!===========================================================================
  subroutine RunMonteCarlo
    use VarPrecision
    use MPI

    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use Debug, only: Debug_DumpNeiList
    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use MoveClassDef, only: MCMove
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use Output_DumpCoords, only: Output_DumpData
    use ParallelVar, only: myid, ierror, nout
    use RandomGen, only: sgrnd, grnd, ListRNG
    use SimControl, only: nMoves, nCycles

    implicit none
 
    logical :: accept
    integer :: i, j, nBoxes
    integer :: iAtom, moveNum, boxNum
    integer(kind=8) :: iCycle, iMove
    character(len=50) :: fileName
    class(MCMove), pointer :: curMove

    iCycle = 0
    iMove = 0

    call Prologue(iCycle, iMove)

    call Trajectory(iCycle, iMove)
    write(nout, *) "============================================"
    write(nout, *) "       Simulation Start!"
    write(nout, *) "============================================"

    nBoxes = size(BoxArray)
    boxNum = 1
    !-------Main Monte Carlo Simulation Loop-------
    do iCycle = 1, nCycles

      !-----Start Move Loop
      do iMove = 1, nMoves
        moveNum = ListRNG(MoveProb)
        curMove => Moves(moveNum) % Move 
        select type( curMove )
          class is (MCMultiBoxMove)
            call curMove % MultiBox (accept)
 
          class is (MCMove)
            if(nBoxes > 1) then
              boxNum = floor(grnd()*nBoxes + 1E0_dp)
            endif
            call curMove % FullMove(BoxArray(boxNum)%box, accept)

        end select
        call Analyze(iCycle, iMove, accept, .true.)
      enddo 
      !------End Move Loop
      if(mod(iCycle, 1000) == 0) then
        write(*,*) iCycle, BoxArray(1)%box%ETotal, (Moves(j)%Move%GetAcceptRate(), j=1, size(Moves))
      endif

      if(mod(iCycle, 100) == 0) then
        call BoxArray(1) % box % NeighList(1) % BuildList
      endif

      call Analyze(iCycle, iMove, accept, .false.)
      call Trajectory(iCycle, iMove)
    enddo
    !-------End of Main Monte Carlo Simulation Loop-------
 
    write(nout,*) "======================================="
    write(nout,*) "     Simulation End"
    write(nout,*) "======================================="

    call Epilogue(iCycle, iMove)

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)       

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Finalize
      enddo
    endif

    call Output_DumpData
      
  end subroutine
!===========================================================================
  subroutine Analyze(iCycle, iMove, accept, moveloop)
    use AnalysisData, only: AnalysisArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    logical, intent(in) :: accept
    logical, intent(in) :: moveloop
    integer :: iAn

    if( allocated(AnalysisArray) ) then
      do iAn = 1, size(AnalysisArray)
        if( AnalysisArray(iAn)%func%perMove .eqv. moveloop) then
          if(mod(iCycle, AnalysisArray(iAn)%func%UpdateFreq) == 0) then
            call AnalysisArray(iAn) % func % Compute(accept)
          endif
        endif
      enddo
    endif
 
  end subroutine
!===========================================================================
  subroutine ScreenOut(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    use SimControl, only: printBox, printAcc
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i

    if(printBox) then
      do i = 1, size(BoxArray)
        if(mod(iCycle, BoxArray(i)%box%maintFreq) == 0) then
          call BoxArray(i) % box % Maintenance
        endif
      enddo
    endif


    if(printAcc) then
      do i = 1, size(Moves)
        if(mod(iCycle, Moves(i)%move%maintFreq) == 0) then
          call Moves(i) % move % Maintenance
        endif
      enddo
    endif

  end subroutine
!===========================================================================
  subroutine Maintenance(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i

    if(mod(iCycle, Sampling%maintFreq) == 0 ) then
      call Sampling % Maintenance
    endif

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        if(mod(iCycle, AnalysisArray(i)%func%maintFreq) == 0) then
          call AnalysisArray(i) % func % Maintenance
        endif
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        if(mod(iCycle, TrajArray(i)%traj%maintFreq) == 0) then
          call TrajArray(i) % traj % Maintenance
        endif
      enddo
    endif

    do i = 1, size(BoxArray)
      if(mod(iCycle, BoxArray(i)%box%maintFreq) == 0) then
        call BoxArray(i) % box % Maintenance
      endif
    enddo


    do i = 1, size(Moves)
      if(mod(iCycle, Moves(i)%move%maintFreq) == 0) then
        call Moves(i) % move % Maintenance
      endif
    enddo

  end subroutine
!===========================================================================
  subroutine Trajectory(iCycle, iMove)
    use TrajData, only: TrajArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: iTraj

    if( allocated(TrajArray) ) then
      do iTraj = 1, size(TrajArray)
        if(mod(iCycle, TrajArray(iTraj)%traj%outfreq) == 0) then
          call TrajArray(iTraj) % traj % WriteFrame
        endif
      enddo
    endif

  end subroutine
!===========================================================================
  subroutine Epilogue(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    use ParallelVar, only: nout
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i

    call Sampling % Epilogue
    
    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Epilogue
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Epilogue
      enddo
    endif

    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Epilogue
    enddo
    write(nout, *) "-----------------------"

    do i = 1, size(Moves)
      call Moves(i) % move % Epilogue
    enddo

  end subroutine
!===========================================================================
  subroutine Prologue(iCycle, iMove)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i

    call Sampling % Prologue

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Prologue
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Prologue
      enddo
    endif

    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Prologue
    enddo

    do i = 1, size(Moves)
      call Moves(i) % move % Prologue
    enddo

  end subroutine

!===========================================================================
end module
!===========================================================================
