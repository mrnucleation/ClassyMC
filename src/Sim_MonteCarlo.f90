!===========================================================================
module SimMonteCarlo
  use ParallelVar, only: myid, ierror, nout
!===========================================================================
contains
!===========================================================================
  subroutine RunMonteCarlo
    use VarPrecision
    use MPI

    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use CommonSampling, only: Sampling
    use Common_MolInfo, only: nAtomTypes, nMolTypes, MolData
    use Debug, only: Debug_DumpNeiList
    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use MoveClassDef, only: MCMove
    use MultiBoxMoveDef, only: MCMultiBoxMove
    use Output_DumpCoords, only: Output_DumpData
    use RandomGen, only: sgrnd, grnd, ListRNG
    use SimControl, only: nMoves, nCycles
    use Units, only: outEngUnit

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
    write(nout, *) "============================================"
    write(nout, *) "       Simulation Start!"
    write(nout, *) "============================================"

    flush(nout)
    nBoxes = size(BoxArray)
    boxNum = 1

    call Analyze(iCycle, iMove, accept, .true.)
    call Analyze(iCycle, iMove, accept, .false.)
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
        call Sampling%UpdateStatistics(accept)

        if(accept) then
          call Update(iCycle, iMove, accept)
        endif

        call Analyze(iCycle, iMove, accept, .true.)
      enddo 
      !------End Move Loop
      if(mod(iCycle, 1000) == 0) then
        write(nout, *) iCycle, BoxArray(1)%box%ETotal,  BoxArray(1)%box%NMol,(Moves(j)%Move%GetAcceptRate(), j=1, size(Moves))
        flush(nout)
      endif

!      if(mod(iCycle, 100) == 0) then
!        call BoxArray(1) % box % NeighList(1) % BuildList
!      endif

      call Analyze(iCycle, iMove, accept, .false.)
      call Maintenance(iCycle, iMove)
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
          if((mod(iCycle, AnalysisArray(iAn)%func%UpdateFreq) == 0) .or.  AnalysisArray(iAn)%func%perMove) then
!            write(*,*) iCycle, iMove, AnalysisArray(iAn)%func%perMove, moveloop
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

!    if(printBox) then
      do i = 1, size(BoxArray)
        if(mod(iCycle, BoxArray(i)%box%maintFreq) == 0) then
          call BoxArray(i) % box % Maintenance
        endif
      enddo
!    endif


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
!      if(mod(iCycle, BoxArray(i)%box%maintFreq) == 0) then
        call BoxArray(i) % box % Maintenance
!      endif
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
  subroutine Update(iCycle, iMove, accept)
    use AnalysisData, only: AnalysisArray
    use BoxData, only: BoxArray
    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    use CommonSampling, only: Sampling
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    logical, intent(in) :: accept
    integer :: i

    if( .not. accept) then
      return
    endif

    call Sampling % Update

    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Update
      enddo
    endif

    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Update
      enddo
    endif

    do i = 1, size(EnergyCalculator)
      call EnergyCalculator(i)%method%Update
    enddo


    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Update
    enddo

    do i = 1, size(Moves)
      call Moves(i) % move % Update
    enddo

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
    use CommonSampling, only: Sampling
    use ForcefieldData, only: EnergyCalculator
    use MCMoveData, only: Moves, MoveProb
    use TrajData, only: TrajArray
    implicit none
    integer(kind=8), intent(in) :: iCycle, iMove
    integer :: i


    if( allocated(AnalysisArray) ) then
      do i = 1, size(AnalysisArray)
        call AnalysisArray(i) % func % Prologue
      enddo
    endif

    call Sampling % Prologue
!    write(nout, *) "Traj"
    if( allocated(TrajArray) ) then
      do i = 1, size(TrajArray)
        call TrajArray(i) % traj % Prologue
      enddo
    endif

!    write(nout,*) "ECalc"
    do i = 1, size(EnergyCalculator)
      call EnergyCalculator(i)%method%Prologue
    enddo


!    write(nout,*) "Box"
    do i = 1, size(BoxArray)
      call BoxArray(i) % box % Prologue

    enddo

!    write(nout,*) "Moves"
!    flush(nout)
    do i = 1, size(Moves)
      call Moves(i) % move % Prologue
    enddo

  end subroutine

!===========================================================================
end module
!===========================================================================
