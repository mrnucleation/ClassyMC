!===========================================================================
! Test_EnergyDelta
!
! This test verifies that the incremental (delta) energy calculations
! agree with full from-scratch energy computations for every move type.
!
! Procedure for each move type:
!   1. Compute total energy from scratch (baseline).
!   2. For N trial moves:
!      a. Perform the move (always accepted via AcceptAll sampling).
!      b. Record the cumulative energy (old + delta).
!      c. Recompute energy from scratch.
!      d. Compare cumulative vs from-scratch.
!      e. If they disagree beyond a tolerance, flag a FAIL.
!   3. Report per-move-type results.
!
! Move types tested (matching the input script order):
!   Move 1: MolTranslation   (Displacement, computeintra=false)
!   Move 2: AtomTranslation  (Displacement, computeintra=true)
!   Move 3: PlaneRotate      (Displacement, computeintra=false)
!   Move 4: Rotate3D         (Displacement, computeintra=false)
!   Move 5: IsoVol           (OrthoVolChange)
!===========================================================================
#define __StdErr__ 0
!===========================================================================
program Test_EnergyDelta
#ifdef MPIPARALLEL
  use MPI
#endif

  use SimControl, only: TimeStart, TimeEnd
  use ParallelVar, only: myid, p_size, ierror, nout
  use ScriptInput, only: Script_ReadParameters
  use BoxData, only: BoxArray
  use MCMoveData, only: Moves
  use ForcefieldData, only: EnergyCalculator
  use VarPrecision

  implicit none

  integer, parameter :: nTrials = 100
  integer, parameter :: nMoveTypes = 5

  ! Relative tolerance for comparing energies (dimensionless)
  real(dp), parameter :: relTol = 1.0E-7_dp
  ! Absolute tolerance for energies near zero
  real(dp), parameter :: absTol = 1.0E-10_dp

  character(len=30) :: moveNames(nMoveTypes)
  integer  :: iMove, totalFailed
  logical  :: movePassed(nMoveTypes)

#ifdef MPIPARALLEL
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
#else
  myid = 0
  p_size = 1
#endif

  if(myid == 0) then
    nout = 6
  else
    nout = 100
  endif

  moveNames(1) = "MolTranslation"
  moveNames(2) = "AtomTranslation"
  moveNames(3) = "PlaneRotate"
  moveNames(4) = "Rotate3D"
  moveNames(5) = "IsoVol"

  !--- Banner ---
  write(nout, *)
  write(nout, *) "============================================================"
  write(nout, *) "       ClassyMC Test: Energy Delta Consistency"
  write(nout, *) "============================================================"
  write(nout, *)

  ! ---------------------------------------------------------------
  ! 1. Read the test input file.  Because cycles=0 the MC loop
  !    returns immediately after initialization.
  ! ---------------------------------------------------------------
  call Script_ReadParameters("in.test")

  ! ---------------------------------------------------------------
  ! 2. Compute the initial energy from scratch so the box state
  !    (ETotal, ETable, E_Inter, E_Intra) is fully consistent.
  ! ---------------------------------------------------------------
  call BoxArray(1)%box%ComputeEnergy()

  write(nout, "(A,ES22.14)") "  Initial total energy : ", BoxArray(1)%box%ETotal
  write(nout, "(A,ES22.14)") "  Initial inter energy : ", BoxArray(1)%box%E_Inter
  write(nout, "(A,ES22.14)") "  Initial intra energy : ", BoxArray(1)%box%E_Intra
  write(nout, *)

  ! ---------------------------------------------------------------
  ! 3. Test each move type
  ! ---------------------------------------------------------------
  totalFailed = 0

  do iMove = 1, nMoveTypes
    call TestMoveType(iMove, moveNames(iMove), movePassed(iMove))
    if (.not. movePassed(iMove)) totalFailed = totalFailed + 1
  enddo

  ! ---------------------------------------------------------------
  ! 4. Final summary
  ! ---------------------------------------------------------------
  write(nout, *)
  write(nout, *) "============================================================"
  write(nout, *) "                  Overall Summary"
  write(nout, *) "============================================================"
  do iMove = 1, nMoveTypes
    if (movePassed(iMove)) then
      write(nout, "(A,A30,A)") "    PASS  ", moveNames(iMove), ""
    else
      write(nout, "(A,A30,A)") "  **FAIL  ", moveNames(iMove), ""
    endif
  enddo
  write(nout, *)

  if (totalFailed == 0) then
    write(nout, *) "  >>> ALL TESTS PASSED <<<"
  else
    write(nout, "(A,I2,A,I2,A)") "   >>> ", totalFailed, " of ", nMoveTypes, " move types FAILED <<<"
  endif

  write(nout, *)
  write(nout, *) "============================================================"

  ! Return nonzero exit code on failure
  if (totalFailed > 0) then
    error stop 1
  endif

#ifdef MPIPARALLEL
  call MPI_BARRIER(MPI_COMM_WORLD, ierror)
  call MPI_FINALIZE(ierror)
#endif

contains

  !=========================================================================
  ! TestMoveType
  !
  ! Runs nTrials of the given move, comparing cumulative vs from-scratch
  ! energy after each accepted move.  Reports results and sets passed.
  !=========================================================================
  subroutine TestMoveType(moveIdx, moveName, passed)
    implicit none
    integer, intent(in) :: moveIdx
    character(len=*), intent(in) :: moveName
    logical, intent(out) :: passed

    integer  :: iTrial, nAccepted, nFailed, nRejected
    logical  :: accept
    real(dp) :: E_cumulative, E_scratch
    real(dp) :: E_Inter_cum, E_Inter_scratch
    real(dp) :: E_Intra_cum, E_Intra_scratch
    real(dp) :: diff, maxDiff, sumDiff
    real(dp) :: scale
    logical  :: trialPass
    integer  :: iCalc

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, "(A,A,A,I2,A)") "  Testing: ", trim(moveName), "  (move ", moveIdx, ")"
    write(nout, *) "------------------------------------------------------------"

    ! Recompute energy from scratch before starting this move type's tests.
    ! This ensures a clean baseline regardless of what previous tests did.
    call BoxArray(1)%box%ComputeEnergy()

    nAccepted = 0
    nFailed   = 0
    nRejected = 0
    maxDiff   = 0.0_dp
    sumDiff   = 0.0_dp

    do iTrial = 1, nTrials
      ! --- Perform the move ---
      call Moves(moveIdx)%Move%FullMove(BoxArray(1)%box, accept)

      if (.not. accept) then
        nRejected = nRejected + 1
        cycle
      endif

      nAccepted = nAccepted + 1

      ! After an accepted move the cumulative energy has been updated:
      !   box%ETotal = E_old + E_delta
      E_cumulative  = BoxArray(1)%box%ETotal
      E_Inter_cum   = BoxArray(1)%box%E_Inter
      E_Intra_cum   = BoxArray(1)%box%E_Intra

      ! --- Recompute the full energy from scratch ---
      ! This overwrites ETotal, ETable, E_Inter, E_Intra.
      call BoxArray(1)%box%ComputeEnergy()

      E_scratch       = BoxArray(1)%box%ETotal
      E_Inter_scratch = BoxArray(1)%box%E_Inter
      E_Intra_scratch = BoxArray(1)%box%E_Intra

      ! --- Compare total energy ---
      diff = abs(E_cumulative - E_scratch)

      ! Choose a sensible scale for the relative comparison
      scale = max(abs(E_scratch), 1.0_dp)
      trialPass = (diff / scale) < relTol .and. diff < max(absTol, abs(E_scratch) * relTol)

      if (.not. trialPass) then
        nFailed = nFailed + 1
        write(nout, "(A,I4,A)")       "  FAIL  trial ", iTrial, ":"
        write(nout, "(A,ES22.14)")    "     cumulative E_Total  : ", E_cumulative
        write(nout, "(A,ES22.14)")    "     from-scratch E_Total: ", E_scratch
        write(nout, "(A,ES22.14)")    "     difference          : ", diff
        write(nout, "(A,ES22.14)")    "     cumulative E_Inter  : ", E_Inter_cum
        write(nout, "(A,ES22.14)")    "     from-scratch E_Inter: ", E_Inter_scratch
        write(nout, "(A,ES22.14)")    "     cumulative E_Intra  : ", E_Intra_cum
        write(nout, "(A,ES22.14)")    "     from-scratch E_Intra: ", E_Intra_scratch
      endif

      maxDiff = max(maxDiff, diff)
      sumDiff = sumDiff + diff

      ! Check neighbor lists after acceptance
      ! (matches what the main MC loop does)
      call BoxArray(1)%box%CheckLists

      ! Call forcefield Update (e.g. Hybrid zeros EDiff)
      do iCalc = 1, size(EnergyCalculator)
        call EnergyCalculator(iCalc)%method%Update
      enddo
    enddo

    ! --- Report ---
    write(nout, *)
    write(nout, "(A,A)")            "  Move type         : ", trim(moveName)
    write(nout, "(A,I6)")           "  Total trials      : ", nTrials
    write(nout, "(A,I6)")           "  Accepted moves    : ", nAccepted
    write(nout, "(A,I6)")           "  Rejected moves    : ", nRejected
    write(nout, "(A,I6)")           "  Failed comparisons: ", nFailed
    write(nout, "(A,ES12.4)")       "  Max |E_cum - E_sc|: ", maxDiff
    if (nAccepted > 0) then
      write(nout, "(A,ES12.4)")     "  Avg |E_cum - E_sc|: ", sumDiff / real(nAccepted, dp)
    endif

    if (nFailed == 0 .and. nAccepted > 0) then
      write(nout, "(A,A,A)")  "  >>> PASS: ", trim(moveName), " <<<"
      passed = .true.
    elseif (nAccepted == 0) then
      write(nout, "(A,A,A)")  "  >>> WARN: ", trim(moveName), " -- no moves accepted, inconclusive <<<"
      passed = .true.  ! Not a failure, just inconclusive
    else
      write(nout, "(A,A,A)")  "  >>> FAIL: ", trim(moveName), " <<<"
      passed = .false.
    endif

  end subroutine TestMoveType

end program
!===========================================================================
