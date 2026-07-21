!===========================================================================
! Test_ClusterCriteria
!
! Verifies the cluster-connectivity constraints (DistCriteria and
! DistCriteriaAllMol) by constructing deterministic atom geometries and
! checking that Deletion perturbations are correctly accepted or rejected.
!
! Three geometries are tested for each constraint type:
!
!   1. Straight chain  -- atoms linearly spaced just under rCut.
!      Deleting an interior atom breaks the chain (reject); deleting an
!      end atom leaves a shorter connected chain (accept).
!
!   2. Bridged cluster -- two dense sub-clusters connected by a single
!      bridge atom.  Deleting the bridge disconnects them (reject);
!      deleting any other atom leaves the cluster connected (accept).
!
!   3. Fully connected -- every pair within rCut, so removing any single
!      atom still leaves a connected cluster (accept for all).
!
! The test uses 1-atom molecules so each molecule maps directly to a
! single position, isolating the topology / BFS logic from multi-atom
! bookkeeping.
!===========================================================================
program test_cluster_criteria
  use VarPrecision
  use CoordinateTypes, only: Deletion, Addition, Displacement
  use BoxData, only: BoxArray
  use SimpleSimBox, only: SimpleBox
  use OrthoBoxDef, only: OrthoBox
  use Constrain_DistanceCriteria, only: DistCriteria
  use Constrain_DistanceCriteriaAllMol, only: DistCriteriaAllMol
  use ParallelVar, only: nout

  implicit none

  integer, parameter :: nTests = 11
  logical :: results(nTests)
  integer :: totalFailed, i
  character(len=40) :: testNames(nTests)

  nout = 6

  testNames(1) = "DistCriteria: Chain"
  testNames(2) = "DistCriteria: Bridge"
  testNames(3) = "DistCriteria: Fully Connected"
  testNames(4) = "DistCriteria: Broken + Addition"
  testNames(5) = "AllMolDistCriteria: Chain"
  testNames(6) = "AllMolDistCriteria: Bridge"
  testNames(7) = "AllMolDistCriteria: Fully Connected"
  testNames(8) = "AllMolDistCriteria: Broken + Addition"
  testNames(9) = "AllMolDistCriteria: Mixed rCut Chain"
  testNames(10) = "AllMolDistCriteria: Alternating Chain"
  testNames(11) = "DistCriteria: MC Loop Safety Check"

  write(nout, *)
  write(nout, *) "============================================================"
  write(nout, *) "       ClassyMC Test: Cluster Criteria Constraints"
  write(nout, *) "============================================================"

  call test_distcrit_chain(results(1))
  call test_distcrit_bridge(results(2))
  call test_distcrit_dense(results(3))
  call test_distcrit_addition(results(4))
  call test_allmol_chain(results(5))
  call test_allmol_bridge(results(6))
  call test_allmol_dense(results(7))
  call test_allmol_addition(results(8))
  call test_allmol_mixed_rcut_chain(results(9))
  call test_allmol_alternating_chain(results(10))
  call test_distcrit_mc_loop(results(11))

  write(nout, *)
  write(nout, *) "============================================================"
  write(nout, *) "                  Overall Summary"
  write(nout, *) "============================================================"
  totalFailed = 0
  do i = 1, nTests
    if(results(i)) then
      write(nout, "(A,A)") "   PASS  ", trim(testNames(i))
    else
      write(nout, "(A,A)") " **FAIL  ", trim(testNames(i))
      totalFailed = totalFailed + 1
    endif
  enddo
  write(nout, *)

  if(totalFailed == 0) then
    write(nout, *) "  >>> ALL TESTS PASSED <<<"
  else
    write(nout, "(A,I2,A)") "   >>> ", totalFailed, " test(s) FAILED <<<"
    error stop 1
  endif

contains

!===========================================================================
! SetupSingleTypeSystem
!
! Creates a minimal box with one molecule type (1 atom per molecule).
! After return the caller must fill box%atoms with desired positions.
!===========================================================================
  subroutine SetupSingleTypeSystem(NMolMaxIn, NMol0In, boxL)
    use Common_MolInfo, only: nMolTypes, nAtomTypes, MolData, AtomData
    implicit none
    integer, intent(in) :: NMolMaxIn, NMol0In
    real(dp), intent(in) :: boxL

    class(SimpleBox), pointer :: box => null()
    character(len=80) :: dimline
    integer :: lineStat

    nMolTypes = 1
    nAtomTypes = 1
    if(allocated(MolData)) deallocate(MolData)
    if(allocated(AtomData)) deallocate(AtomData)
    allocate(MolData(1:1))
    allocate(AtomData(1:1))
    AtomData(1)%Symb = "A"
    AtomData(1)%mass = 1.0E0_dp
    MolData(1)%nAtoms = 1
    allocate(MolData(1)%atomType(1:1))
    MolData(1)%atomType(1) = 1

    if(allocated(BoxArray)) deallocate(BoxArray)
    allocate(BoxArray(1:1))
    allocate(OrthoBox :: BoxArray(1)%box)
    box => BoxArray(1)%box
    box%boxID = 1

    call box%AllocateMolBound()
    box%NMolMin(1) = 0
    box%NMol(1) = NMol0In
    box%NMolMax(1) = NMolMaxIn
    call box%Constructor()
    box%nMolTotal = NMol0In

    write(dimline, "(A,3(F12.5,1x))") "dimension ", boxL, boxL, boxL
    call box%LoadDimension(dimline, lineStat)

  end subroutine

!===========================================================================
! SetupTwoTypeSystem
!
! Creates a minimal box with two molecule types (1 atom per molecule each).
! nMol1 molecules of type 1 and nMol2 molecules of type 2.
!===========================================================================
  subroutine SetupTwoTypeSystem(NMolMax1, NMol1, NMolMax2, NMol2, boxL)
    use Common_MolInfo, only: nMolTypes, nAtomTypes, MolData, AtomData
    implicit none
    integer, intent(in) :: NMolMax1, NMol1, NMolMax2, NMol2
    real(dp), intent(in) :: boxL

    class(SimpleBox), pointer :: box => null()
    character(len=80) :: dimline
    integer :: lineStat

    nMolTypes = 2
    nAtomTypes = 2
    if(allocated(MolData)) deallocate(MolData)
    if(allocated(AtomData)) deallocate(AtomData)
    allocate(MolData(1:2))
    allocate(AtomData(1:2))
    AtomData(1)%Symb = "A"
    AtomData(1)%mass = 1.0E0_dp
    AtomData(2)%Symb = "B"
    AtomData(2)%mass = 1.0E0_dp
    MolData(1)%nAtoms = 1
    allocate(MolData(1)%atomType(1:1))
    MolData(1)%atomType(1) = 1
    MolData(2)%nAtoms = 1
    allocate(MolData(2)%atomType(1:1))
    MolData(2)%atomType(1) = 2

    if(allocated(BoxArray)) deallocate(BoxArray)
    allocate(BoxArray(1:1))
    allocate(OrthoBox :: BoxArray(1)%box)
    box => BoxArray(1)%box
    box%boxID = 1

    call box%AllocateMolBound()
    box%NMolMin(1) = 0
    box%NMol(1) = NMol1
    box%NMolMax(1) = NMolMax1
    box%NMolMin(2) = 0
    box%NMol(2) = NMol2
    box%NMolMax(2) = NMolMax2
    call box%Constructor()
    box%nMolTotal = NMol1 + NMol2

    write(dimline, "(A,3(F12.5,1x))") "dimension ", boxL, boxL, boxL
    call box%LoadDimension(dimline, lineStat)

  end subroutine

!===========================================================================
! CheckDeletion
!
! Calls DiffCheck on the constraint with a Deletion perturbation for
! molecule subIndx of type molType.  Returns the accept/reject result.
!===========================================================================
  subroutine CheckDeletion(crit, molType, subIndx, accept)
    use ConstraintTemplate, only: constraint
    implicit none
    class(constraint), intent(inout) :: crit
    integer, intent(in) :: molType, subIndx
    logical, intent(out) :: accept

    class(SimpleBox), pointer :: box => null()
    type(Deletion) :: delDisp(1)
    integer :: molIndx

    box => BoxArray(1)%box
    molIndx = box%MolGlobalIndx(molType, subIndx)
    delDisp(1)%molType = molType
    delDisp(1)%molIndx = molIndx
    delDisp(1)%atmIndx = box%MolStartIndx(molIndx)
    call crit%DiffCheck(box, delDisp, accept)

  end subroutine

!===========================================================================
! CheckDeletionByGlobal
!
! Like CheckDeletion but takes the global molecule index directly (for
! AllMol tests where we know the global index).
!===========================================================================
  subroutine CheckDeletionByGlobal(crit, globalIndx, accept)
    use ConstraintTemplate, only: constraint
    implicit none
    class(constraint), intent(inout) :: crit
    integer, intent(in) :: globalIndx
    logical, intent(out) :: accept

    class(SimpleBox), pointer :: box => null()
    type(Deletion) :: delDisp(1)
    integer :: mType

    box => BoxArray(1)%box
    mType = box%MolType(box%MolStartIndx(globalIndx))
    delDisp(1)%molType = mType
    delDisp(1)%molIndx = globalIndx
    delDisp(1)%atmIndx = box%MolStartIndx(globalIndx)
    call crit%DiffCheck(box, delDisp, accept)

  end subroutine

!===========================================================================
! CheckAddition
!
! Calls DiffCheck on the constraint with an Addition perturbation for a
! new molecule of type molType at position (x, y, z).  The new molecule
! occupies sub-index NMol(molType)+1, i.e. the first empty slot.
!===========================================================================
  subroutine CheckAddition(crit, molType, x, y, z, accept)
    use ConstraintTemplate, only: constraint
    implicit none
    class(constraint), intent(inout) :: crit
    integer, intent(in) :: molType
    real(dp), intent(in) :: x, y, z
    logical, intent(out) :: accept

    class(SimpleBox), pointer :: box => null()
    type(Addition) :: addDisp(1)
    integer :: newSubIndx, molIndx

    box => BoxArray(1)%box
    newSubIndx = box%NMol(molType) + 1
    molIndx = box%MolGlobalIndx(molType, newSubIndx)
    addDisp(1)%molType = molType
    addDisp(1)%molIndx = molIndx
    addDisp(1)%atmIndx = box%MolStartIndx(molIndx)
    addDisp(1)%x_new = x
    addDisp(1)%y_new = y
    addDisp(1)%z_new = z
    call crit%DiffCheck(box, addDisp, accept)

  end subroutine

!===========================================================================
!                     DistCriteria Tests
!===========================================================================
! test_distcrit_chain
!
!   Five atoms in a line, each 2.9 apart with rCut = 3.0:
!
!     1 ---- 2 ---- 3 ---- 4 ---- 5
!     0     2.9    5.8    8.7   11.6
!
!   Deleting any interior atom (2, 3, 4) breaks the chain.
!   Deleting either end (1, 5) leaves a connected sub-chain.
!===========================================================================
  subroutine test_distcrit_chain(passed)
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol = 5
    real(dp), parameter :: rCut = 3.0E0_dp
    real(dp), parameter :: spacing = 2.9E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteria) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iMol

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: DistCriteria - Straight Chain"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSingleTypeSystem(nMol, nMol, boxL)
    box => BoxArray(1)%box

    do iMol = 1, nMol
      box%atoms(1, iMol) = (iMol - 1) * spacing
      box%atoms(2, iMol) = 0.0E0_dp
      box%atoms(3, iMol) = 0.0E0_dp
    enddo

    crit%molType = 1
    crit%atomNum = 1
    crit%rCut = rCut
    crit%rCutSq = rCut * rCut
    call crit%Constructor(1)
    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for a valid chain"
      nFail = nFail + 1
    endif

    ! --- End deletions should be accepted ---
    call CheckDeletion(crit, 1, 1, accept)
    if(.not. accept) then
      write(nout, *) "  FAIL: deleting end atom 1 should be accepted"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete end atom 1 -> accepted"
    endif

    call CheckDeletion(crit, 1, nMol, accept)
    if(.not. accept) then
      write(nout, *) "  FAIL: deleting end atom 5 should be accepted"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete end atom 5 -> accepted"
    endif

    ! --- Interior deletions should be rejected ---
    do iMol = 2, nMol - 1
      call CheckDeletion(crit, 1, iMol, accept)
      if(accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting interior atom ", iMol, " should be rejected"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete interior atom ", iMol, " -> rejected"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: DistCriteria Chain <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: DistCriteria Chain <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_distcrit_bridge
!
!   Two dense sub-clusters connected only through a single bridge atom.
!   rCut = 6.0
!
!     Sub-cluster A:               Bridge:        Sub-cluster B:
!       1 (0, 0, 0)                  4 (4.5, 0, 0)  5 (9, 0, 0)
!       2 (1, 0, 0)                                   6 (10, 0, 0)
!       3 (0.5, 0.866, 0)                             7 (9.5, 0.866, 0)
!
!   All A-A and B-B pairs are within rCut (triangles with ~1.0 sides).
!   Bridge (4) is within rCut of all A and all B atoms.
!   No direct A-B connection exists (nearest A-B pair is ~8.0 apart).
!
!   Deleting the bridge disconnects A from B (reject).
!   Deleting any other atom leaves the cluster connected (accept).
!===========================================================================
  subroutine test_distcrit_bridge(passed)
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol = 7
    real(dp), parameter :: rCut = 6.0E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteria) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iMol

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: DistCriteria - Bridged Cluster"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSingleTypeSystem(nMol, nMol, boxL)
    box => BoxArray(1)%box

    ! Sub-cluster A (triangle)
    box%atoms(1:3, 1) = [0.0_dp,   0.0_dp,    0.0_dp]
    box%atoms(1:3, 2) = [1.0_dp,   0.0_dp,    0.0_dp]
    box%atoms(1:3, 3) = [0.5_dp,   0.866_dp,  0.0_dp]
    ! Bridge
    box%atoms(1:3, 4) = [4.5_dp,   0.0_dp,    0.0_dp]
    ! Sub-cluster B (triangle)
    box%atoms(1:3, 5) = [9.0_dp,   0.0_dp,    0.0_dp]
    box%atoms(1:3, 6) = [10.0_dp,  0.0_dp,    0.0_dp]
    box%atoms(1:3, 7) = [9.5_dp,   0.866_dp,  0.0_dp]

    crit%molType = 1
    crit%atomNum = 1
    crit%rCut = rCut
    crit%rCutSq = rCut * rCut
    call crit%Constructor(1)
    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for bridged cluster"
      nFail = nFail + 1
    endif

    ! --- Deleting the bridge should be rejected ---
    call CheckDeletion(crit, 1, 4, accept)
    if(accept) then
      write(nout, *) "  FAIL: deleting bridge atom 4 should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete bridge atom 4 -> rejected"
    endif

    ! --- Deleting any non-bridge atom should be accepted ---
    do iMol = 1, nMol
      if(iMol == 4) cycle
      call CheckDeletion(crit, 1, iMol, accept)
      if(.not. accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting non-bridge atom ", iMol, " should be accepted"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete non-bridge atom ", iMol, " -> accepted"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: DistCriteria Bridge <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: DistCriteria Bridge <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_distcrit_dense
!
!   Five atoms all within rCut = 5.0 of each other:
!
!     1 (0, 0, 0)     2 (2, 0, 0)     3 (1, 2, 0)
!     4 (3, 1, 0)     5 (1, 1, 0)
!
!   Max pairwise distance is 1-4 = sqrt(10) ~ 3.16, well within rCut.
!   Every pair is connected.  Removing any single atom leaves a fully
!   connected cluster (accept for all deletions).
!===========================================================================
  subroutine test_distcrit_dense(passed)
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol = 5
    real(dp), parameter :: rCut = 5.0E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteria) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iMol

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: DistCriteria - Fully Connected"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSingleTypeSystem(nMol, nMol, boxL)
    box => BoxArray(1)%box

    box%atoms(1:3, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    box%atoms(1:3, 2) = [2.0_dp, 0.0_dp, 0.0_dp]
    box%atoms(1:3, 3) = [1.0_dp, 2.0_dp, 0.0_dp]
    box%atoms(1:3, 4) = [3.0_dp, 1.0_dp, 0.0_dp]
    box%atoms(1:3, 5) = [1.0_dp, 1.0_dp, 0.0_dp]

    crit%molType = 1
    crit%atomNum = 1
    crit%rCut = rCut
    crit%rCutSq = rCut * rCut
    call crit%Constructor(1)
    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for dense cluster"
      nFail = nFail + 1
    endif

    ! --- Every deletion should be accepted ---
    do iMol = 1, nMol
      call CheckDeletion(crit, 1, iMol, accept)
      if(.not. accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting atom ", iMol, " should be accepted in dense cluster"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete atom ", iMol, " -> accepted"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: DistCriteria Dense <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: DistCriteria Dense <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_distcrit_addition
!
!   21 atoms from a simulation snapshot.  rCut = 1.01 (atom 17's nearest
!   neighbour is atom 6 at r=1.008; rCut must exceed this so that the
!   only break is between atoms 15 and 21 at r~1.27).
!
!   Tests:
!     - CheckInitialConstraint detects the broken cluster (reject).
!     - Adding a particle at the bridging position repairs the cluster
!       (accept).
!     - Adding a particle far from the gap does not repair it (reject).
!===========================================================================
  subroutine test_distcrit_addition(passed)
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol = 21
    integer, parameter :: nMolMax = 22
    real(dp), parameter :: rCut = 1.01E0_dp
    real(dp), parameter :: boxL = 20.0E0_dp

    type(DistCriteria) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: DistCriteria - Broken Cluster + Addition"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSingleTypeSystem(nMolMax, nMol, boxL)
    box => BoxArray(1)%box

    box%atoms(1:3,  1) = [ 0.419828907_dp,  3.084384358_dp, -0.696230034_dp]
    box%atoms(1:3,  2) = [-0.295913559_dp,  3.360484286_dp,  0.461977592_dp]
    box%atoms(1:3,  3) = [-0.295913559_dp,  3.360484286_dp,  0.461977592_dp]
    box%atoms(1:3,  4) = [-0.706711957_dp, -3.476433077_dp, -0.128838491_dp]
    box%atoms(1:3,  5) = [-0.544655113_dp, -1.027652238_dp,  0.315535769_dp]
    box%atoms(1:3,  6) = [ 0.427186735_dp, -2.125512257_dp, -0.406045023_dp]
    box%atoms(1:3,  7) = [-0.035134706_dp,  3.262694251_dp, -0.434337423_dp]
    box%atoms(1:3,  8) = [ 0.513069276_dp, -1.085374677_dp,  0.032627594_dp]
    box%atoms(1:3,  9) = [-0.758296760_dp, -2.511078959_dp, -0.142632491_dp]
    box%atoms(1:3, 10) = [-0.586434384_dp,  3.707181376_dp,  0.171187432_dp]
    box%atoms(1:3, 11) = [-0.530289580_dp, -0.069874442_dp,  0.163751815_dp]
    box%atoms(1:3, 12) = [ 0.223309618_dp,  0.288483741_dp,  0.213603682_dp]
    box%atoms(1:3, 13) = [ 0.090837114_dp, -2.529872010_dp, -0.545825420_dp]
    box%atoms(1:3, 14) = [-0.159075635_dp, -1.485773468_dp,  0.015596988_dp]
    box%atoms(1:3, 15) = [-0.432797679_dp,  2.418002946_dp, -0.388961245_dp]
    box%atoms(1:3, 16) = [ 0.513069276_dp, -1.085374677_dp,  0.032627594_dp]
    box%atoms(1:3, 17) = [ 1.145130606_dp, -2.613800444_dp,  0.105924577_dp]
    box%atoms(1:3, 18) = [ 0.430365456_dp,  1.025512795_dp,  0.238697703_dp]
    box%atoms(1:3, 19) = [-0.161693157_dp, -3.342407056_dp, -0.444230094_dp]
    box%atoms(1:3, 20) = [ 0.687332507_dp, -0.719665828_dp,  0.556729571_dp]
    box%atoms(1:3, 21) = [ 0.056786595_dp,  1.565591091_dp,  0.416862313_dp]

    crit%molType = 1
    crit%atomNum = 1
    crit%rCut = rCut
    crit%rCutSq = rCut * rCut
    call crit%Constructor(1)

    ! --- Initial check should FAIL (broken cluster) ---
    call crit%CheckInitialConstraint(box, initOK)
    if(initOK) then
      write(nout, *) "  FAIL: initial check should reject the broken cluster"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: initial check correctly rejects broken cluster"
    endif

    ! --- Adding a repair atom that bridges the gap should ACCEPT ---
    call CheckAddition(crit, 1, -0.188005542_dp, 1.991797019_dp, 0.013950534_dp, accept)
    if(.not. accept) then
      write(nout, *) "  FAIL: adding repair atom should be accepted"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: adding repair atom -> accepted"
    endif

    ! --- Adding a far-away atom that doesn't bridge the gap should REJECT ---
    call CheckAddition(crit, 1, 100.0_dp, 0.0_dp, 0.0_dp, accept)
    if(accept) then
      write(nout, *) "  FAIL: adding far-away atom should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: adding far-away atom -> rejected"
    endif

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: DistCriteria Addition <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: DistCriteria Addition <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
!                 DistCriteriaAllMol Tests
!===========================================================================
! test_allmol_chain
!
!   Five 1-atom molecules across two types arranged in a chain.
!   Type 1: 3 molecules, type 2: 2 molecules.
!   rCut(1) = rCut(2) = 3.0  ->  all pair cutoffs = 3.0
!
!   Atom layout (global molecule index in parentheses):
!
!     t1-m1(1) --- t2-m1(4) --- t1-m2(2) --- t2-m2(5) --- t1-m3(3)
!       0.0          2.9          5.8          8.7          11.6
!
!   NMolMax = [3, 2], NMol = [3, 2]
!   Global indices: type 1 -> 1,2,3; type 2 -> 4,5
!   Atom indices match global mol indices (1 atom per molecule).
!
!   Deleting ends (global 1 or 3) -> accept.
!   Deleting interior (global 4, 2, or 5) -> reject.
!===========================================================================
  subroutine test_allmol_chain(passed)
    use Common_MolInfo, only: nMolTypes
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol1 = 3, nMol2 = 2
    real(dp), parameter :: rCutVal = 3.0E0_dp
    real(dp), parameter :: spacing = 2.9E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteriaAllMol) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iType, jType
    real(dp) :: rMix

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: AllMolDistCriteria - Straight Chain"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupTwoTypeSystem(nMol1, nMol1, nMol2, nMol2, boxL)
    box => BoxArray(1)%box

    ! Place atoms: chain order is global 1, 4, 2, 5, 3
    box%atoms(1:3, 1) = [0.0_dp,          0.0_dp, 0.0_dp]  ! t1-m1
    box%atoms(1:3, 2) = [2*spacing,       0.0_dp, 0.0_dp]  ! t1-m2
    box%atoms(1:3, 3) = [4*spacing,       0.0_dp, 0.0_dp]  ! t1-m3
    box%atoms(1:3, 4) = [1*spacing,       0.0_dp, 0.0_dp]  ! t2-m1
    box%atoms(1:3, 5) = [3*spacing,       0.0_dp, 0.0_dp]  ! t2-m2

    ! Configure the constraint
    call crit%Constructor(1)
    crit%rCut(1) = rCutVal
    crit%rCut(2) = rCutVal
    crit%atomNumType(1) = 1
    crit%atomNumType(2) = 1
    do iType = 1, nMolTypes
      do jType = 1, nMolTypes
        rMix = 0.5E0_dp * (crit%rCut(iType) + crit%rCut(jType))
        crit%rCutPairSq(iType, jType) = rMix * rMix
      enddo
    enddo

    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for chain"
      nFail = nFail + 1
    endif

    ! --- End deletions (global 1 and 3) should be accepted ---
    call CheckDeletionByGlobal(crit, 1, accept)
    if(.not. accept) then
      write(nout, *) "  FAIL: deleting end (global 1) should be accepted"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete end global 1 -> accepted"
    endif

    call CheckDeletionByGlobal(crit, 3, accept)
    if(.not. accept) then
      write(nout, *) "  FAIL: deleting end (global 3) should be accepted"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete end global 3 -> accepted"
    endif

    ! --- Interior deletions (global 4, 2, 5) should be rejected ---
    call CheckDeletionByGlobal(crit, 4, accept)
    if(accept) then
      write(nout, *) "  FAIL: deleting interior (global 4) should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete interior global 4 -> rejected"
    endif

    call CheckDeletionByGlobal(crit, 2, accept)
    if(accept) then
      write(nout, *) "  FAIL: deleting interior (global 2) should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete interior global 2 -> rejected"
    endif

    call CheckDeletionByGlobal(crit, 5, accept)
    if(accept) then
      write(nout, *) "  FAIL: deleting interior (global 5) should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete interior global 5 -> rejected"
    endif

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: AllMol Chain <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: AllMol Chain <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_allmol_bridge
!
!   Seven 1-atom molecules (4 of type 1, 3 of type 2), two dense
!   sub-clusters connected through a single bridge atom.
!   rCut(1) = rCut(2) = 6.0  -> pair cutoff 6.0 for all pairs.
!
!   NMolMax = [4, 3], NMol = [4, 3]
!   Global indices: type 1 -> 1..4, type 2 -> 5..7
!
!   Sub-cluster A (globals 1, 2, 5):
!     t1-m1(1): (0, 0, 0)    t1-m2(2): (1, 0, 0)    t2-m1(5): (0.5, 0.866, 0)
!
!   Bridge (global 3):
!     t1-m3(3): (4.5, 0, 0)
!
!   Sub-cluster B (globals 4, 6, 7):
!     t1-m4(4): (9, 0, 0)    t2-m2(6): (10, 0, 0)   t2-m3(7): (9.5, 0.866, 0)
!
!   Deleting bridge (global 3) -> reject.
!   Deleting any other -> accept.
!===========================================================================
  subroutine test_allmol_bridge(passed)
    use Common_MolInfo, only: nMolTypes
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol1 = 4, nMol2 = 3
    real(dp), parameter :: rCutVal = 6.0E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteriaAllMol) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iGlob, bridgeGlob, iType, jType
    real(dp) :: rMix
    integer :: testMols(7)

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: AllMolDistCriteria - Bridged Cluster"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupTwoTypeSystem(nMol1, nMol1, nMol2, nMol2, boxL)
    box => BoxArray(1)%box

    ! Sub-cluster A
    box%atoms(1:3, 1) = [0.0_dp,  0.0_dp,   0.0_dp]   ! t1-m1, global 1
    box%atoms(1:3, 2) = [1.0_dp,  0.0_dp,   0.0_dp]   ! t1-m2, global 2
    ! Bridge
    box%atoms(1:3, 3) = [4.5_dp,  0.0_dp,   0.0_dp]   ! t1-m3, global 3
    ! Sub-cluster B
    box%atoms(1:3, 4) = [9.0_dp,  0.0_dp,   0.0_dp]   ! t1-m4, global 4
    ! Type 2 atoms (global 5,6,7)
    box%atoms(1:3, 5) = [0.5_dp,  0.866_dp, 0.0_dp]   ! t2-m1, global 5 (A)
    box%atoms(1:3, 6) = [10.0_dp, 0.0_dp,   0.0_dp]   ! t2-m2, global 6 (B)
    box%atoms(1:3, 7) = [9.5_dp,  0.866_dp, 0.0_dp]   ! t2-m3, global 7 (B)

    call crit%Constructor(1)
    crit%rCut(1) = rCutVal
    crit%rCut(2) = rCutVal
    crit%atomNumType(1) = 1
    crit%atomNumType(2) = 1
    do iType = 1, nMolTypes
      do jType = 1, nMolTypes
        rMix = 0.5E0_dp * (crit%rCut(iType) + crit%rCut(jType))
        crit%rCutPairSq(iType, jType) = rMix * rMix
      enddo
    enddo

    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for bridged cluster"
      nFail = nFail + 1
    endif

    bridgeGlob = 3
    testMols = [1, 2, 3, 4, 5, 6, 7]

    ! --- Deleting the bridge should be rejected ---
    call CheckDeletionByGlobal(crit, bridgeGlob, accept)
    if(accept) then
      write(nout, *) "  FAIL: deleting bridge (global 3) should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: delete bridge global 3 -> rejected"
    endif

    ! --- Deleting any non-bridge should be accepted ---
    do iGlob = 1, 7
      if(iGlob == bridgeGlob) cycle
      call CheckDeletionByGlobal(crit, iGlob, accept)
      if(.not. accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting non-bridge (global ", iGlob, ") should be accepted"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete non-bridge global ", iGlob, " -> accepted"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: AllMol Bridge <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: AllMol Bridge <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_allmol_dense
!
!   Five 1-atom molecules (3 of type 1, 2 of type 2) all within rCut = 5.0.
!   Every pair is connected; deleting any single atom still leaves a fully
!   connected cluster.
!
!   NMolMax = [3, 2], NMol = [3, 2]
!   Global indices: type 1 -> 1,2,3; type 2 -> 4,5
!
!     t1-m1(1): (0, 0, 0)     t1-m2(2): (2, 0, 0)    t1-m3(3): (1, 2, 0)
!     t2-m1(4): (3, 1, 0)     t2-m2(5): (1, 1, 0)
!===========================================================================
  subroutine test_allmol_dense(passed)
    use Common_MolInfo, only: nMolTypes
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol1 = 3, nMol2 = 2
    real(dp), parameter :: rCutVal = 5.0E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteriaAllMol) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iGlob, iType, jType
    real(dp) :: rMix

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: AllMolDistCriteria - Fully Connected"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupTwoTypeSystem(nMol1, nMol1, nMol2, nMol2, boxL)
    box => BoxArray(1)%box

    box%atoms(1:3, 1) = [0.0_dp, 0.0_dp, 0.0_dp]   ! t1-m1
    box%atoms(1:3, 2) = [2.0_dp, 0.0_dp, 0.0_dp]   ! t1-m2
    box%atoms(1:3, 3) = [1.0_dp, 2.0_dp, 0.0_dp]   ! t1-m3
    box%atoms(1:3, 4) = [3.0_dp, 1.0_dp, 0.0_dp]   ! t2-m1
    box%atoms(1:3, 5) = [1.0_dp, 1.0_dp, 0.0_dp]   ! t2-m2

    call crit%Constructor(1)
    crit%rCut(1) = rCutVal
    crit%rCut(2) = rCutVal
    crit%atomNumType(1) = 1
    crit%atomNumType(2) = 1
    do iType = 1, nMolTypes
      do jType = 1, nMolTypes
        rMix = 0.5E0_dp * (crit%rCut(iType) + crit%rCut(jType))
        crit%rCutPairSq(iType, jType) = rMix * rMix
      enddo
    enddo

    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for dense cluster"
      nFail = nFail + 1
    endif

    ! --- Every deletion should be accepted ---
    do iGlob = 1, nMol1 + nMol2
      call CheckDeletionByGlobal(crit, iGlob, accept)
      if(.not. accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting atom (global ", iGlob, ") should be accepted"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete atom global ", iGlob, " -> accepted"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: AllMol Dense <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: AllMol Dense <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_allmol_addition
!
!   Same 21-atom broken-cluster geometry as test_distcrit_addition, but
!   tested through the DistCriteriaAllMol constraint with a single
!   molecule type so that the AllMol addition path is exercised.
!===========================================================================
  subroutine test_allmol_addition(passed)
    use Common_MolInfo, only: nMolTypes
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol = 21
    integer, parameter :: nMolMax = 22
    real(dp), parameter :: rCut = 1.01E0_dp
    real(dp), parameter :: boxL = 20.0E0_dp

    type(DistCriteriaAllMol) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: AllMolDistCriteria - Broken Cluster + Addition"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSingleTypeSystem(nMolMax, nMol, boxL)
    box => BoxArray(1)%box

    box%atoms(1:3,  1) = [ 0.419828907_dp,  3.084384358_dp, -0.696230034_dp]
    box%atoms(1:3,  2) = [-0.295913559_dp,  3.360484286_dp,  0.461977592_dp]
    box%atoms(1:3,  3) = [-0.295913559_dp,  3.360484286_dp,  0.461977592_dp]
    box%atoms(1:3,  4) = [-0.706711957_dp, -3.476433077_dp, -0.128838491_dp]
    box%atoms(1:3,  5) = [-0.544655113_dp, -1.027652238_dp,  0.315535769_dp]
    box%atoms(1:3,  6) = [ 0.427186735_dp, -2.125512257_dp, -0.406045023_dp]
    box%atoms(1:3,  7) = [-0.035134706_dp,  3.262694251_dp, -0.434337423_dp]
    box%atoms(1:3,  8) = [ 0.513069276_dp, -1.085374677_dp,  0.032627594_dp]
    box%atoms(1:3,  9) = [-0.758296760_dp, -2.511078959_dp, -0.142632491_dp]
    box%atoms(1:3, 10) = [-0.586434384_dp,  3.707181376_dp,  0.171187432_dp]
    box%atoms(1:3, 11) = [-0.530289580_dp, -0.069874442_dp,  0.163751815_dp]
    box%atoms(1:3, 12) = [ 0.223309618_dp,  0.288483741_dp,  0.213603682_dp]
    box%atoms(1:3, 13) = [ 0.090837114_dp, -2.529872010_dp, -0.545825420_dp]
    box%atoms(1:3, 14) = [-0.159075635_dp, -1.485773468_dp,  0.015596988_dp]
    box%atoms(1:3, 15) = [-0.432797679_dp,  2.418002946_dp, -0.388961245_dp]
    box%atoms(1:3, 16) = [ 0.513069276_dp, -1.085374677_dp,  0.032627594_dp]
    box%atoms(1:3, 17) = [ 1.145130606_dp, -2.613800444_dp,  0.105924577_dp]
    box%atoms(1:3, 18) = [ 0.430365456_dp,  1.025512795_dp,  0.238697703_dp]
    box%atoms(1:3, 19) = [-0.161693157_dp, -3.342407056_dp, -0.444230094_dp]
    box%atoms(1:3, 20) = [ 0.687332507_dp, -0.719665828_dp,  0.556729571_dp]
    box%atoms(1:3, 21) = [ 0.056786595_dp,  1.565591091_dp,  0.416862313_dp]

    call crit%Constructor(1)
    crit%rCut(1) = rCut
    crit%atomNumType(1) = 1
    crit%rCutPairSq(1, 1) = rCut * rCut

    ! --- Initial check should FAIL (broken cluster) ---
    call crit%CheckInitialConstraint(box, initOK)
    if(initOK) then
      write(nout, *) "  FAIL: initial check should reject the broken cluster"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: initial check correctly rejects broken cluster"
    endif

    ! --- Adding a repair atom that bridges the gap should ACCEPT ---
    call CheckAddition(crit, 1, -0.188005542_dp, 1.991797019_dp, 0.013950534_dp, accept)
    if(.not. accept) then
      write(nout, *) "  FAIL: adding repair atom should be accepted"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: adding repair atom -> accepted"
    endif

    ! --- Adding a far-away atom that doesn't bridge the gap should REJECT ---
    call CheckAddition(crit, 1, 100.0_dp, 0.0_dp, 0.0_dp, accept)
    if(accept) then
      write(nout, *) "  FAIL: adding far-away atom should be rejected"
      nFail = nFail + 1
    else
      write(nout, *) "    ok: adding far-away atom -> rejected"
    endif

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: AllMol Addition <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: AllMol Addition <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_allmol_mixed_rcut_chain
!
!   Seven 1-atom molecules (4 of type 1, 3 of type 2) in a chain where
!   each link uses a DIFFERENT pair cutoff.
!
!   rCut(1) = 2.0,  rCut(2) = 4.0
!   Pair cutoffs:  t1-t1 = 2.0,  t1-t2 = 3.0,  t2-t2 = 4.0
!
!   NMolMax = [4, 3], NMol = [4, 3]
!   Global indices: type 1 -> 1,2,3,4; type 2 -> 5,6,7
!
!   Chain order (spatial, with link type and distance):
!
!     t1(G1)  --2.5--  t2(G5)  --3.5--  t2(G6)  --2.5--  t1(G2)  --1.5--  t1(G3)  --2.5--  t2(G7)  --2.5--  t1(G4)
!      0.0               2.5              6.0              8.5              10.0              12.5              15.0
!            t1-t2<3.0         t2-t2<4.0        t1-t2<3.0        t1-t1<2.0        t1-t2<3.0        t1-t2<3.0
!
!   No non-adjacent pairs are within their respective cutoff.
!
!   Deleting either end (G1 or G4) -> accept.
!   Deleting any interior atom (G5, G6, G2, G3, G7) -> reject.
!===========================================================================
  subroutine test_allmol_mixed_rcut_chain(passed)
    use Common_MolInfo, only: nMolTypes
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol1 = 4, nMol2 = 3
    real(dp), parameter :: rCut1 = 2.0E0_dp
    real(dp), parameter :: rCut2 = 4.0E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteriaAllMol) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iGlob, iType, jType
    real(dp) :: rMix

    integer, parameter :: nEnds = 2
    integer, parameter :: nInterior = 5
    integer :: ends(nEnds), interior(nInterior)
    integer :: k

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: AllMolDistCriteria - Mixed rCut Chain"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupTwoTypeSystem(nMol1, nMol1, nMol2, nMol2, boxL)
    box => BoxArray(1)%box

    ! Type 1 atoms (globals 1-4)
    box%atoms(1:3, 1) = [ 0.0_dp, 0.0_dp, 0.0_dp]   ! chain end
    box%atoms(1:3, 2) = [ 8.5_dp, 0.0_dp, 0.0_dp]
    box%atoms(1:3, 3) = [10.0_dp, 0.0_dp, 0.0_dp]
    box%atoms(1:3, 4) = [15.0_dp, 0.0_dp, 0.0_dp]   ! chain end
    ! Type 2 atoms (globals 5-7)
    box%atoms(1:3, 5) = [ 2.5_dp, 0.0_dp, 0.0_dp]
    box%atoms(1:3, 6) = [ 6.0_dp, 0.0_dp, 0.0_dp]
    box%atoms(1:3, 7) = [12.5_dp, 0.0_dp, 0.0_dp]

    call crit%Constructor(1)
    crit%rCut(1) = rCut1
    crit%rCut(2) = rCut2
    crit%atomNumType(1) = 1
    crit%atomNumType(2) = 1
    do iType = 1, nMolTypes
      do jType = 1, nMolTypes
        rMix = 0.5E0_dp * (crit%rCut(iType) + crit%rCut(jType))
        crit%rCutPairSq(iType, jType) = rMix * rMix
      enddo
    enddo

    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for mixed-rCut chain"
      nFail = nFail + 1
    endif

    ends     = [1, 4]
    interior = [5, 6, 2, 3, 7]

    ! --- End deletions should be accepted ---
    do k = 1, nEnds
      iGlob = ends(k)
      call CheckDeletionByGlobal(crit, iGlob, accept)
      if(.not. accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting end (global ", iGlob, ") should be accepted"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete end global ", iGlob, " -> accepted"
      endif
    enddo

    ! --- Interior deletions should be rejected ---
    do k = 1, nInterior
      iGlob = interior(k)
      call CheckDeletionByGlobal(crit, iGlob, accept)
      if(accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting interior (global ", iGlob, ") should be rejected"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete interior global ", iGlob, " -> rejected"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: AllMol Mixed rCut Chain <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: AllMol Mixed rCut Chain <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_allmol_alternating_chain
!
!   Seven 1-atom molecules (4 type-1, 3 type-2) with types strictly
!   alternating along the chain so that EVERY link is a t1-t2 pair.
!
!   rCut(1) = 2.0,  rCut(2) = 4.0
!   Pair cutoffs:  t1-t1 = 2.0,  t1-t2 = 3.0,  t2-t2 = 4.0
!
!   NMolMax = [4, 3], NMol = [4, 3]
!   Global indices: type 1 -> 1,2,3,4; type 2 -> 5,6,7
!
!   Chain order (all links are t1-t2 at distance 2.5 < 3.0):
!
!     t1(G1) --2.5-- t2(G5) --2.5-- t1(G2) --2.5-- t2(G6) --2.5-- t1(G3) --2.5-- t2(G7) --2.5-- t1(G4)
!      0.0             2.5            5.0             7.5            10.0            12.5            15.0
!
!   Non-adjacent same-type pairs (spacing >= 5.0) exceed their cutoff:
!     t1-t1: 5.0 > 2.0,  t2-t2: 5.0 > 4.0
!
!   Deleting either end (G1 or G4) -> accept.
!   Deleting any interior atom (G5, G2, G6, G3, G7) -> reject.
!===========================================================================
  subroutine test_allmol_alternating_chain(passed)
    use Common_MolInfo, only: nMolTypes
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol1 = 4, nMol2 = 3
    real(dp), parameter :: rCut1 = 2.0E0_dp
    real(dp), parameter :: rCut2 = 4.0E0_dp
    real(dp), parameter :: spacing = 2.5E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteriaAllMol) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK
    integer :: nFail, iGlob, iType, jType
    real(dp) :: rMix

    integer, parameter :: nEnds = 2
    integer, parameter :: nInterior = 5
    integer :: ends(nEnds), interior(nInterior)
    integer :: k

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: AllMolDistCriteria - Alternating Type Chain"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupTwoTypeSystem(nMol1, nMol1, nMol2, nMol2, boxL)
    box => BoxArray(1)%box

    ! Chain: t1 t2 t1 t2 t1 t2 t1, evenly spaced at 2.5
    ! Type 1 atoms (globals 1-4)
    box%atoms(1:3, 1) = [ 0.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos  0.0
    box%atoms(1:3, 2) = [ 2.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos  5.0
    box%atoms(1:3, 3) = [ 4.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos 10.0
    box%atoms(1:3, 4) = [ 6.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos 15.0
    ! Type 2 atoms (globals 5-7)
    box%atoms(1:3, 5) = [ 1.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos  2.5
    box%atoms(1:3, 6) = [ 3.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos  7.5
    box%atoms(1:3, 7) = [ 5.0_dp * spacing, 0.0_dp, 0.0_dp]  ! pos 12.5

    call crit%Constructor(1)
    crit%rCut(1) = rCut1
    crit%rCut(2) = rCut2
    crit%atomNumType(1) = 1
    crit%atomNumType(2) = 1
    do iType = 1, nMolTypes
      do jType = 1, nMolTypes
        rMix = 0.5E0_dp * (crit%rCut(iType) + crit%rCut(jType))
        crit%rCutPairSq(iType, jType) = rMix * rMix
      enddo
    enddo

    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass for alternating chain"
      nFail = nFail + 1
    endif

    ends     = [1, 4]
    interior = [5, 2, 6, 3, 7]

    ! --- End deletions should be accepted ---
    do k = 1, nEnds
      iGlob = ends(k)
      call CheckDeletionByGlobal(crit, iGlob, accept)
      if(.not. accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting end (global ", iGlob, ") should be accepted"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete end global ", iGlob, " -> accepted"
      endif
    enddo

    ! --- Interior deletions should be rejected ---
    do k = 1, nInterior
      iGlob = interior(k)
      call CheckDeletionByGlobal(crit, iGlob, accept)
      if(accept) then
        write(nout, "(A,I2,A)") "  FAIL: deleting interior (global ", iGlob, ") should be rejected"
        nFail = nFail + 1
      else
        write(nout, "(A,I2,A)") "    ok: delete interior global ", iGlob, " -> rejected"
      endif
    enddo

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: AllMol Alternating Chain <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: AllMol Alternating Chain <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
! test_distcrit_mc_loop
!
!   Simulates the Monte Carlo loop: repeatedly propose random displacement
!   moves, run DiffCheck/Update on the constraint, and after every accepted
!   move call CheckInitialConstraint from scratch to verify the incremental
!   topology stays consistent with the full calculation.
!
!   Uses a 5-atom chain with rCut = 3.0 and max displacement = 0.3.
!   This gives a healthy mix of accepted and rejected moves.
!
!   If the from-scratch check ever disagrees with the incremental state,
!   a topology drift has occurred and the test fails.
!===========================================================================
  subroutine test_distcrit_mc_loop(passed)
    use RandomGen, only: grnd, sgrnd
    use ParallelVar, only: nout
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nMol = 5
    integer, parameter :: nTrials = 500
    real(dp), parameter :: rCut = 3.0E0_dp
    real(dp), parameter :: spacing = 2.9E0_dp
    real(dp), parameter :: maxDisp = 0.3E0_dp
    real(dp), parameter :: boxL = 50.0E0_dp

    type(DistCriteria) :: crit
    class(SimpleBox), pointer :: box => null()
    logical :: accept, initOK, verifyOK
    integer :: nFail, nAccepted, nRejected
    integer :: iTrial, iMol, molIndx, atomIndx
    integer :: savedNout
    real(dp) :: dx, dy, dz
    type(Displacement) :: disp(1)

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: DistCriteria - MC Loop Safety Check"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    nAccepted = 0
    nRejected = 0

    call sgrnd(12345)

    call SetupSingleTypeSystem(nMol, nMol, boxL)
    box => BoxArray(1)%box

    do iMol = 1, nMol
      box%atoms(1, iMol) = (iMol - 1) * spacing
      box%atoms(2, iMol) = 0.0E0_dp
      box%atoms(3, iMol) = 0.0E0_dp
    enddo

    crit%molType = 1
    crit%atomNum = 1
    crit%rCut = rCut
    crit%rCutSq = rCut * rCut
    call crit%Constructor(1)
    call crit%CheckInitialConstraint(box, initOK)
    if(.not. initOK) then
      write(nout, *) "  FAIL: initial constraint check should pass"
      nFail = nFail + 1
      passed = .false.
      return
    endif

    ! Suppress verbose output from CheckInitialConstraint during the loop
    savedNout = nout
    open(unit=99, file='/dev/null', status='old', action='write')

    do iTrial = 1, nTrials
      iMol = int(nMol * grnd()) + 1
      if(iMol > nMol) iMol = nMol

      dx = maxDisp * (2.0_dp * grnd() - 1.0_dp)
      dy = maxDisp * (2.0_dp * grnd() - 1.0_dp)
      dz = maxDisp * (2.0_dp * grnd() - 1.0_dp)

      molIndx = box%MolGlobalIndx(1, iMol)
      atomIndx = box%MolStartIndx(molIndx)

      disp(1)%molType = 1
      disp(1)%molIndx = molIndx
      disp(1)%atmIndx = atomIndx
      disp(1)%x_new = box%atoms(1, atomIndx) + dx
      disp(1)%y_new = box%atoms(2, atomIndx) + dy
      disp(1)%z_new = box%atoms(3, atomIndx) + dz

      call crit%DiffCheck(box, disp, accept)

      if(accept) then
        call crit%Update(.true.)
        box%atoms(1, atomIndx) = disp(1)%x_new
        box%atoms(2, atomIndx) = disp(1)%y_new
        box%atoms(3, atomIndx) = disp(1)%z_new
        nAccepted = nAccepted + 1

        ! Verify from scratch (suppress output)
        nout = 99
        call crit%CheckInitialConstraint(box, verifyOK)
        nout = savedNout
        if(.not. verifyOK) then
          write(nout, "(A,I5,A)") "  FAIL: constraint drift at trial ", &
            iTrial, " (from-scratch check failed after accepted move)"
          nFail = nFail + 1
          exit
        endif
      else
        nRejected = nRejected + 1
      endif
    enddo

    close(99)

    write(nout, "(A,I5,A,I5,A,I5)") "    Trials: ", nTrials, &
      "  Accepted: ", nAccepted, "  Rejected: ", nRejected

    if(nFail == 0) then
      write(nout, *) "  >>> PASS: DistCriteria MC Loop Safety Check <<<"
      passed = .true.
    else
      write(nout, "(A,I4)") "   Failed checks: ", nFail
      write(nout, *) "  >>> FAIL: DistCriteria MC Loop Safety Check <<<"
      passed = .false.
    endif

  end subroutine

!===========================================================================
end program
!===========================================================================
