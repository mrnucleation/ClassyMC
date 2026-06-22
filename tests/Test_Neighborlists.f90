!===========================================================================
! Test_Neighborlists
!
! Test cases which verify that the box's neighborlists do not become
! corrupted as various operations are applied to the box.  The biggest
! operation where the lists can get messed up is molecule deletion.
!
! When a molecule is deleted, ClassyMC fills the empty slot by taking the
! top molecule of the same type (the molecule of that type with the highest
! array index currently in the box) and copying it into the slot that was
! occupied by the deleted molecule.  This means every neighborlist must be
! updated such that:
!   1. All references to the deleted molecule's atoms are removed from the
!      lists of its former neighbors.
!   2. The top molecule's per-atom lists are moved into the deleted
!      molecule's atom slots.
!   3. Every atom that had the top molecule's atoms as neighbors has those
!      entries re-indexed to point at the new (shifted) array locations.
!
! Verification strategy (used after every deletion):
!   - A reference list is computed by brute force from the box's current
!     coordinates using the same pair rules as the list builder (skip
!     inactive molecules, skip same-molecule pairs, rsq <= rCutSq).
!   - The incrementally updated list must match the reference list for
!     every active atom (compared as sets, since ordering may differ).
!   - Inactive (vacated) atom slots must have empty lists, otherwise stale
!     entries can corrupt future additions.
!
! Tests:
!   test_deletion_singletype : One molecule type in the box.  Deletes the
!                              first, middle, and top molecules, then
!                              randomly deletes down to a single molecule.
!   test_deletion_multitype  : Several molecule types with differing atom
!                              counts, so the molecule indexing of each
!                              type block is offset.  Deletes from every
!                              type, completely empties the middle type,
!                              then randomly deletes from the remainder.
!===========================================================================
program test_neighborlist
    use VarPrecision
    use CoordinateTypes, only: Perturbation, SingleMol, Deletion
    use BoxData, only: BoxArray
    use SimpleSimBox, only: SimpleBox
    use OrthoBoxDef, only: OrthoBox
    use RSqListDef, only: RSqList
    use RandomGen, only: grnd, sgrnd
    use ParallelVar, only: nout

    implicit none

    logical :: passSingle, passMulti
    integer :: totalFailed

    nout = 6
    call sgrnd(42)

    write(nout, *)
    write(nout, *) "============================================================"
    write(nout, *) "       ClassyMC Test: Neighborlist Integrity"
    write(nout, *) "============================================================"

    call test_deletion_singletype(passSingle)
    call test_deletion_multitype(passMulti)

    write(nout, *)
    write(nout, *) "============================================================"
    write(nout, *) "                  Overall Summary"
    write(nout, *) "============================================================"
    totalFailed = 0
    if(passSingle) then
      write(nout, *) "   PASS  Deletion (Single Type)"
    else
      write(nout, *) " **FAIL  Deletion (Single Type)"
      totalFailed = totalFailed + 1
    endif
    if(passMulti) then
      write(nout, *) "   PASS  Deletion (Multi Type)"
    else
      write(nout, *) " **FAIL  Deletion (Multi Type)"
      totalFailed = totalFailed + 1
    endif
    write(nout, *)

    if(totalFailed == 0) then
      write(nout, *) "  >>> ALL TESTS PASSED <<<"
    else
      write(nout, "(A,I2,A)") "   >>> ", totalFailed, " test(s) FAILED <<<"
      error stop 1
    endif

  contains

!===========================================================================
! test_deletion_singletype
!
! Box contains a single molecule type (3 atoms per molecule).  Exercises
! the swap-with-top bookkeeping for the simplest type-indexing case.
!===========================================================================
  subroutine test_deletion_singletype(passed)
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nAtomsPerMol(1:1) = [3]
    integer, parameter :: NMolMax(1:1) = [30]
    integer, parameter :: NMol0(1:1)   = [20]
    integer, parameter :: NMolMin(1:1) = [1]
    real(dp), parameter :: boxL = 12.0E0_dp
    real(dp), parameter :: rCut = 4.0E0_dp

    class(SimpleBox), pointer :: testBox => null()
    integer :: nFail, subIndx

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: Deletion (Single Molecule Type)"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSystem(nAtomsPerMol, NMolMax, NMol0, NMolMin, boxL, rCut)
    testBox => BoxArray(1)%box

    !Sanity check that the freshly built list agrees with brute force.
    call VerifyNeighborList("initial build", nFail)

    !Delete a molecule from the middle of the type block.
    call DeleteAndCheck(1, 10, "delete middle molecule", nFail)

    !Delete the first molecule of the type block.
    call DeleteAndCheck(1, 1, "delete first molecule", nFail)

    !Delete the top molecule.  Here molIndx == topIndx, so the
    !swap-with-top step must safely collapse onto itself.
    call DeleteAndCheck(1, testBox%NMol(1), "delete top molecule", nFail)

    !Randomly delete molecules until we hit the minimum bound.  This
    !repeatedly mixes first/middle/top deletions at shrinking box sizes.
    do while(testBox%NMol(1) > NMolMin(1))
      subIndx = RandomActiveMol(1)
      call DeleteAndCheck(1, subIndx, "random deletion", nFail)
    enddo

    write(nout, *)
    write(nout, "(A,I4)") "   Failed checks: ", nFail
    if(nFail == 0) then
      write(nout, *) "  >>> PASS: Deletion (Single Type) <<<"
      passed = .true.
    else
      write(nout, *) "  >>> FAIL: Deletion (Single Type) <<<"
      passed = .false.
    endif

  end subroutine
!===========================================================================
! test_deletion_multitype
!
! Box contains three molecule types with different atom counts so that
! each type block has a different molecule/atom index offset.  Deleting
! from one type must shift only molecules within that type's block and
! leave the other blocks' neighbor entries intact.
!===========================================================================
  subroutine test_deletion_multitype(passed)
    implicit none
    logical, intent(out) :: passed

    integer, parameter :: nAtomsPerMol(1:3) = [1, 3, 2]
    integer, parameter :: NMolMax(1:3) = [15, 12, 10]
    integer, parameter :: NMol0(1:3)   = [10,  8,  6]
    integer, parameter :: NMolMin(1:3) = [ 0,  0,  0]
    real(dp), parameter :: boxL = 12.0E0_dp
    real(dp), parameter :: rCut = 4.0E0_dp

    class(SimpleBox), pointer :: testBox => null()
    integer :: nFail, iType, subIndx, nActive

    write(nout, *)
    write(nout, *) "------------------------------------------------------------"
    write(nout, *) "  Testing: Deletion (Multiple Molecule Types)"
    write(nout, *) "------------------------------------------------------------"

    nFail = 0
    call SetupSystem(nAtomsPerMol, NMolMax, NMol0, NMolMin, boxL, rCut)
    testBox => BoxArray(1)%box

    !Sanity check that the freshly built list agrees with brute force.
    call VerifyNeighborList("initial build", nFail)

    !Scripted deletions hitting first/middle/top slots of each type block.
    call DeleteAndCheck(2, 1,               "delete first molecule of type 2", nFail)
    call DeleteAndCheck(1, testBox%NMol(1), "delete top molecule of type 1",   nFail)
    call DeleteAndCheck(3, 3,               "delete middle molecule of type 3", nFail)
    call DeleteAndCheck(2, testBox%NMol(2), "delete top molecule of type 2",   nFail)
    call DeleteAndCheck(1, 1,               "delete first molecule of type 1", nFail)

    !Completely empty the middle type.  The other type blocks must remain
    !consistent even when an entire block in the middle becomes inactive.
    do while(testBox%NMol(2) > 0)
      subIndx = RandomActiveMol(2)
      call DeleteAndCheck(2, subIndx, "delete type 2", nFail)
    enddo

    !Randomly delete from the remaining types until only a couple of
    !molecules are left in the box.
    do
      nActive = sum(testBox%NMol(:))
      if(nActive <= 2) exit
      iType = RandomActiveType()
      subIndx = RandomActiveMol(iType)
      call DeleteAndCheck(iType, subIndx, "random multitype deletion", nFail)
    enddo

    write(nout, *)
    write(nout, "(A,I4)") "   Failed checks: ", nFail
    if(nFail == 0) then
      write(nout, *) "  >>> PASS: Deletion (Multi Type) <<<"
      passed = .true.
    else
      write(nout, *) "  >>> FAIL: Deletion (Multi Type) <<<"
      passed = .false.
    endif

  end subroutine
!===========================================================================
! SetupSystem
!
! Builds a fresh test system from scratch:
!   - Defines the molecule/atom type data in Common_MolInfo.
!   - Creates an OrthoBox in BoxArray(1) with the given molecule bounds.
!   - Assigns random coordinates to every active molecule.
!   - Attaches an RSqList neighborlist and builds it.
!===========================================================================
  subroutine SetupSystem(nAtomsPerMol, NMolMaxIn, NMol0, NMolMinIn, boxL, rCut)
    use Common_MolInfo, only: nMolTypes, nAtomTypes, MolData, AtomData
    implicit none
    integer, intent(in) :: nAtomsPerMol(:), NMolMaxIn(:), NMol0(:), NMolMinIn(:)
    real(dp), intent(in) :: boxL, rCut

    class(SimpleBox), pointer :: box => null()
    character(len=80) :: dimline
    integer :: iType, iMol, iAtom, atmType
    integer :: molIndx, molStart, lineStat

    !--- Molecule/atom type definitions ---
    nMolTypes = size(nAtomsPerMol)
    nAtomTypes = sum(nAtomsPerMol)
    if(allocated(MolData)) deallocate(MolData)
    if(allocated(AtomData)) deallocate(AtomData)
    allocate(MolData(1:nMolTypes))
    allocate(AtomData(1:nAtomTypes))
    AtomData(:)%Symb = "A"
    AtomData(:)%mass = 1.0E0_dp

    atmType = 0
    do iType = 1, nMolTypes
      MolData(iType)%nAtoms = nAtomsPerMol(iType)
      allocate(MolData(iType)%atomType(1:nAtomsPerMol(iType)))
      do iAtom = 1, nAtomsPerMol(iType)
        atmType = atmType + 1
        MolData(iType)%atomType(iAtom) = atmType
      enddo
    enddo

    !--- Box creation ---
    if(allocated(BoxArray)) deallocate(BoxArray)
    allocate(BoxArray(1:1))
    allocate(OrthoBox :: BoxArray(1)%box)
    box => BoxArray(1)%box
    box%boxID = 1

    call box%AllocateMolBound()
    box%NMolMin(:) = NMolMinIn(:)
    box%NMol(:)    = NMol0(:)
    box%NMolMax(:) = NMolMaxIn(:)
    call box%Constructor()

    write(dimline, "(A,3(F12.5,1x))") "dimension ", boxL, boxL, boxL
    call box%LoadDimension(dimline, lineStat)

    !--- Random coordinates for the active molecules ---
    do iType = 1, nMolTypes
      do iMol = 1, NMol0(iType)
        molIndx = box%MolGlobalIndx(iType, iMol)
        molStart = box%MolStartIndx(molIndx)
        box%atoms(1, molStart) = (grnd() - 0.5E0_dp)*boxL
        box%atoms(2, molStart) = (grnd() - 0.5E0_dp)*boxL
        box%atoms(3, molStart) = (grnd() - 0.5E0_dp)*boxL
        !Place the remaining atoms of the molecule near the first one.
        do iAtom = 2, nAtomsPerMol(iType)
          box%atoms(1, molStart+iAtom-1) = box%atoms(1, molStart) &
                                         + 0.6E0_dp*(grnd() - 0.5E0_dp)
          box%atoms(2, molStart+iAtom-1) = box%atoms(2, molStart) &
                                         + 0.6E0_dp*(grnd() - 0.5E0_dp)
          box%atoms(3, molStart+iAtom-1) = box%atoms(3, molStart) &
                                         + 0.6E0_dp*(grnd() - 0.5E0_dp)
        enddo
      enddo
    enddo

    !--- Neighborlist creation ---
    box%nLists = 1
    allocate(RSqList :: box%NeighList(1:1))
    call box%NeighList(1)%Constructor(1, rCut)
    call box%NeighList(1)%BuildList(1)

  end subroutine
!===========================================================================
! DeleteAndCheck
!
! Deletes the subIndx-th molecule of type iType the same way an accepted
! MC deletion move does (via box%DeleteMol, which shifts the top molecule
! of that type into the vacated slot and updates all neighborlists), then
! verifies the lists against a brute-force reference.
!===========================================================================
  subroutine DeleteAndCheck(iType, subIndx, label, nFail)
    implicit none
    integer, intent(in) :: iType, subIndx
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nFail

    type(Deletion) :: delMove
    class(SimpleBox), pointer :: box => null()
    integer :: molIndx

    box => BoxArray(1)%box
    molIndx = box%MolGlobalIndx(iType, subIndx)

    !Fill out the perturbation the same way MCMove_Delete does.  The
    !neighborlist update itself is triggered by DeleteMol on acceptance.
    delMove%molType = iType
    delMove%molIndx = molIndx
    delMove%atmIndx = box%MolStartIndx(molIndx)

    call box%DeleteMol(molIndx)
    call VerifyNeighborList(label, nFail)

  end subroutine
!===========================================================================
! VerifyNeighborList
!
! Recomputes the neighborlist by brute force from the box's current
! coordinates and compares it (as per-atom sets) against the incrementally
! maintained list.  Inactive atom slots must have empty lists.
!===========================================================================
  subroutine VerifyNeighborList(label, nFail)
    implicit none
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nFail

    class(SimpleBox), pointer :: box => null()
    integer, allocatable :: refList(:,:), refNNei(:)
    integer :: iAtom, jAtom, iNei, actN
    real(dp) :: rx, ry, rz, rsq, rCutSq
    logical :: okay, atomOK

    box => BoxArray(1)%box
    rCutSq = box%NeighList(1)%rCutSq

    allocate(refList(1:box%nMaxAtoms, 1:box%nMaxAtoms))
    allocate(refNNei(1:box%nMaxAtoms))
    refList = 0
    refNNei = 0

    !--- Brute force reference using the same pair rules as the builder ---
    do iAtom = 1, box%nMaxAtoms-1
      if(.not. box%IsActive(iAtom)) cycle
      do jAtom = iAtom+1, box%nMaxAtoms
        if(.not. box%IsActive(jAtom)) cycle
        if(box%MolIndx(iAtom) == box%MolIndx(jAtom)) cycle
        rx = box%atoms(1, iAtom) - box%atoms(1, jAtom)
        ry = box%atoms(2, iAtom) - box%atoms(2, jAtom)
        rz = box%atoms(3, iAtom) - box%atoms(3, jAtom)
        call box%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq <= rCutSq) then
          refNNei(iAtom) = refNNei(iAtom) + 1
          refList(refNNei(iAtom), iAtom) = jAtom
          refNNei(jAtom) = refNNei(jAtom) + 1
          refList(refNNei(jAtom), jAtom) = iAtom
        endif
      enddo
    enddo

    !--- Compare against the incrementally updated list ---
    okay = .true.
    do iAtom = 1, box%nMaxAtoms
      actN = box%NeighList(1)%nNeigh(iAtom)

      !Vacated/inactive slots must not retain stale neighbor entries.
      if(.not. box%IsActive(iAtom)) then
        if(actN /= 0) then
          okay = .false.
          write(nout, "(3A)")  "  FAIL (", trim(label), ")"
          write(nout, "(A,I5,A,I4)") "    Inactive atom ", iAtom, &
                       " still has neighbors! nNeigh =", actN
        endif
        cycle
      endif

      atomOK = (actN == refNNei(iAtom))
      if(atomOK) then
        !Set comparison: every reference neighbor must be in the actual
        !list and vice versa (counts already match).
        do iNei = 1, refNNei(iAtom)
          if(.not. any(box%NeighList(1)%list(1:actN, iAtom) == refList(iNei, iAtom))) then
            atomOK = .false.
            exit
          endif
        enddo
        if(atomOK) then
          do iNei = 1, actN
            if(.not. any(refList(1:refNNei(iAtom), iAtom) == box%NeighList(1)%list(iNei, iAtom))) then
              atomOK = .false.
              exit
            endif
          enddo
        endif
      endif

      if(.not. atomOK) then
        okay = .false.
        write(nout, "(3A)")  "  FAIL (", trim(label), ")"
        write(nout, "(A,I5)")        "    Atom: ", iAtom
        write(nout, "(A,1000(I5,1x))") "    Expected: ", &
              (refList(iNei, iAtom), iNei = 1, refNNei(iAtom))
        write(nout, "(A,1000(I5,1x))") "    Actual:   ", &
              (box%NeighList(1)%list(iNei, iAtom), iNei = 1, actN)
      endif
    enddo

    if(okay) then
      write(nout, "(3A,1000(I3,1x))") "    ok: ", label, "   NMol:", box%NMol(:)
    else
      nFail = nFail + 1
    endif

    deallocate(refList)
    deallocate(refNNei)

  end subroutine
!===========================================================================
! RandomActiveMol
!
! Returns a random molecule sub-index (1..NMol) for the given type.
!===========================================================================
  function RandomActiveMol(iType) result(subIndx)
    implicit none
    integer, intent(in) :: iType
    integer :: subIndx

    subIndx = floor(grnd()*BoxArray(1)%box%NMol(iType)) + 1
    subIndx = min(subIndx, int(BoxArray(1)%box%NMol(iType)))

  end function
!===========================================================================
! RandomActiveType
!
! Returns a random molecule type that still has molecules which can be
! deleted (NMol > NMolMin).
!===========================================================================
  function RandomActiveType() result(iType)
    use Common_MolInfo, only: nMolTypes
    implicit none
    integer :: iType

    do
      iType = floor(grnd()*nMolTypes) + 1
      iType = min(iType, nMolTypes)
      if(BoxArray(1)%box%NMol(iType) > BoxArray(1)%box%NMolMin(iType)) exit
    enddo

  end function
!===========================================================================
end program
!===========================================================================
