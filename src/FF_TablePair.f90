!================================================================================
! The easy pair forcefield class is designed to help users quickly add a 
! pairwise forcefield into the Classy code base. This provides a little less "fine control"
! compared to writing a separate pair style from scratch, but for most standard pairwise potentials
! it should yield the correct statistics.  However, similar computational efficiency to from scratch 
! is not guarenteed.
! To use apply the "extends(Pair_TableCut)" 
!================================================================================
module FF_Pair_TableCut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_TableCut
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)

    reap(dp), allocatable :: l
!    real(dp) :: rCut, rCutSq
    contains
!      procedure, pass :: Constructor => Constructor_TableCut
      procedure, pass :: DetailedECalc => Detailed_TableCut
      procedure, pass :: DiffECalc => DiffECalc_TableCut
      procedure, pass :: ShiftECalc_Single => Shift_TableCut_Single
      procedure, pass :: ShiftECalc_Multi => Shift_TableCut_Multi
      procedure, pass :: NewECalc => New_TableCut
      procedure, pass :: OldECalc => Old_TableCut
!      procedure, pass :: ProcessIO => ProcessIO_TableCut
      procedure, pass :: GetCutOff => GetCutOff_TableCut
      procedure, pass :: ComputePair => GetCutOff_TableCut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_TableCut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 0.5E0_dp
    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

!    IF (AllocateStat /= 0) STOP "Allocation in the Pair/Cut Pair Style"

  end subroutine
!============================================================================
  subroutine DiffECalc_TableCut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
!    class(displacement), intent(in) :: disp(:)
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

    select type(disp)
      class is(DisplacementNew)
         call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)

      class is(Addition)
!         write(*,*) size(tempList)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

!      class is(Displacement)
!        stop
      class default
        write(*,*) "Unknown Perturbation Type Encountered by the TableCut Pair Style."
    end select


  end subroutine
  !===================================================================================
  subroutine Detailed_TableCut(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: Pair
    real(dp) :: E_Pair
    real(dp) :: rmin_ij      

    E_Pair = 0E0
    curbox%ETable = 0E0
    accept = .true.
    do iAtom = 1, curbox%nMaxAtoms-1
      atmType1 = curbox % AtomType(iAtom)
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
        cycle
      endif
      do jAtom = iAtom+1, curbox%nMaxAtoms
        if( curbox%MolSubIndx(jAtom) > curbox%NMol(curbox%MolType(jAtom)) ) then
          cycle
        endif
        if( curbox%MolIndx(jAtom) == curbox%MolIndx(iAtom)  ) then
          cycle
        endif

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx**2 + ry**2 + rz**2
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self % rMinTable(atmType1, atmType2)          
          if(rsq < rmin_ij) then
            write(*,*) sqrt(rsq)
            write(*,*) iAtom, jAtom
            write(*,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
            write(*,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
            write(*,*) "ERROR! Overlaping atoms found in the current configuration!"
          endif 
          Pair = self%ComputePair(rsq, atmType1, atmType2)
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + Pair
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + Pair 
        endif
      enddo
    enddo
  
!    write(nout,*) "Lennard-Jones Energy:", E_Pair
      
    E_T = E_Pair    
  end subroutine
  !=====================================================================
  subroutine Shift_TableCut_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(DisplacementNew), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: Pair
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)

!      write(*,*) iAtom, curbox%NeighList(1)%nNeigh(iAtom)
!      write(*,*) iAtom, curbox%NeighList(1)%list(:, iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        atmType2 = curbox % AtomType(jAtom)
        if(rsq < self%rCutSq) then
          rmin_ij = self % rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif 

          Pair = self%ComputePair(rsq, atmType1, atmType2)
          E_Diff = E_Diff + Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + Pair
        endif

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          Pair = self%ComputePair(rsq, atmType1, atmType2)
          E_Diff = E_Diff - Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - Pair
        endif

      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine Shift_TableCut_Multi(self, curbox, disp, E_Diff)
    implicit none
      class(Pair_TableCut), intent(inout) :: self
      class(SimBox), intent(inout) :: curbox
      type(displacement), intent(in) :: disp(:)
      real(dp), intent(inout) :: E_Diff
   
  end subroutine
  !=====================================================================
  subroutine New_TableCut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: Pair
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.

!    write(*,*) "Length:", size(tempNNei)
!    write(*,*) "Length:", size(tempList)
!    write(*,*)
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
!      write(*,*) iAtom
      atmType1 = curbox % AtomType(iAtom)

      listIndx = disp(iDisp)%listIndex
      maxNei = tempNNei(listIndx)
      do jNei = 1, maxNei
        jAtom = tempList(jNei, listIndx)
        rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self%rMinTable(atmType2, atmType1)          
        
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
          Pair = self%ComputePair(rsq, atmType1, atmType2)
          E_Diff = E_Diff + Pair
!          write(*,*) iAtom, jAtom, rsq, Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + Pair
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_TableCut(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
!    type(displacement), intent(in) :: disp(:)
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2, globIndx
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: Pair
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

    globIndx = curBox % MolGlobalIndx(disp(1)%molType, disp(1)%molIndx)
    call curBox % GetMolData(globIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
      atmType1 = curbox % AtomType(iAtom) 
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)

        rx = curbox % atoms(1, iAtom) - curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom) - curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom) - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          Pair = self%ComputePair(rsq, atmType1, atmType2)
          E_Diff = E_Diff - Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - Pair
        endif
      enddo
    enddo
  end subroutine

  !=============================================================================+
  function GetCutOff_TableCut(self) result(rCut)
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
  !=============================================================================+
  function ComputePair_TableCut(self, rsq, atmType1, atmType2) result(energy)
    implicit none
    class(Pair_TableCut), intent(inout) :: self
    integer, intent(in) :: atmType1, atmType2
    real(dp), intent(in) :: rsq
    real(dp) :: energy

  end function
  !=====================================================================
end module
!=====================================================================
