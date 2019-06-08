!================================================================================
module FF_EasyPair_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: EasyPair_Cut
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor => Constructor_EasyPair_Cut
      procedure, pass :: PairFunction => PairFunction_EasyPair_Cut
      procedure, pass :: DetailedECalc => Detailed_EasyPair_Cut
      procedure, pass :: DiffECalc => DiffECalc_EasyPair_Cut
      procedure, pass :: ShiftECalc_Single => Shift_EasyPair_Cut_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_EasyPair_Cut_Multi
      procedure, pass :: NewECalc => New_EasyPair_Cut
      procedure, pass :: OldECalc => Old_EasyPair_Cut
      procedure, pass :: OrthoVolECalc => OrthoVol_EasyPair_Cut
      procedure, pass :: ProcessIO => ProcessIO_EasyPair_Cut
!      procedure, pass :: Prologue => Prologue_EasyPair_Cut
      procedure, pass :: GetCutOff => GetCutOff_EasyPair_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EasyPair_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer :: AllocateStat

    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "Allocation in the EasyPair_Cut Pair Style"

  end subroutine
  !=============================================================================+
  function PairFunction_EasyPair_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    E_Pair = 0E0_dp

  end function
  !===================================================================================
  subroutine Detailed_EasyPair_Cut(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: E_Pair
    real(dp) :: E_Total
    real(dp) :: rmin_ij      

    real(dp), pointer :: atoms(:,:) => null()

    call curbox%GetCoordinates(atoms)

    E_Total = 0E0
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

        rx = atoms(1, iAtom)  -  atoms(1, jAtom)
        ry = atoms(2, iAtom)  -  atoms(2, jAtom)
        rz = atoms(3, iAtom)  -  atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx**2 + ry**2 + rz**2
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self % rMinTable(atmType1, atmType2)          
          if(rsq < rmin_ij) then
            write(*,*) sqrt(rsq)
            write(*,*) iAtom, jAtom
            write(*,*) atoms(1,iAtom), atoms(2,iAtom), atoms(3,iAtom)
            write(*,*) atoms(1,jAtom), atoms(2,jAtom), atoms(3,jAtom)
            write(*,*) "ERROR! Overlaping atoms found in the current configuration!"
          endif 
          E_Total = E_Total + LJ
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + LJ
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_Pair 
        endif
      enddo
    enddo
  
!    write(nout,*) "Lennard-Jones Energy:", E_Total
      
    E_T = E_Total    
  end subroutine
!============================================================================
  subroutine DiffECalc_EasyPair_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept
    real(dp) :: E_Half

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp

    select type(disp)
      class is(Displacement)
         call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)

      class is(Addition)
!         write(*,*) size(tempList)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

      class is(OrthoVolChange)
         call self % OrthoVolECalc( curbox, disp, E_Diff, accept)

      class default
        write(*,*) "Unknown Perturbation Type Encountered by the EasyPair_Cut Pair Style."
        stop

    end select


  end subroutine

  !=====================================================================
  subroutine Shift_EasyPair_Cut_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)

!      write(*,*) iAtom, curbox%NeighList(1)%nNeigh(iAtom)
!      write(*,*) iAtom, curbox%NeighList(1)%list(:, iAtom)
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)

        rx = disp(iDisp)%x_new  - atoms(1, jAtom)
        ry = disp(iDisp)%y_new  - atoms(2, jAtom)
        rz = disp(iDisp)%z_new  - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) sqrt(rsq)

        atmType2 = curbox % AtomType(jAtom)
        if(rsq < self%rCutSq) then
          rmin_ij = self % rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif 
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
          E_Diff = E_Diff + E_Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + E_Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + E_Pair
        endif

        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
!        write(*,*) sqrt(rsq)
        if(rsq < self%rCutSq) then
          if(iAtom == jAtom) then
            write(*,*) "NeighborList Error!", iAtom, jAtom
          endif
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
          E_Diff = E_Diff - E_Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - E_Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - E_Pair
        endif

      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine New_EasyPair_Cut(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: E_Pair
    real(dp) :: rmin_ij      

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)


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
        rx = disp(iDisp)%x_new - atoms(1, jAtom)
        ry = disp(iDisp)%y_new - atoms(2, jAtom)
        rz = disp(iDisp)%z_new - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          atmType2 = curbox % AtomType(jAtom)
          rmin_ij = self%rMinTable(atmType2, atmType1)          
        
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
          E_Diff = E_Diff + E_Pair
!          write(*,*) iAtom, jAtom, rsq, E_Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + E_Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + E_Pair
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_EasyPair_Cut(self, curbox, disp, E_Diff)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: E_Pair, E_Pair2
    real(dp) :: rmin_ij      
    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)


    E_Diff = 0E0_dp

!    globIndx = curBox % MolGlobalIndx(disp(1)%molType, disp(1)%molIndx)
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
      atmType1 = curbox % AtomType(iAtom) 
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)
        ep = self % epsTable(atmType1, atmType2)
        sig_sq = self % sigTable(atmType1, atmType2)          
!        rmin_ij = self % rMinTable(atmType1, atmType2)          

        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
!          write(*,*) iAtom, jAtom, rsq, E_Pair
          E_Diff = E_Diff - E_Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - E_Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - E_Pair
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine OrthoVol_EasyPair_Cut(self, curbox, disp, E_Diff, accept)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(OrthoVolChange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iAtom, jAtom, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    integer :: molIndx1, molIndx2
    real(dp) :: dxi, dyi, dzi
    real(dp) :: dxj, dyj, dzj
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: rx2, ry2, rz2, rsq2
    real(dp) :: ep, sig_sq, ratio
    real(dp) :: LJ, LJ2
    real(dp) :: rmin_ij      

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)
    E_Diff = 0E0_dp
    do iAtom = 1, curbox%nMaxAtoms-1
      if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
        cycle
      endif
      atmType1 = curbox % AtomType(iAtom)
      molIndx1 = curbox % MolIndx(iAtom)
      dxi = curbox % centerMass(1, molIndx1) * (disp(1)%xScale-1E0_dp)
      dyi = curbox % centerMass(2, molIndx1) * (disp(1)%yScale-1E0_dp)
      dzi = curbox % centerMass(3, molIndx1) * (disp(1)%zScale-1E0_dp)
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)
        if(jAtom <= iAtom) then
          cycle
        endif
        molIndx2 = curbox % MolIndx(jAtom)
        dxj = curbox % centerMass(1,molIndx2) * (disp(1)%xScale-1E0_dp)
        dyj = curbox % centerMass(2,molIndx2) * (disp(1)%yScale-1E0_dp)
        dzj = curbox % centerMass(3,molIndx2) * (disp(1)%zScale-1E0_dp)

        rx =  atoms(1, iAtom) + dxi - atoms(1, jAtom) - dxj
        ry =  atoms(2, iAtom) + dyi - atoms(2, jAtom) - dyj
        rz =  atoms(3, iAtom) + dzi - atoms(3, jAtom) - dzj

        call curbox%BoundaryNew(rx, ry, rz, disp)
        rsq = rx*rx + ry*ry + rz*rz

        atmType2 = curbox % AtomType(jAtom)
        if(rsq < self%rCutSq) then
          rmin_ij = self%rMinTable(atmType2, atmType1)          
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
          E_Diff = E_Diff + E_Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) + E_Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + E_Pair
        endif


      enddo
    enddo

    E_Diff = E_Diff - curbox%ETotal
    curbox % dETable = curbox%dETable - curbox % ETable

  end subroutine
  !=============================================================================+
  function GetCutOff_EasyPair_Cut(self) result(rCut)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
  !=====================================================================
  subroutine Prologue_EasyPair_Cut(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer :: i, j


  end subroutine
  !=====================================================================
end module
!=====================================================================
