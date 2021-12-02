!================================================================================
! Easy Pair Template for general pair forcefields. The purpose of this class is
! to act as an inheritable object that can allow the user to quickly set up
! any standard pair distance based forcefields by simply creating a child object
! from this class. By overwriting the PairFunction class and parameter settings
! the procedures for all cut-off based forcefields are automatically defined.
!================================================================================
module FF_EasyPair_Cut
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: EasyPair_Cut
    logical :: usetailcorrection = .false.
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor => Constructor_EasyPair_Cut
      procedure, pass :: PairFunction => PairFunction_EasyPair_Cut
      procedure, pass :: TailCorrection => TailCorrection_EasyPair_Cut
      procedure, pass :: DetailedECalc => Detailed_EasyPair_Cut
      procedure, pass :: DiffECalc => DiffECalc_EasyPair_Cut
      procedure, pass :: ShiftECalc_Single => Shift_EasyPair_Cut_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_EasyPair_Cut_Multi
      procedure, pass :: NewECalc => New_EasyPair_Cut
      procedure, pass :: OldECalc => Old_EasyPair_Cut
      procedure, pass :: OrthoVolECalc => OrthoVol_EasyPair_Cut
      procedure, pass :: AtomExchange => AtomExchange_EasyPair_Cut
      procedure, pass :: SinglePair => SinglePair_EasyPair_Cut
      procedure, pass :: ManyBody => ManyBody_EasyPair_Cut
      procedure, pass :: ProcessIO => ProcessIO_EasyPair_Cut
!      procedure, pass :: Prologue => Prologue_EasyPair_Cut
      procedure, pass :: GetCutOff => GetCutOff_EasyPair_Cut
  end type

  contains
  !=============================================================================+
  subroutine Constructor_EasyPair_Cut(self)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer :: AllocateStat

    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) error STOP "Allocation in the EasyPair_Cut Pair Style"

  end subroutine
  !======================================================================+
  function PairFunction_EasyPair_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    E_Pair = 0E0_dp

  end function
  !======================================================================+
  subroutine TailCorrection_EasyPair_Cut(self, curbox, E_Corr, disp) 
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(Perturbation), intent(in), optional :: disp(:)
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(out) :: E_Corr

    E_Corr = 0E0_dp

  end subroutine
  !=======================================================================
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
    real(dp) :: E_Pair
    real(dp) :: E_Total, E_Corr
    real(dp) :: rmin_ij      

    real(dp), pointer :: atoms(:,:) => null()

    call curbox%GetCoordinates(atoms)

    E_Total = 0E0_dp
    curbox%ETable = 0E0_dp
    accept = .true.
    do iAtom = 1, curbox%nMaxAtoms-1
      atmType1 = curbox % AtomType(iAtom)
      if( .not. curbox%IsActive(iAtom) ) then
        cycle
      endif
      do jAtom = iAtom+1, curbox%nMaxAtoms
        if( .not. curbox%IsActive(jAtom) ) then
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
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2) 
          E_Total = E_Total + E_Pair
          curbox%ETable(iAtom) = curbox%ETable(iAtom) + E_Pair
          curbox%ETable(jAtom) = curbox%ETable(jAtom) + E_Pair 
        endif
      enddo
    enddo
  
    write(nout,*) "Total Pair Energy:", E_Total
      

    if( self%usetailcorrection ) then
      E_Corr = 0E0_dp
      call self%TailCorrection( curbox, E_Corr )
      E_Total = E_Total + E_Corr
      write(nout,*) "Total Tail Corrections:", E_Corr
    endif
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
    real(dp) :: E_Half, E_Corr

    accept = .true.
    curbox % dETable = 0E0_dp
    E_Diff = 0E0_dp
!       write(*,*) self%RCut, self%RCutSq

    select type(disp)
      class is(Displacement)
        if(disp(1)%newlist) then
          call self % ShiftECalc_Single(curbox, disp, E_Diff, accept, tempList, tempNNei)
        else
          call self % ShiftECalc_Single(curbox, disp, E_Diff, accept)
        endif

      class is(Addition)
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         call self % OldECalc(curbox, disp, E_Diff)

      class is(OrthoVolChange)
         call self % OrthoVolECalc( curbox, disp, E_Diff, accept)

      class is(AtomExchange)
         call self % AtomExchange( curbox, disp, E_Diff, accept)

      class default
        write(0,*) "Unknown Perturbation Type Encountered by the EasyPair_Cut Pair Style."
        error stop

    end select

    if( self%usetailcorrection ) then
      E_Corr = 0E0_dp
      call self%TailCorrection( curbox, E_Corr, disp )
      E_Diff = E_Diff + E_Corr
    endif

  end subroutine
  !=====================================================================
  subroutine Shift_EasyPair_Cut_Single(self, curbox, disp, E_Diff, accept, tempList, tempNNei)
    USE OMP_LIB
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Displacement), intent(in) :: disp(:)
    integer, intent(in), optional, target :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: E_Pair
    real(dp) :: rmin_ij

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    if( present(tempList) .and. present(tempNNei) ) then
      neighlist => tempList
      nNeigh => tempNNei
    else
      call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
      call curbox%GetCoordinates(atoms)
    endif

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)

        rx = disp(iDisp)%x_new  - atoms(1, jAtom)
        ry = disp(iDisp)%y_new  - atoms(2, jAtom)
        rz = disp(iDisp)%z_new  - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
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
        if(rsq < self%rCutSq) then
!          if(iAtom == jAtom) then
!            write(*,*) "NeighborList Error!", iAtom, jAtom
!          endif
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
    real(dp) :: E_Pair, E_Pair2
    real(dp) :: rmin_ij      
    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)


    E_Diff = 0E0_dp
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
      atmType1 = curbox % AtomType(iAtom) 
      do jNei = 1, nNeigh(iAtom)
        jAtom = neighlist(jNei, iAtom)

        atmType2 = curbox % AtomType(jAtom)
        rx = atoms(1, iAtom) - atoms(1, jAtom)
        ry = atoms(2, iAtom) - atoms(2, jAtom)
        rz = atoms(3, iAtom) - atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
          E_Diff = E_Diff - E_Pair
          curbox % dETable(iAtom) = curbox % dETable(iAtom) - E_Pair
          curbox % dETable(jAtom) = curbox % dETable(jAtom) - E_Pair
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine OrthoVol_EasyPair_Cut(self, curbox, disp, E_Diff, accept)
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(OrthoVolChange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iAtom, jAtom, maxNei, listIndx, jNei
    integer :: iType, typeStart, typeEnd
    integer :: atmType1, atmType2
    integer :: molIndx1, molIndx2
    real(dp) :: iTemp(1:3), jTemp(1:3)
    real(dp) :: dxi, dyi, dzi
    real(dp) :: dxj, dyj, dzj
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: rx2, ry2, rz2, rsq2
    real(dp) :: rmin_ij
    real(dp) :: E_Pair

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)
    E_Diff = 0E0_dp
   do iType = 1, nMolTypes
      call curbox%GetTypeAtoms(iType, typeStart, typeEnd)
      if(typeStart < 1) cycle

      do iAtom = typeStart, typeEnd
        atmType1 = curbox % AtomType(iAtom)
        molIndx1 = curbox % MolIndx(iAtom)
        dxi = curbox % centerMass(1, molIndx1) * (disp(1)%xScale-1E0_dp)
        dyi = curbox % centerMass(2, molIndx1) * (disp(1)%yScale-1E0_dp)
        dzi = curbox % centerMass(3, molIndx1) * (disp(1)%zScale-1E0_dp)
        iTemp(1:3) = atoms(1:3,iAtom) - curbox % centerMass(1:3, molIndx1)
        call curbox%Boundary(iTemp(1), iTemp(2), iTemp(3))
        iTemp(1:3) = iTemp(1:3) + curbox % centerMass(1:3, molIndx1)
        do jNei = 1, nNeigh(iAtom)
          jAtom = neighlist(jNei, iAtom)
          if(jAtom <= iAtom) then
            cycle
          endif
          molIndx2 = curbox % MolIndx(jAtom)
          dxj = curbox % centerMass(1,molIndx2) * (disp(1)%xScale-1E0_dp)
          dyj = curbox % centerMass(2,molIndx2) * (disp(1)%yScale-1E0_dp)
          dzj = curbox % centerMass(3,molIndx2) * (disp(1)%zScale-1E0_dp)
          jTemp(1:3) = atoms(1:3,jAtom) - curbox % centerMass(1:3, molindx2)
          call curbox%Boundary(jTemp(1), jTemp(2), jTemp(3))
          jTemp(1:3) = jTemp(1:3) + curbox % centerMass(1:3, molindx2)
          rx =  iTemp(1) + dxi - jTemp(1) - dxj
          ry =  iTemp(2) + dyi - jTemp(2) - dyj
          rz =  iTemp(3) + dzi - jTemp(3) - dzj

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
    enddo

    E_Diff = E_Diff - curbox%E_Inter
    curbox % dETable = curbox%dETable - curbox % ETable

  end subroutine
  !=====================================================================
  subroutine AtomExchange_EasyPair_Cut(self, curbox, disp, E_Diff, accept)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(AtomExchange), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtomNew, iAtomOld, jAtom, remLen, jNei
    integer :: newType1, oldType1
    integer :: atmType2, globIndx
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq, r
    real(dp) :: E_Pair_New, E_Pair_Old
    real(dp) :: rmin_ij      

    integer, pointer :: nNeigh(:) => null()
    integer, pointer :: neighlist(:,:) => null()
    real(dp), pointer :: atoms(:,:) => null()

    call curbox%Neighlist(1)%GetListArray(neighlist, nNeigh)
    call curbox%GetCoordinates(atoms)

    E_Diff = 0E0_dp


    iAtomNew = disp(1) % newAtmIndx
    iAtomOld = disp(1) % oldAtmIndx
    newType1 = curbox % AtomType(iAtomNew)
    oldType1 = curbox % AtomType(iAtomOld) 
    do jNei = 1, nNeigh(iAtomOld)
        jAtom = neighlist(jNei, iAtomOld)
        atmType2 = curbox % AtomType(jAtom)

        rx = curbox % atoms(1, iAtomOld) - curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtomOld) - curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtomOld) - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < self%rCutSq) then
          rmin_ij = self%rMinTable(atmType2, newType1)
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
          E_Pair_new = self%PairFunction(rsq, newType1, atmType2)
          E_Pair_Old = self%PairFunction(rsq, oldType1, atmType2)

          E_Diff = E_Diff + E_Pair_New - E_Pair_Old
          curbox % dETable(iAtomOld) = curbox % dETable(iAtomOld) + E_Pair_New - E_Pair_Old
          curbox % dETable(jAtom) = curbox % dETable(jAtom) + E_Pair_New - E_Pair_Old
        endif
    enddo
  end subroutine
  !=============================================================================+
  function GetCutOff_EasyPair_Cut(self) result(rCut)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    real(dp) :: rCut

    rCut = self%rCut
  end function
!=============================================================================+
  function SinglePair_EasyPair_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)
  end function
!=============================================================================+
  function SinglePair_Approx_EasyPair_Cut(self, rsq, atmtype1, atmtype2) result(E_Pair)
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    integer, intent(in) :: atmtype1, atmtype2
    real(dp), intent(in) :: rsq
    real(dp) :: E_Pair

    E_Pair = self%PairFunction(rsq, atmtype1, atmtype2)

  end function
!=============================================================================+
  subroutine ManyBody_EasyPair_Cut(self, curbox, atmtype1, pos1, atmtypes, posN, E_Many, accept )
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    integer, intent(in) :: atmtype1
    integer, intent(in) :: atmtypes(:)
    real(dp), intent(in) :: pos1(:)
    real(dp), intent(in) :: posN(:,:)
    logical, intent(out) :: accept
    real(dp), intent(out) :: E_Many


    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType2
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: E_Pair
    real(dp) :: rmin_ij      

    accept = .true.
    E_Many = 0E0_dp
    do jAtom = 1, size(posN, 2)
      atmType2 = atmtypes(jAtom)
      rx = pos1(1) - posN(1, jAtom)
      ry = pos1(2) - posN(2, jAtom)
      rz = pos1(3) - posN(3, jAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      if(rsq < self%rCutSq) then
          rmin_ij = self%rMinTable(atmType2, atmtype1)
          if(rsq < rmin_ij) then
            accept = .false.
            return
          endif
        E_Pair = self % PairFunction(rsq, atmtype1, atmtype2)
        E_Many = E_Many + E_Pair
      endif
    enddo

  end subroutine
  !=====================================================================
  subroutine ProcessIO_EasyPair_Cut(self, line)
    use Input_Format, only: maxLineLen
    implicit none
    class(EasyPair_Cut), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line

  end subroutine
  !=====================================================================
end module
!=====================================================================
