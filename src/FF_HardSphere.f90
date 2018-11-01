!================================================================================
module FF_HardSphere
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_HardSphere
    real(dp), allocatable :: rMin(:)
    real(dp), allocatable :: rMinTable(:,:)
    contains
      procedure, pass :: Constructor => Constructor_HardSphere
      procedure, pass :: DetailedECalc => Detailed_HardSphere
      procedure, pass :: DiffECalc => DiffECalc_HardSphere
      procedure, pass :: ShiftECalc_Single => Shift_HardSphere_Single
      procedure, pass :: NewECalc => New_HardSphere
      procedure, pass :: ProcessIO => ProcessIO_HardSphere
      procedure, pass :: Prologue => Prologue_HardSphere
      procedure, pass :: GetCutOff => GetCutOff_HardSphere
  end type

  contains
  !=============================================================================+
  subroutine Constructor_HardSphere(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%rMin(1:nAtomTypes), stat = AllocateStat)
    allocate(self%rMinTable(1:nAtomTypes, 1:nAtomTypes), stat = AllocateStat)

    self%rMin = 1.0E0_dp

    self%rMinTable = 0.5E0_dp
    self%rCut = 5E0_dp
    self%rCutSq = 5E0_dp**2

    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_HardSphere(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: E_LJ
    real(dp) :: rmin_ij      

    E_LJ = 0E0
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
        atmType2 = curbox % AtomType(jAtom)
        rmin_ij = self % rMinTable(atmType1,atmType2)          

        rx = curbox % atoms(1, iAtom)  -  curbox % atoms(1, jAtom)
        ry = curbox % atoms(2, iAtom)  -  curbox % atoms(2, jAtom)
        rz = curbox % atoms(3, iAtom)  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx**2 + ry**2 + rz**2
        if(rsq < rmin_ij) then
!          write(*,*) sqrt(rsq)
!          write(*,*) iAtom, jAtom
!          write(*,*) curbox%atoms(1,iAtom), curbox%atoms(2,iAtom), curbox%atoms(3,iAtom)
!          write(*,*) curbox%atoms(1,jAtom), curbox%atoms(2,jAtom), curbox%atoms(3,jAtom)
!          write(*,*) "ERROR! Overlaping atoms found in the current configuration!"
        endif
      enddo
    enddo
  
      
  end subroutine
!============================================================================
  subroutine DiffECalc_HardSphere(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
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
         call self % NewECalc(curbox, disp, tempList, tempNNei, E_Diff, accept)

      class is(Deletion)
         return

      class default
        write(*,*) "Unknown Perturbation Type."
    end select


  end subroutine
  !=====================================================================
  subroutine Shift_HardSphere_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
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

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      do jNei = 1, curbox%NeighList(1)%nNeigh(iAtom)
        jAtom = curbox%NeighList(1)%list(jNei, iAtom)
        atmType2 = curbox % AtomType(jAtom)
        rmin_ij = self % rMinTable(atmType2, atmType1)          

        rx = disp(iDisp)%x_new  -  curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new  -  curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new  -  curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < rmin_ij) then
          accept = .false.
          return
        endif 
      enddo
    enddo
 
  end subroutine
  !=====================================================================
  subroutine New_HardSphere(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig_sq
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.

    do iDisp = 1, dispLen
!      if(.not. disp(iDisp)%newAtom) then
!        cycle
!      endif
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      listIndx = disp(iDisp)%listIndex
!      if(disp(iDisp)%newlist) then
        maxNei = tempNNei(listIndx)
!      else
!        maxNei = curbox%NeighList(1)%nNeigh(listIndx)
!      endif

      do jNei = 1, maxNei
!        if(disp(iDisp)%newlist) then
          jAtom = tempList(jNei, listIndx)
!        else
!          jAtom = curbox%NeighList(1)%list(jNei, listIndx)
!        endif
        if( any(jAtom == disp(:)%atmIndx) ) then
          cycle
        endif

        atmType2 = curbox % AtomType(jAtom)
        rmin_ij = self%rMinTable(atmType2, atmType1)          

        rx = disp(iDisp)%x_new - curbox % atoms(1, jAtom)
        ry = disp(iDisp)%y_new - curbox % atoms(2, jAtom)
        rz = disp(iDisp)%z_new - curbox % atoms(3, jAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        if(rsq < rmin_ij) then
          accept = .false.
          return
        endif
      enddo
    enddo
  end subroutine
  !=====================================================================
  subroutine ProcessIO_HardSphere(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: rMin
  

    call GetXCommand(line, command, 1, lineStat)

    call CountCommands(line, nPar)
    select case(nPar)
      case(2)
        read(line, *) type1, rMin
        self%rMin(type1) = rMin
        do jType = 1, nAtomTypes
          self%rMinTable(type1, jType) = (rMin + self%rMin(jType))**2
          self%rMinTable(jType, type1) = (rMin + self%rMin(jType))**2
        enddo
      case(3)
        read(line, *) type1, type2, rMin
        self%rMinTable(type1, type2) = rMin**2
        self%rMinTable(type2, type1) = rMin**2

      case default
        lineStat = -1
    end select

  end subroutine
  !=============================================================================+
  function GetCutOff_HardSphere(self) result(rCut)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
    real(dp) :: rCut
    integer :: i, j

    rCut = -1E0_dp
    do i = 1, nAtomTypes
      do j = 1, nAtomTypes
        if(rCut < self%rMinTable(i,j)) then
          rCut = self%rMinTable(i,j)
        endif
      enddo
    enddo

    rCut = sqrt(rCut)

  end function
  !=====================================================================
  subroutine Prologue_HardSphere(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_HardSphere), intent(inout) :: self
    integer :: i, j

!    do i = 1, nAtomTypes
!      write(*,*) (self%rMinTable(i,j), j=1,nAtomTypes)
!    enddo



  end subroutine
  !=====================================================================
end module
!=====================================================================
