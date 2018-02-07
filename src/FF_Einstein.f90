!================================================================================
! Einstein Crystal Potential
module FF_Einstein
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_Einstein
!    real(dp) :: rCut, rCutSq
    real(dp) :: kSpring
    real(dp), allocatable :: initPos(:,:)
    contains
      procedure, pass :: Constructor => Constructor_Pair_Einstein
      procedure, pass :: DetailedECalc => Detailed_Pair_Einstein
      procedure, pass :: ShiftECalc_Single => Shift_Pair_Einstein_Single
      procedure, pass :: NewECalc => New_Pair_Einstein
      procedure, pass :: OldECalc => Old_Pair_Einstein
      procedure, pass :: ProcessIO => ProcessIO_Pair_Einstein
  end type

  contains
  !=============================================================================+
  subroutine Constructor_Pair_Einstein(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_Einstein), intent(inout) :: self
!    integer :: AllocateStat

!    IF (AllocateStat /= 0) STOP "*** Not enough memory ***"

  end subroutine
  !===================================================================================
  subroutine Detailed_Pair_Einstein(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_Einstein), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iAtom
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: E_ein

    accept = .true.
    E_Ein = 0E0_dp
    if(.not. allocated(self%initPos) ) then
      allocate( self%initPos(1:3, 1:curbox%nMaxAtoms) ) 
      self%initPos = 0E0_dp
      do iAtom = 1, curbox%nMaxAtoms
        if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
          cycle
        endif
        self%initPos(1, iAtom) = curBox%atoms(1, iAtom)
        self%initPos(2, iAtom) = curBox%atoms(2, iAtom)
        self%initPos(3, iAtom) = curBox%atoms(3, iAtom)
      enddo

    else
      do iAtom = 1, curbox%nMaxAtoms
        if( curbox%MolSubIndx(iAtom) > curbox%NMol(curbox%MolType(iAtom)) ) then
          cycle
        endif
        rx = self%initPos(1, iAtom) - curBox%atoms(1, iAtom)
        ry = self%initPos(2, iAtom) - curBox%atoms(2, iAtom)
        rz = self%initPos(3, iAtom) - curBox%atoms(3, iAtom)
        call curbox%Boundary(rx, ry, rz)
        rsq = rx*rx + ry*ry + rz*rz
        E_Ein = E_Ein + self%kSpring*rsq
      enddo
    endif
      
    write(nout, *) "Einstein Crystal Energy:", E_Ein
      
    E_T = E_Ein
  end subroutine
  !=====================================================================
  subroutine Shift_Pair_Einstein_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_Einstein), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
    real(dp) :: rx, ry, rz, rsq
    

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    curbox%dETable = 0E0_dp
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx

      rx = disp(iDisp)%x_new  -  self % initPos(1, iAtom)
      ry = disp(iDisp)%y_new  -  self % initPos(2, iAtom) 
      rz = disp(iDisp)%z_new  -  self % initPos(3, iAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      E_Diff = E_Diff + self%kSpring*rsq

      iAtom = disp(iDisp)%oldAtmIndx
      rx = curBox % atoms(1, iAtom)  -  self % initPos(1, iAtom)
      ry = curBox % atoms(2, iAtom)  -  self % initPos(2, iAtom)
      rz = curBox % atoms(3, iAtom)  -  self % initPos(3, iAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      E_Diff = E_Diff - self%kSpring*rsq

    enddo
 
  end subroutine
  !=====================================================================
  subroutine New_Pair_Einstein(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_Einstein), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx

      rx = disp(iDisp)%x_new  -  self % initPos(1, iAtom)
      ry = disp(iDisp)%y_new  -  self % initPos(2, iAtom) 
      rz = disp(iDisp)%z_new  -  self % initPos(3, iAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      E_Diff = E_Diff + self%kSpring*rsq
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_Pair_Einstein(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_Einstein), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, dispLen
    real(dp) :: rx, ry, rz, rsq

    E_Diff = 0E0_dp
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%oldAtmIndx

      rx = curBox % atoms(1, iAtom)  -  self % initPos(1, iAtom)
      ry = curBox % atoms(2, iAtom)  -  self % initPos(2, iAtom)
      rz = curBox % atoms(3, iAtom)  -  self % initPos(3, iAtom)
      call curbox%Boundary(rx, ry, rz)
      rsq = rx*rx + ry*ry + rz*rz
      E_Diff = E_Diff - self%kSpring*rsq
    enddo

  end subroutine
  !=====================================================================
  subroutine ProcessIO_Pair_Einstein(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: GetAllCommands, GetXCommand, maxLineLen
    implicit none
    class(Pair_Einstein), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30), allocatable :: parlist(:)
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2
    real(dp) :: k
  
    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("kspring")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) k
        self % kSpring = 0.5E0_dp * k

      case default
        lineStat = -1
    end select

  end subroutine
  !=====================================================================
end module
!=====================================================================
