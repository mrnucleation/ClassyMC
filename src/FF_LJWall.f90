!================================================================================
module FF_Pair_LJ_Wall
  use Template_ForceField, only: ForceField
  use VarPrecision
  use Template_SimBox, only: SimBox
  use CoordinateTypes

  type, extends(forcefield) :: Pair_LJ_Wall
    real(dp), allocatable :: eps_star(:)
    real(dp), allocatable :: sig_star(:)
    real(dp) :: zWallPos

!    real(dp) :: rCut, rCutSq
    contains
      procedure, pass :: Constructor => Constructor_LJ_Wall
      procedure, pass :: DetailedECalc => Detailed_LJ_Wall
      procedure, pass :: DiffECalc => DiffECalc_LJ_Wall
      procedure, pass :: ShiftECalc_Single => Shift_LJ_Wall_Single
!      procedure, pass :: ShiftECalc_Multi => Shift_LJ_Wall_Multi
      procedure, pass :: NewECalc => New_LJ_Wall
      procedure, pass :: OldECalc => Old_LJ_Wall
      procedure, pass :: OrthoVolECalc => OrthoVol_LJ_Wall
      procedure, pass :: ProcessIO => ProcessIO_LJ_Wall
      procedure, pass :: Prologue => Prologue_LJ_Wall
      procedure, pass :: GetCutOff => GetCutOff_LJ_Wall
  end type

  real(dp), parameter :: coeff1 = 2.0E0_dp/15E0_dp

  contains
  !=============================================================================+
  subroutine Constructor_LJ_Wall(self)
    use Common_MolInfo, only: nAtomTypes
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    integer :: AllocateStat

    allocate(self%eps_star(1:nAtomTypes), stat = AllocateStat)
    allocate(self%sig_star(1:nAtomTypes), stat = AllocateStat)

    self%eps_star = 4E0_dp
    self%sig_star = 1E0_dp

    self%rCut = 0E0_dp
    self%rCutSq = 0E0_dp

    IF (AllocateStat /= 0) STOP "Allocation in the LJ/Wall Pair Style"

  end subroutine
  !===================================================================================
  subroutine Detailed_LJ_Wall(self, curbox, E_T, accept)
    use ParallelVar, only: nout
    use Common_MolInfo, only: nMolTypes
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    real(dp), intent(inOut) :: E_T
    logical, intent(out) :: accept
    integer :: iType, jType, iAtom, jAtom
    integer :: iLow, iUp, jLow, jUp
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig
    real(dp) :: LJ
    real(dp) :: E_LJ
    real(dp) :: rmin_ij      

    E_LJ = 0E0
!    curbox%ETable = 0E0
    accept = .true.
    do iAtom = 1, curbox%nMaxAtoms
      atmType1 = curbox % AtomType(iAtom)
      if( curbox%IsActive(iAtom) ) then
        cycle
      endif
      rz = curbox % atoms(3, iAtom) - self%zWallPos
!      call curbox%Boundary(rx, ry, rz)
      call curbox%Boundary(rz=rz)
      if(rz < self%rCutSq) then
      endif
      ep = self % eps_star(atmType1)
      sig = self % sig_star(atmType1)
      LJ = (sig/rz)**3
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_LJ = E_LJ + LJ
    enddo
  
    write(nout,*) "Lennard-Jones Wall Energy:", E_LJ
      
    E_T = E_LJ    
  end subroutine
!============================================================================
  subroutine DiffECalc_LJ_Wall(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    class(simBox), intent(inout) :: curbox
    class(Perturbation), intent(inout) :: disp(:)
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
        write(*,*) "Unknown Perturbation Type Encountered by the LJ_Wall Pair Style."
        stop

    end select


  end subroutine

  !=====================================================================
  subroutine Shift_LJ_Wall_Single(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Displacement), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jNei, jAtom, dispLen
!    integer :: maxIndx, minIndx
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig
    real(dp) :: LJ

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      rz = disp(iDisp)%z_new  -  self%zWallPos
      call curbox%Boundary(rz=rz)
      if(rz < 0.0) then
        accept = .false.
        return
      endif 
      ep = self % eps_star(atmType1)
      sig = self % sig_star(atmType1)
      LJ = (sig/rz)
      LJ = LJ*LJ*LJ
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_Diff = E_Diff + LJ

      rz = curbox % atoms(3, iAtom)  -  self%zWallPos
      call curbox%Boundary(rz=rz)
      LJ = (sig/rz)
      LJ = LJ*LJ*LJ
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_Diff = E_Diff - LJ
    enddo
 
  end subroutine
  !=====================================================================
  subroutine New_LJ_Wall(self, curbox, disp, tempList, tempNNei, E_Diff, accept)
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Addition), intent(in) :: disp(:)
    integer, intent(in) :: tempList(:,:), tempNNei(:)
    real(dp), intent(inOut) :: E_Diff
    logical, intent(out) :: accept

    integer :: iDisp, iAtom, jAtom, dispLen, maxNei, listIndx, jNei
    integer :: atmType1, atmType2
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig
    real(dp) :: LJ
    real(dp) :: rmin_ij      

    dispLen = size(disp)
    E_Diff = 0E0_dp
    accept = .true.

!    write(*,*) "Length:", size(tempNNei)
!    write(*,*) "Length:", size(tempList)
!    write(*,*)
    do iDisp = 1, dispLen
      iAtom = disp(iDisp)%atmIndx
      atmType1 = curbox % AtomType(iAtom)
      rz = disp(iDisp)%z_new - self%zWallPos
      call curbox%Boundary(rz=rz)
      if(rz < 0.0) then
        accept = .false.
        return
      endif
      ep = self % eps_star(atmType1)
      sig = self % sig_star(atmType1)
      LJ = (sig/rz)
      LJ = LJ*LJ*LJ
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_Diff = E_Diff + LJ
    enddo
  end subroutine
  !=====================================================================
  subroutine Old_LJ_Wall(self, curbox, disp, E_Diff)
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    class(SimBox), intent(inout) :: curbox
    type(Deletion), intent(in) :: disp(:)
    real(dp), intent(inOut) :: E_Diff
    integer :: iDisp, iAtom, jAtom, remLen, jNei
    integer :: atmType1, atmType2
    integer :: molEnd, molStart
    real(dp) :: rx, ry, rz, rsq
    real(dp) :: ep, sig
    real(dp) :: LJ, LJ2
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp

!    globIndx = curBox % MolGlobalIndx(disp(1)%molType, disp(1)%molIndx)
    call curBox % GetMolData(disp(1)%molIndx, molEnd=molEnd, molStart=molStart)

    do iAtom = molStart, molEnd
      atmType1 = curbox % AtomType(iAtom) 
      rz = curbox % atoms(3, iAtom) - self%zWallPos
      call curbox%Boundary(rz=rz)
      ep = self % eps_star(atmType1)
      sig = self % sig_star(atmType1)
      LJ = (sig/rz)
      LJ = LJ*LJ*LJ
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_Diff = E_Diff - LJ
    enddo
  end subroutine
  !=====================================================================
  subroutine OrthoVol_LJ_Wall(self, curbox, disp, E_Diff, accept)
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
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
    real(dp) :: ep, sig, ratio
    real(dp) :: LJ, LJ2
    real(dp) :: rmin_ij      

    E_Diff = 0E0_dp
    do iAtom = 1, curbox%nMaxAtoms
      if( curbox%IsActive(iAtom) ) then
        cycle
      endif
      atmType1 = curbox % AtomType(iAtom)
      molIndx1 = curbox % MolIndx(iAtom)
      dzi = curbox % centerMass(3, molIndx1) * (disp(1)%zScale-1E0_dp)
      rz =  curbox % atoms(3, iAtom) + dzi - self%zWallPos
      call curbox%BoundaryNew(rz=rz, disp=disp)
      if(rz < 0.0) then
        accept = .false.
        return
      endif
      ep = self % eps_star(atmType1)
      sig = self % sig_star(atmType1)
      LJ = (sig/rz)
      LJ = LJ*LJ*LJ
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_Diff = E_Diff + LJ

      rz = curbox % atoms(3, iAtom) - self%zWallPos
      call curbox%Boundary(rz=rz)
      LJ = (sig/rz)
      LJ = LJ*LJ*LJ
      LJ = ep * LJ * (coeff1*LJ*LJ - 1E0_dp)
      E_Diff = E_Diff - LJ
    enddo


  end subroutine
  !=====================================================================
  subroutine ProcessIO_LJ_Wall(self, line)
    use Common_MolInfo, only: nAtomTypes
    use Input_Format, only: CountCommands, GetXCommand
    use Input_Format, only: maxLineLen
    use Units, only: inEngUnit, inLenUnit
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    character(len=maxLineLen), intent(in) :: line
    character(len=30) :: command
    logical :: param = .false.
    integer :: jType, lineStat
    integer :: type1, type2, nPar
    real(dp) :: ep, sig, zWall, rMin
  

    call GetXCommand(line, command, 1, lineStat)

    select case(trim(adjustl(command)))
      case("wallposition")
        call GetXCommand(line, command, 2, lineStat)
        read(command, *) zWall
        self % zWallPos = zWall

      case default
        param = .true.
    end select


    if(param) then
!      call GetAllCommands(line, parlist, nPar, lineStat)
      call CountCommands(line, nPar)
      select case(nPar)
        case(4)
          read(line, *) type1, ep, sig
          ep = ep * inEngUnit
          sig = sig * inLenUnit

          self%eps_star(type1) = ep 
          self%sig_star(type1) = sig 
        case default
          lineStat = -1
      end select
    endif

  end subroutine
  !=============================================================================+
  function GetCutOff_LJ_Wall(self) result(rCut)
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    real(dp) :: rCut

    rCut = 0.0E0_dp
  end function
  !=====================================================================
  subroutine Prologue_LJ_Wall(self)
    use Common_MolInfo, only: nAtomTypes
    use ParallelVar, only: nout
    implicit none
    class(Pair_LJ_Wall), intent(inout) :: self
    integer :: i, j


  end subroutine
  !=====================================================================
end module
!=====================================================================
